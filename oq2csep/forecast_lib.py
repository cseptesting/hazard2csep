import shapely
from csep.core.regions import CartesianGrid2D
from csep.core.forecasts import GriddedForecast
from openquake.hazardlib.source.complex_fault import ComplexFaultSource
from openquake.hazardlib.source.simple_fault import SimpleFaultSource
from openquake.hazardlib.source.area import AreaSource
from openquake.hazardlib.source.point import PointSource
from openquake.hazardlib.source.multi_point import MultiPointSource
from openquake.hazardlib.geo.surface import (ComplexFaultSurface,
                                             SimpleFaultSurface)
from openquake.hazardlib.mfd import (TruncatedGRMFD, EvenlyDiscretizedMFD,
                                     TaperedGRMFD)

import numpy as np
from oq2csep import region_lib
import matplotlib.pyplot as plt
from os import path
from shapely.geometry import Polygon, Point, MultiPoint
import geopandas as gpd
import pandas as pd
from oq2csep import sm_lib

import logging
log = logging.getLogger('oq2csepLogger')

max_depth = 50


def get_rate_simple_fault(sources, buffer_shp=None):
    """
    Get the rates and polygons from simple fault sources.
    :param sources: list of SimpleFaultSources
    :param buffer_shp: pending
    :return:
    """
    polygons = []
    rates = []
    for src in sources:

        # If fault is entirely below max depth, skip
        if src.upper_seismogenic_depth >= max_depth:
            continue

        # Approximation to get the fault spacing from trace resolution
        mesh_spacing = np.min(np.sum(np.diff(src.fault_trace.coo[:, :2],
                                             axis=0)**2, axis=1)**0.5) * 110.

        surf = SimpleFaultSurface.from_fault_data(
            src.fault_trace,
            src.upper_seismogenic_depth,
            # cut by max depth if fault extends
            min(src.lower_seismogenic_depth, max_depth),
            src.dip,
            # minimum resolution of ~ dh/2
            max(mesh_spacing, 5))

        rates.append(src.get_annual_occurrence_rates())
        polygons.append(shapely.geometry.Polygon(
                [(i, j) for i, j in zip(*surf.get_surface_boundaries())]))

    return polygons, rates


def get_rate_complex_fault(sources):

    polygons = []
    rates = []
    for src in sources:

        # Select automating spacing according to fault edges resolution
        rupt_spacing = np.hstack([
            np.sum(np.diff(i.coo[:, :2],axis=0)**2, axis=1)**0.5 * 110.
            for i in src.edges]).min()


        surf = ComplexFaultSurface.from_fault_data(
            src.edges,
            # minimum resolution of ~ dh/2
            max(rupt_spacing, 5))

        rates.append(src.get_annual_occurrence_rates())
        polygons.append(shapely.geometry.Polygon(
            [(i, j) for i, j in zip(*surf.get_surface_boundaries())]))

    return polygons, rates


def get_rate_area_source(sources):

    polygons = []
    rates = []
    for src in sources:
        poly = Polygon([(i, j) for i, j in zip(src.polygon.lons,
                                               src.polygon.lats)])
        rates.append(src.get_annual_occurrence_rates())
        polygons.append(poly)

    return polygons, rates


def get_rate_point_source(sources):

    points = []
    rates = []
    for src in sources:
        point = Point(src.location.longitude, src.location.latitude)
        rates.append(src.get_annual_occurrence_rates())
        points.append(point)

    return points, rates


def get_rate_mpoint_source(sources):

    points = []
    rates = []
    for src in sources:
        for point_src in src:
            point = Point(point_src.location.longitude,
                          point_src.location.latitude)
            rates.append(point_src.get_annual_occurrence_rates())
            points.append(point)
    return points, rates


def intersect_geom2grid(geometries, csep_grid):

    geom2grid_map = []
    if all(isinstance(geom, Polygon) for geom in geometries):
        for geom in geometries:
            if isinstance(geom, Polygon):
                # Get index map from polygons to cells
                idxs = np.argwhere(csep_grid.intersects(geom))
                # Get fraction of area of each cell intersecting polygon
                frac = np.array([csep_grid.iloc[i].intersection(geom).area
                                 for i in idxs]) / geom.area
                geom2grid_map.append([idxs, frac])

    elif all(isinstance(geom, Point) for geom in geometries):

        points = gpd.GeoDataFrame(geometry=geometries)
        # Get index map from point_sources to grid cells
        points_in_reg = gpd.tools.sjoin(points, csep_grid,
                                        predicate='within',
                                        how='left')
        # Return index or empty array if point lies outside grid
        geom2grid_map = [[np.array(int(i)), 1] if not np.isnan(i)
                         else [np.array([]), np.array([])]
                         for i in points_in_reg.index_right]
    else:
        raise TypeError('Invalid or non-uniform geometry types')
    return geom2grid_map


def project_mfd(rates, magnitudes):

    mags = np.array([i[0] for i in rates]).round(2)
    rate = np.array([i[1] for i in rates])

    idxs = np.digitize(mags, magnitudes)
    grid_rates = np.zeros(len(magnitudes))
    for i, r in zip(idxs, rate):
        if i > 0:
            grid_rates[i-1] += r

    return grid_rates

def return_rates(sources, region=None, min_mag=4.7, max_mag=8.1, dm=0.2):

    magnitudes = sm_lib.cleaner_range(min_mag, max_mag, dm).round(1)


    # Get source by type
    point_srcs = [src for src in sources if isinstance(src, PointSource)]
    mpoint_srcs = [src for src in sources if isinstance(src, MultiPointSource)]
    area_srcs = [src for src in sources if isinstance(src, AreaSource)]
    sf_srcs = [src for src in sources if isinstance(src, SimpleFaultSource)]
    cf_srcs = [src for src in sources if isinstance(src, ComplexFaultSource)]

    # Load region
    if region:
        csep_grid = [shapely.geometry.Polygon(
            [i.points[0], i.points[3], i.points[2], i.points[1]])
            for i in region.polygons]
    # Make region if not given
    else:
        _, region = region_lib.parse_region(sources, fill=False)
        csep_grid = [shapely.geometry.Polygon(
            [i.points[0], i.points[3], i.points[2], i.points[1]])
            for i in region.polygons]

    # Initialize forecast data
    csep_gdf = gpd.GeoDataFrame({'geometry': csep_grid})
    forecast_data = np.zeros((len(csep_grid), len(magnitudes)))  # todo magnitudes bins

    # Pair source type with functions
    src2func_map = ((sf_srcs, get_rate_simple_fault),
                    (cf_srcs, get_rate_complex_fault),
                    (area_srcs, get_rate_area_source),
                    (point_srcs, get_rate_point_source),
                    (mpoint_srcs, get_rate_mpoint_source))

    log.info(f'Processing {len(sources)} sources')
    for src_grp, func in src2func_map:
        # Get rates and polygons per src type
        polygons, rates = func(src_grp)

        if len(src_grp) > 0:
            log.info(f'Intersecting {len(polygons)}'
                     f' {src_grp[0].__class__.__name__} with CSEP grid')
        poly2csep = intersect_geom2grid(polygons, csep_gdf)

        # Allocate rates to grid cells
        for i, ((indices, frac), rate) in enumerate(zip(poly2csep, rates)):
            if indices.size > 0:
                total_rate = np.sum([j for k, j in rate])  # todo

                print([i[0] for i in rate])
                total_rate = total_rate * frac
                rate_mags = project_mfd(rate, magnitudes)

                forecast_data[indices] += rate_mags

    a = GriddedForecast(data=forecast_data.reshape(-1, 1), region=region,
                        magnitudes=[5.])
    return a


if __name__ == '__main__':

    dir_module = path.join(path.dirname(__file__), '..')
    eshm20_sm = path.join(dir_module, 'eshm_test', 'eshm20', 'oq_computational',
                 'oq_configuration_eshm20_v12e_region_main/source_models')
    eshm20_subd_interface = [path.join(eshm20_sm, 'interface_v12b', i)
                             for i in ['CaA_IF2222222_M40.xml',
                                       'CyA_IF2222222_M40.xml',
                                       'GiA_IF2222222_M40.xml',
                                       'HeA_IF2222222_M40.xml']]
    fsbg = path.join(eshm20_sm, 'fsm_v09', 'fs_ver09e_model_aGR_SRL_ML_fMthr.xml')
    asm = path.join(eshm20_sm, 'asm_v12e', 'asm_ver12e_winGT_fs017_hi_abgrs_maxmag_upp.xml')
    ssm = path.join(eshm20_sm, 'ssm_v09',
               'seis_ver12b_fMthr_asm_ver12e_winGT_fs017_agbrs_point.xml')


    eshm13_sm = path.join(dir_module, 'eshm_test', 'eshm13',
                          'SHARE_OQ_input_20140807')
    seifa = path.join(eshm13_sm, 'seifa_model_test.xml')
    fsbg = path.join(eshm13_sm, 'faults_backg_source_model_test.xml')


    reg = path.join(dir_module, 'eshm_test', 'regions', 'region_final.txt')
    csep_reg = CartesianGrid2D.from_origins(np.loadtxt(reg))


    sm = sm_lib.parse_source_model([seifa])
    srcs = sm_lib.parse_srcs(sm)
    data = return_rates(srcs, region=csep_reg, min_mag=4.7, max_mag=8.1, dm=0.2)
    ax = data.plot(plot_args={'region_border': False})

    # ax.plot(*data.region.midpoints().T, '.', ms=10)
    plt.show()

    # if reg is given
    #     reg = read_region()
    # else:
    #     reg_get_region(sm)
    #
    # # parse sms
    # # parse srcs
    # # parse rates
    # # resample to mws
    # # project to mesh
    #
    #
    # mesh_fn = 'mesh_europe.csv'
    # mesh = np.genfromtxt(mesh_fn, delimiter=',')
    # mesh_shp = MultiPoint([Point(i[0],i[1]) for i in mesh])
    # rupt_spacing = 3
    # area_mesh = 8
    # max_depth = 35
    # reload = False
    # mw_bin = 0.1
    # for filename in files:
    #
    #     if reload is False:
    #
    #         parser = sourceconverter.SourceConverter(area_source_discretization=area_mesh,
    #                                                  width_of_mfd_bin=mw_bin)
    #         if isinstance(filename, list):
    #             source_list = filename
    #             filename = filename[-1]
    #             source_list.extend(subd)
    #             a = nrml.read_source_models(source_list, parser)
    #         else:
    #             source_list = [filename, 'asm_v12e/asm_ver12e_winGT_fs017_twingr.xml']
    #             source_list.extend(subd)
    #             a = nrml.read_source_models(source_list, parser)
    #         src_models = []
    #         for i in a:
    #             src_models.append(i)
    #         Rates = []
    #         Magnitudes = []
    #         Cells = []
    #         SrcType = []
    #         srcs = []
    #
    #         for smm in src_models:
    #             for sgc in smm:
    #                 for src in sgc:
    #                     srcs.append(src)
    #         for i, src in enumerate(srcs):
    #             if i % 50 == 0:
    #                 print('src', i)
    #             cc, mm, rr, st = return_rates(src, mesh_shp)
    #             if len(cc) > 0:
    #                 Cells.extend(cc)
    #                 Magnitudes.extend(mm)
    #                 Rates.extend(rr)
    #                 SrcType.extend(st)
    #             else:
    #                 print('ho')
    #         with open(path.splitext(filename)[0] + '.obj', 'wb') as pick:
    #             pickle.dump((Cells, Magnitudes, Rates, SrcType), pick)
    #     else:
    #         if isinstance(filename, list):
    #             filename = filename[-1]
    #         with open(path.splitext(filename)[0] + '.obj', 'rb') as pick:
    #             Cells, Magnitudes, Rates, SrcType = pickle.load( pick)
    #
    #     np.savetxt(path.splitext(filename)[0] + '_original.csv', np.array(Cells), delimiter=',')
    #
    #
    #     plt.scatter([i[0] for i in Cells], [i[1] for i in Cells],
    #                 c = np.log10([np.sum(i) for i in Rates]), vmin=-8, vmax=-3, s=0.01)
    #     plt.show()
    #
    #     data = project2mesh(mesh_fn, np.array(Cells), Magnitudes, Rates, SrcType)
    #     plt.scatter((data[:, 0]+data[:, 1])/2., (data[:, 2]+data[:, 3])/2.,
    #                 c=np.log10(data[:, 6:].sum(axis=1)), s=0.04)
    #     plt.show()
    #     mm = []
    #     for i in Magnitudes:
    #         mm.extend(i)
    #     min_mag = min(mm)
    #     max_mag = max(mm)
    #     dm = 0.2
    #     magnitudes = cleaner_range(min_mag, max_mag, dm).round(1)
    #
    #     if min(magnitudes) == 4.6:
    #         magnitudes += 0.1
    #     elif max(magnitudes) == 4.5:
    #         magnitudes += 0.2
    #     elif max(magnitudes) == 4.55:
    #         magnitudes += 0.15
    #
    #     header = 'lon_min,lon_max,lat_min,lat_max,depth_min,depth_max,' + ','.join([f'{i :.1f}' for i in magnitudes[:data.shape[1] -6]])
    #     fn_csv = path.splitext(filename)[0] + '.csv'
    #     np.savetxt(f'{fn_csv}', data, delimiter=',', header=header)

