import csep.models
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
from openquake.hazardlib.geo.mesh import RectangularMesh


from multiprocessing import Pool
import numpy as np
from hazard2csep import region_lib
from shapely.geometry import Polygon, Point
import geopandas as gpd
from hazard2csep import sm_lib
import pandas

import logging
log = logging.getLogger('hazard2csepLogger')



def get_rate_simple_fault(sources, max_depth=200, *args, **kwargs):
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

        if src.lower_seismogenic_depth <= max_depth:
            surf = SimpleFaultSurface.from_fault_data(
                src.fault_trace,
                src.upper_seismogenic_depth,
                src.lower_seismogenic_depth,
                src.dip,
                # minimum resolution of ~ dh/2
                max(mesh_spacing, 5))
            rate = src.get_annual_occurrence_rates()

        else:
            # Surface needs to be cropped and rates scaled
            surf_0 = SimpleFaultSurface.from_fault_data(
                src.fault_trace,
                src.upper_seismogenic_depth,
                src.lower_seismogenic_depth,
                src.dip,
                max(mesh_spacing, 5))

            surf = SimpleFaultSurface.from_fault_data(
                src.fault_trace,
                src.upper_seismogenic_depth,
                # cut by max depth if fault extends
                max_depth,
                src.dip,
                max(mesh_spacing, 5))
            area_factor = surf.get_area() / surf_0.get_area()

            rate = [(i, j * area_factor) for (i, j)
                    in src.get_annual_occurrence_rates()]

        rates.append(rate)
        polygons.append(shapely.geometry.Polygon(
                [(i, j) for i, j in zip(*surf.get_surface_boundaries())]))

    return polygons, rates


def get_rate_complex_fault(sources, max_depth=200, *args, **kwargs):

    polygons = []
    rates = []
    for src in sources:

        # Select automating spacing according to fault edges resolution
        # rupt_spacing = np.hstack([
        #     np.sum(np.diff(i.coo[:, :2],axis=0)**2, axis=1)**0.5 * 110.
        #     for i in src.edges]).min()
        rupt_spacing = 2.

        surf_0 = ComplexFaultSurface.from_fault_data(
            src.edges,
            # minimum resolution of ~ dh/2
            max(rupt_spacing, 2))

        depths = np.array([i[0] for i in surf_0.mesh.depths])

        depth_idx = np.argwhere(depths <= max_depth).squeeze()
        rec_mesh = RectangularMesh(surf_0.mesh.lons[depth_idx],
                                   surf_0.mesh.lats[depth_idx],
                                   surf_0.mesh.depths[depth_idx])
        surf = ComplexFaultSurface(rec_mesh)
        area_factor = surf.get_area()/surf_0.get_area()
        rates.append([(i, j * area_factor) for (i, j)
                      in src.get_annual_occurrence_rates()])
        polygons.append(shapely.geometry.Polygon(
            [(i, j) for i, j in zip(*surf.get_surface_boundaries())]))

    return polygons, rates


def get_rate_area_source(sources, *args, **kwargs):

    polygons = []
    rates = []
    for src in sources:
        poly = Polygon([(i, j) for i, j in zip(src.polygon.lons,
                                               src.polygon.lats)])
        rates.append(src.get_annual_occurrence_rates())
        polygons.append(poly)

    return polygons, rates


def get_rate_point_source(sources, *args, **kwargs):

    points = []
    rates = []
    for src in sources:
        point = Point(src.location.longitude, src.location.latitude)
        rates.append(src.get_annual_occurrence_rates())
        points.append(point)

    return points, rates


def get_rate_mpoint_source(sources, *args, **kwargs):

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
                geom2grid_map.append([idxs.squeeze(), frac])

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


def return_rates(sources, region=None, min_mag=4.7, max_mag=8.1, dm=0.2,
                 max_depth=200, dh=0.1):

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
        _, region = region_lib.make_region(sources, fill=False, dh=dh)
        csep_grid = [shapely.geometry.Polygon(
            [i.points[0], i.points[3], i.points[2], i.points[1]])
            for i in region.polygons]

    # Initialize forecast data
    csep_gdf = gpd.GeoDataFrame({'geometry': csep_grid})
    forecast_data = np.zeros((len(csep_grid), len(magnitudes)))

    # Pair source type with functions
    src2func_map = ((sf_srcs, get_rate_simple_fault),
                    (cf_srcs, get_rate_complex_fault),
                    (area_srcs, get_rate_area_source),
                    (point_srcs, get_rate_point_source),
                    (mpoint_srcs, get_rate_mpoint_source))
    print(csep_gdf)

    log.info(f'Processing {len(sources)} sources')
    for src_grp, func in src2func_map:
        # Get rates and polygons per src type
        polygons, rates = func(src_grp, max_depth=max_depth)

        if len(src_grp) > 0:
            log.info(f'Intersecting {len(polygons)}'
                     f' {src_grp[0].__class__.__name__} with CSEP grid')
        poly2csep = intersect_geom2grid(polygons, csep_gdf)

        # Allocate rates to grid cells
        for i, ((indices, frac), rate) in enumerate(zip(poly2csep, rates)):
            if indices.size > 0:
                rate_mags = project_mfd(rate, magnitudes)
                forecast_data[indices] += (frac * rate_mags).squeeze()

    a = GriddedForecast(data=forecast_data, region=region,
                        magnitudes=magnitudes)
    return a


def read_forecast(filename):

    def is_mag(num):
        try:
            m = float(num)
            if -1 < m < 12.:
                return True
            else:
                return False
        except ValueError:
            return False

    with open(filename, 'r') as file_:
        line = file_.readline()
        line = file_.readline()
        print(line.split(','))
        if len(line.split(',')) > 3:
            sep = ','
        else:
            sep = ' '

    data = pandas.read_csv(filename, sep=sep,
                           )
    data.columns = [i.strip() for i in data.columns]
    magnitudes = np.array([float(i) for i in data.columns if is_mag(i)])
    rates = data[[i for i in data.columns if is_mag(i)]].to_numpy()

    all_polys = data[
        ['lon_min', 'lon_max', 'lat_min', 'lat_max']].to_numpy()
    bboxes = [((i[0], i[2]), (i[0], i[3]), (i[1], i[3]), (i[1], i[2]))
              for i in all_polys]
    dh = float(all_polys[0, 3] - all_polys[0, 2])

    try:
        poly_mask = data['mask']
    except KeyError:
        poly_mask = np.ones(len(bboxes))

    region = CartesianGrid2D(
        [csep.models.Polygon(bbox) for bbox in bboxes], dh, mask=poly_mask)

    return GriddedForecast(data=rates, region=region, magnitudes=magnitudes)


def write_forecast(forecast, dest='forecast.csv',
                   fmt="csep", depths=(0, 30), dh=0.1):

    log.info(f'Writing forecast to {dest}')

    data = forecast.data
    magnitudes = forecast.magnitudes
    points = forecast.region.origins()

    cells = np.vstack((points[:, 0], points[:, 0] + dh, points[:, 1],
                       points[:, 1] + dh)).T
    depths = np.vstack((depths[0]*np.ones(points.shape[0]),
                        depths[1]*np.ones(points.shape[0]))).T

    if fmt == "csep":
        header = (
                'lon_min, lon_max, lat_min, lat_max, depth_min, depth_max, ' +
                ','.join([str(i) for i in magnitudes]))
        np.savetxt(dest, np.hstack((cells, depths, data)),
                   fmt=6*['%.1f'] + len(magnitudes)*['%.16e'],
                   header=header, comments='', delimiter=',')
    log.info(' > Total events: %.4f' % np.sum(data))


if __name__ == '__main__':

    pass

