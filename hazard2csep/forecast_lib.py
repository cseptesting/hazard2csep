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


def get_rate_simple_fault(sources, max_depth=200, buffer=0, *args, **kwargs):
    """
    Get the rates and polygons from simple fault sources.
    :param sources: list of SimpleFaultSources
    :param buffer_shp: pending
    :return:
    """
    polygons = []
    rates = []
    bin_width = []
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
            src.lower_seismogenic_depth,
            src.dip,
            # minimum resolution of ~ dh/2
            max(mesh_spacing, 5))
        area_factor = 1

        if src.lower_seismogenic_depth >= max_depth:
            area_0 = surf.get_area()
            surf = SimpleFaultSurface.from_fault_data(
                src.fault_trace,
                src.upper_seismogenic_depth,
                # cut by max depth if fault extends
                max_depth,
                src.dip,
                max(mesh_spacing, 5))
            area_factor = surf.get_area() / area_0

        rate = [(i, j * area_factor) for (i, j)
                in src.get_annual_occurrence_rates()]

        rates.append(rate)

        if not hasattr(src, 'buffer_polygon'):
            polygons.append(shapely.geometry.Polygon(
                    [(i, j) for i, j in zip(*surf.get_surface_boundaries())]))
        else:
            polygons.append(src.buffer_polygon)
        bin_width.append(src.mfd.bin_width)
    if isinstance(buffer, (int, float)) and buffer != 0:
        polygons = [i.buffer(buffer) for i in polygons]

    return polygons, rates, bin_width


def get_rate_complex_fault(sources, max_depth=200, *args, **kwargs):

    polygons = []
    rates = []
    bin_width = []
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
        bin_width.append(src.mfd.bin_width)
    return polygons, rates, bin_width


def get_rate_area_source(sources, max_depth=200, *args, **kwargs):
    polygons = []
    rates = []
    bin_width = []
    for src in sources:
        if max_depth < src.upper_seismogenic_depth:
            continue
        elif (src.upper_seismogenic_depth < max_depth <
              src.lower_seismogenic_depth):
            factor = (max_depth - src.upper_seismogenic_depth) / (
                 src.lower_seismogenic_depth - src.upper_seismogenic_depth)
        else:
            factor = 1

        poly = Polygon([(i, j) for i, j in zip(src.polygon.lons,
                                               src.polygon.lats)])
        rates.append([(i, j * factor) for (i, j)
                      in src.get_annual_occurrence_rates()])
        polygons.append(poly)
        bin_width.append(src.mfd.bin_width)
    return polygons, rates, bin_width


def get_rate_point_source(sources, max_depth=200, *args, **kwargs):

    points = []
    rates = []
    bin_width = []

    for src in sources:
        if max_depth < src.upper_seismogenic_depth:
            continue
        elif (src.upper_seismogenic_depth < max_depth <
              src.lower_seismogenic_depth):
            factor = (
             src.lower_seismogenic_depth - src.upper_seismogenic_depth) / (
             max_depth - src.upper_seismogenic_depth)
        else:
            factor = 1

        point = Point(src.location.longitude, src.location.latitude)
        rates.append([(i, j * factor) for (i, j)
                      in src.get_annual_occurrence_rates()])
        points.append(point)
        bin_width.append(src.mfd.bin_width)

    return points, rates, bin_width


def get_rate_mpoint_source(sources, *args, **kwargs):

    points = []
    rates = []
    bin_width = []
    for src in sources:
        for point_src in src:
            point = Point(point_src.location.longitude,
                          point_src.location.latitude)
            rates.append(point_src.get_annual_occurrence_rates())
            points.append(point)
            bin_width.append(point_src.mfd.bin_width)
    return points, rates, bin_width


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


# def project_mfd(rates, magnitudes):
#
#     mags = np.array([i[0] for i in rates]).round(2)
#     rate = np.array([i[1] for i in rates])
#
#     idxs = np.digitize(mags, magnitudes - 0.00001)
#     grid_rates = np.zeros(len(magnitudes))
#     for i, r in zip(idxs, rate):
#         if i > 0:
#             grid_rates[i-1] += r
#
#     return grid_rates
# import numpy as np

def project_mfd(rates, magnitudes, bin_width, method="center"):
    """Rebin rates from source bin *centers* onto target bins given by *left edges*.

    Args:
        rates (array-like of shape (N, 2)):
            Sequence of (center, value) pairs. The value is the per-bin total.
            Source centers are assumed to be uniformly spaced
        magnitudes (array-like of shape (M,)):
            Left edges of the *target* bins. Target bins are assumed uniformly spaced.
        method (str):
            - "center": assign each source bin's total to the single target bin that
              contains its center .
            - "overlap": distribute each source bin's total proportionally across all
              target bins it overlaps, preserving totals.

    Returns:
        np.ndarray of shape (M,): totals per target bin.

    Raises:
        ValueError: if spacing is not uniform when method="overlap", or inputs are
                    insufficient to infer widths.
    """
    rates = np.asarray(rates, dtype=float)
    if rates.ndim != 2 or rates.shape[1] != 2:
        raise ValueError("rates must be an array of shape (N, 2) with (center, value).")
    centers = rates[:, 0]
    values  = rates[:, 1]

    mags = np.asarray(magnitudes, dtype=float)
    if mags.ndim != 1 or mags.size < 1:
        raise ValueError("magnitudes must be a 1-D array of left edges.")

    # --- Target grid spacing (required for both methods) ---
    if mags.size < 2:
        if method == "overlap":
            raise ValueError("Need at least two target left edges")
        out = np.zeros(1, dtype=float)
        out[0] = values.sum()
        return out

    dt_all = np.diff(mags)
    dt = dt_all.mean()
    if not np.allclose(dt_all, dt, atol=1e-12, rtol=1e-9):
        raise ValueError("Target left edges must be uniformly spaced.")
    tL = mags
    tR = mags + dt
    M = mags.size
    if method.lower() == "center":

        idxs = np.digitize(centers, tL, right=False) - 1

        valid = (centers >= tL[0]) & (centers < tR[-1])

        centers_valid = centers[valid]
        values_valid = values[valid]
        idxs_valid = idxs[valid]

        out = np.zeros(M, dtype=float)
        for i, v in zip(idxs_valid, values_valid):
            if 0 <= i < M:
                out[i] += v
        return out


    elif method.lower() == "overlap":

        if centers.size < 2:
            ds_all = np.array([bin_width])
            # raise ValueError("Need at least two centers to infer source bin width.")
        else:
            ds_all = np.diff(centers)
        ds = ds_all.mean()
        if not np.allclose(ds_all, ds, atol=1e-12, rtol=1e-9):
            raise ValueError("Source centers must be uniformly spaced ")

        sL = centers - ds / 2.0
        sR = centers + ds / 2.0

        out = np.zeros(M, dtype=float)

        for left_j, right_j, vj in zip(sL, sR, values):
            wj = right_j - left_j
            if wj <= 0 or vj == 0:
                continue
            i_start = np.searchsorted(tR, left_j, side="right")
            i_end = np.searchsorted(tL, right_j, side="left") - 1
            if i_start > i_end or i_end < 0 or i_start >= M:
                continue
            i_start = max(i_start, 0)
            i_end = min(i_end, M - 1)
            for i in range(i_start, i_end + 1):
                oL = max(left_j, tL[i])
                oR = min(right_j, tR[i])
                overlap = oR - oL
                if overlap > 0:
                    out[i] += vj * (overlap / wj)

        return out

    else:
        raise ValueError("method must be 'center' or 'overlap'.")


def return_rates(sources, region=None, min_mag=4.7, max_mag=8.1, dm=0.2,
                 max_depth=200, dh=0.1, buffer=0, return_diag=False,):
    diagnostics = {"groups": []}

    # ... inside each source-group block, after you have
    # rebinned_total_rate for that group, do something like:

    # return clean magnitude bin array, defined by their left borders.
    csep_magnitudes = sm_lib.cleaner_range(min_mag, max_mag, dm).round(1)

    # Get source by type
    point_srcs = [src for src in sources if isinstance(src, PointSource)]
    mpoint_srcs = [src for src in sources if isinstance(src, MultiPointSource)]
    area_srcs = [src for src in sources if isinstance(src, AreaSource)]
    sf_srcs = [src for src in sources if isinstance(src, SimpleFaultSource)]
    cf_srcs = [src for src in sources if isinstance(src, ComplexFaultSource)]

    # Get faults' buffer from shapefile
    if isinstance(buffer, str):
        sf_srcs = sm_lib.parse_fault_source_buffer(sf_srcs, buffer_shp=buffer)

    # Load region
    if region:
        csep_grid = [shapely.geometry.Polygon(
            [i.points[0], i.points[3], i.points[2], i.points[1]])
            for i in region.polygons]
    # Make region if not given
    else:
        _, region = region_lib.make_region(sources, fill=False, dh=dh,
                                           fault_buffer=buffer)
        csep_grid = [shapely.geometry.Polygon(
            [i.points[0], i.points[3], i.points[2], i.points[1]])
            for i in region.polygons]
    print("\n[return_rates] region built")
    print("  n_cells:", len(csep_grid))
    print("  csep_mags:", csep_magnitudes)

    # Initialize forecast data
    csep_gdf = gpd.GeoDataFrame({'geometry': csep_grid})
    forecast_data = np.zeros((len(csep_grid), len(csep_magnitudes)))

    # Pair source type with functions
    src2func_map = ((sf_srcs, get_rate_simple_fault),
                    (cf_srcs, get_rate_complex_fault),
                    (area_srcs, get_rate_area_source),
                    (point_srcs, get_rate_point_source),
                    (mpoint_srcs, get_rate_mpoint_source))

    log.info(f'Processing {len(sources)} sources')
    for src_grp, func in src2func_map:
        # Get rates and polygons per src type
        polygons, rates, bin_widths = func(src_grp, max_depth=max_depth, buffer=buffer)
        if not src_grp:
            continue

        print("\n[return_rates] source group:", src_grp[0].__class__.__name__)
        print("  n_sources:", len(src_grp))
        print("  n_polygons:", len(polygons))
        print("  example_bin_width:", bin_widths[0] if bin_widths else None)

        # quick sanity of total original rates in this group
        total_src_rate = sum(
            sum(v for (_, v) in r) for r in rates
        )
        print("  total_src_rate (sum of all MFDs in group):", total_src_rate)
        if len(src_grp) > 0:
            log.info(f'Intersecting {len(polygons)}'
                     f' {src_grp[0].__class__.__name__} with CSEP grid')
        poly2csep = intersect_geom2grid(polygons, csep_gdf)


        # Allocate rates to grid cells
        for i, ((indices, frac), rate, bin_width) in enumerate(
                zip(poly2csep, rates, bin_widths)
        ):
            try:
                if indices.size > 0:
                    centers = np.array([m for (m, _) in rate])
                    vals = np.array([v for (_, v) in rate])

                    print(f"\n  [src #{i}]")
                    print("    indices:", indices)
                     # ("    frac:", frac)

                    try:
                        print("    frac_sum:", frac.sum())
                    except Exception as e:
                        print("    frac_sum couldnt be calculated, because its an int (we go point by point")
                    print("    src_centers(min,max):", centers.min(), centers.max())
                    print("    src_total_rate:", vals.sum())
                    print("    bin_width:", bin_width)

                    rate_mags = project_mfd(rate, csep_magnitudes, bin_width)
                    print("    rebinned_total_rate:", rate_mags.sum())
                    diagnostics["groups"].append({
                        "type": "MultiPointSource",  # or whatever group type
                        "rebinned_total_rate": rate_mags.sum(),
                    })
                    contrib = (frac * rate_mags).squeeze()
                    print("    contrib.shape:", contrib.shape)
                    print("    contrib_total:", contrib.sum())

                    forecast_data[indices] += contrib
            except ValueError:
                print("  ValueError in src #", i, "indices:", indices)
                raise
    a = GriddedForecast(data=forecast_data, region=region,
                        magnitudes=csep_magnitudes)
    print("\n[return_rates] final checks")
    print("  forecast_data.shape:", forecast_data.shape)

    print("  forecast_total_rate:", forecast_data.sum())

    if return_diag:
        return a, diagnostics
    else:
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


# def write_raster(forecast, dest='forecast.tiff', log=True,
#                  mbin=False, depths=(0, 30), dh=0.1):
#     data = np.log(forecast.data) if log else forecast.data
#     magnitudes = forecast.magnitudes
#     points = forecast.region.origins()
#
#     cells = np.vstack((points[:, 0], points[:, 0] + dh, points[:, 1],
#                        points[:, 1] + dh)).T
#     depths = np.vstack((depths[0]*np.ones(points.shape[0]),
#                         depths[1]*np.ones(points.shape[0]))).T


if __name__ == '__main__':

    pass

