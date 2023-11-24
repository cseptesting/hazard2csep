import os
import geopandas as gpd
import pandas as pd
from scipy.stats import mode
from openquake.hazardlib.source.complex_fault import ComplexFaultSource
from openquake.hazardlib.source.simple_fault import SimpleFaultSource
from openquake.hazardlib.source.area import AreaSource
from openquake.hazardlib.source.point import PointSource
from openquake.hazardlib.source.multi_point import MultiPointSource
import numpy as np
import matplotlib.pyplot as plt
import shapely.geometry
import csep
import cartopy
from csep.core.regions import CartesianGrid2D
import logging
from multiprocessing import Pool
from hazard2csep.sm_lib import cleaner_range

log = logging.getLogger('hazard2csepLogger')


def fill_holes(coords,
               # final_region=None,
               dh=0.1):

    log.info(f'Filling holes in the region')
    bounds = np.array([180, 90, -180, -90])

    bounds = np.array([
        *np.min([bounds[:2], coords.min(axis=0)], axis=0),
        *np.max([bounds[2:], coords.max(axis=0)], axis=0)])
    x = cleaner_range(np.round(bounds[0], 1) - dh,
                      np.round(bounds[2], 1) + dh, h=dh)
    y = cleaner_range(np.round(bounds[1], 1) - dh,
                      np.round(bounds[3], 1) + dh, h=dh)

    initial_grid = np.vstack([i.ravel() for i in np.meshgrid(x, y)]).T
    initial_cells = [shapely.geometry.Polygon(
        [(i[0], i[1]), (i[0] + dh, i[1]), (i[0] + dh, i[1] + dh),
         (i[0], i[1] + dh)]) for i in initial_grid]
    initial_region = gpd.GeoSeries(initial_cells)

    # if final_region is None:
    final_cells = [shapely.geometry.Polygon(
        [(i[0], i[1]), (i[0] + dh, i[1]), (i[0] + dh, i[1] + dh),
         (i[0], i[1] + dh)]) for i in coords]
    final_region = gpd.GeoSeries(final_cells)

    tag_from_holes = np.zeros(len(initial_cells))

    # Buffing the region by 1/10 of the cell size
    buff_geodf = gpd.GeoDataFrame(geometry=final_region).buffer(dh/10.)
    # Dissolving the buffered region
    diss_geodf = gpd.GeoDataFrame(geometry=buff_geodf).dissolve()

    # Getting any hole in the region

    inner_polygons = []
    for elem in diss_geodf.geometry:

        int_rings = []
        if isinstance(elem, shapely.geometry.MultiPolygon):
            for poly in elem.geoms:
                int_rings.append(poly.interiors)
        else:
            int_rings.append(elem.interiors)

        for ring in int_rings:
            if ring:
                for j in ring:
                    xy = np.array(j.coords.xy).T
                    poly = shapely.geometry.Polygon(xy)
                    inner_polygons.append(poly)

    for poly in inner_polygons:
        tag_from_holes = np.logical_or(tag_from_holes,
                                       initial_region.intersects(poly))

    if tag_from_holes.any():
        hole_region = initial_region[tag_from_holes]
        hole_coords = np.array([i.exterior.coords[0]
                                for i in hole_region.geometry])
        coords = np.unique(np.vstack((coords, hole_coords)), axis=0)
        coords = coords[np.argsort(coords[:, 1])]

    return coords


def make_region(sources, dh=0.1, fill=False):
    """
    Creates a CSEP region from a Source Model. A Lat/Lon uniform grid is
     created from the boundaries of the SM. Active cells are determined if
     a source touches the cell.
     The final region is conformed from the Union of all sub-models.

    :param sources: list of sources
    :param dh: cell size
    :return:
    """
    log.info('Computing region')
    # Get sources by type
    point_srcs = [src for src in sources if isinstance(src, PointSource)]
    mpoint_srcs = [src for src in sources if isinstance(src, MultiPointSource)]
    area_srcs = [src for src in sources if isinstance(src, AreaSource)]
    sf_srcs = [src for src in sources if isinstance(src, SimpleFaultSource)]
    cf_srcs = [src for src in sources if isinstance(src, ComplexFaultSource)]

    # Set initial bounds
    bounds = np.array([180, 90, -180, -90])

    # Get a point array from PointSources and MultiPointSources
    pointsrc_coords = np.zeros((0, 2))
    if point_srcs:
        pointsrc_coords = np.array([(src.location.longitude,
                                     src.location.latitude)
                                    for src in point_srcs])
    if mpoint_srcs:
        pointsrc_coords = np.vstack((pointsrc_coords,
                                    *[src.mesh.array.T for src in sources if
                                      isinstance(src, MultiPointSource)]))

    # Gets the bounds (min lon, min lat, max lon, max lat) from points
    if pointsrc_coords.size > 0:
        bounds = np.array([
            *np.min([bounds[:2], pointsrc_coords.min(axis=0)], axis=0),
            *np.max([bounds[2:], pointsrc_coords.max(axis=0)], axis=0)])

    # Gets the bounds (min lon, min lat, max lon, max lat) from polygons
    for src in [*area_srcs, *sf_srcs, *cf_srcs]:
        bounds = np.array([
            *np.min([bounds[:2], src.polygon.get_bbox()[:2]], axis=0),
            *np.max([bounds[2:], src.polygon.get_bbox()[2:]], axis=0)])

    log.info('Geographic bounds [lon_min, lat_min, lon_max, lat_max]:')
    log.info(f'\t> {bounds}')

    # Create CSEP Grid from the Source Model's bounds and grid size
    x = cleaner_range(np.round(bounds[0], 1) - dh,
                      np.round(bounds[2], 1) + dh, h=dh)
    y = cleaner_range(np.round(bounds[1], 1) - dh,
                      np.round(bounds[3], 1) + dh, h=dh)
    grid = np.vstack([i.ravel() for i in np.meshgrid(x, y)]).T
    initial_cells = [shapely.geometry.Polygon(
        [(i[0], i[1]), (i[0] + dh, i[1]), (i[0] + dh, i[1] + dh),
         (i[0], i[1] + dh)]) for i in grid]
    initial_region = gpd.GeoSeries(initial_cells)

    # Initialize array of active CSEP cells
    tag_cells = np.zeros(len(grid))

    # Check which cells are touched by Polygon-type Sources
    polygons = [shapely.geometry.Polygon(
        [(i, j) for i, j in zip(src.polygon.lons, src.polygon.lats)])
        for src in [*area_srcs, *sf_srcs, *cf_srcs]]

    if len(polygons) > 0:
        log.info(f'Intersecting region polygons with CSEP region')

    for poly in polygons:
        tag_cells = np.logical_or(tag_cells, initial_region.intersects(poly))

    # Check which cells are touched by Point-type Sources
    if pointsrc_coords.size > 0:
        log.info(f'Intersecting region points with CSEP region')
        df = pd.DataFrame(
            {'lon': pointsrc_coords[:, 0], 'lat': pointsrc_coords[:, 1]})
        df['coords'] = list(zip(df['lon'], df['lat']))
        df['coords'] = df['coords'].apply(shapely.geometry.Point)
        points = gpd.GeoDataFrame(df, geometry='coords')
        cells = gpd.GeoDataFrame(geometry=initial_region)
        points_in_reg = gpd.tools.sjoin(points, cells, predicate='within',
                                        how='left')
        idxs = np.unique(points_in_reg.index_right.to_numpy())

        tag_from_points = np.zeros(len(grid))
        tag_from_points[idxs] = 1
        tag_cells = np.logical_or(tag_cells, tag_from_points)

    # Crops CSEP grid to those cells that were touched by sources
    if tag_cells.any():
        final_region = initial_region[tag_cells]
    else:
        final_region = initial_region

    coords = np.array([i.exterior.coords[0] for i in final_region.geometry])

    # Fill holes in the region
    if fill:
        coords = fill_holes(coords, dh)

    csep_region = CartesianGrid2D.from_origins(coords, dh=dh)

    return coords, csep_region


def intersect_region(region1, *args):
    """
    Intersects multiple regions to get the common cells between them.

    :param region1: (str, CartesianGrid2D) path to region file or
     CartesianGrid2D object
    :param args: (list, str, CartesianGrid2D) paths to region files or
     CartesianGrid2D objects
    """
    regions = args
    other_ = []

    if isinstance(region1, str):
        origins = np.loadtxt(region1)
        dh1 = mode(np.diff(np.unique(origins[:, 0]))).mode
        dh2 = mode(np.diff(np.unique(origins[:, 1]))).mode
        dh = np.min([dh1, dh2])

        region1 = CartesianGrid2D.from_origins(origins, dh=dh)

    for reg in regions:

        if isinstance(reg, str):
            origins = np.loadtxt(reg)
            dh1 = mode(np.diff(np.unique(origins[:, 0]))).mode
            dh2 = mode(np.diff(np.unique(origins[:, 1]))).mode
            dh = np.min([dh1, dh2])
            reg = CartesianGrid2D.from_origins(origins, dh=dh)
        other_.append(reg)

    reg1_poly = [shapely.geometry.Polygon(
        [i.points[0], i.points[3], i.points[2], i.points[1], i.points[0]])
                    for i in region1.polygons]
    reg1_gdf = gpd.GeoSeries(reg1_poly)
    tag_cells = np.ones(len(reg1_poly))
    for i, reg2 in enumerate(other_):

        reg2_points = [shapely.geometry.Point(i) for i in reg2.midpoints()]
        points = gpd.GeoDataFrame(geometry=gpd.GeoSeries(reg2_points))
        points_in_reg = gpd.tools.sjoin(points,
                                        gpd.GeoDataFrame(geometry=reg1_gdf),
                                        predicate='within',
                                        how='left')

        idxs = np.unique(points_in_reg.index_right.to_numpy())
        idxs = idxs[~np.isnan(idxs)].astype(int)

        tag_from_points = np.zeros(len(reg1_poly))
        tag_from_points[idxs] = 1
        tag_cells = np.logical_and(tag_cells, tag_from_points)

    return (region1.origins()[tag_cells],
            CartesianGrid2D.from_origins(region1.origins()[tag_cells]))


def plot_region(grid, fname):

    fig = plt.figure()
    ax = fig.add_subplot(111, projection=cartopy.crs.PlateCarree())
    coords = np.unique(grid[:, 0])
    if coords.size > 1:
        dh = np.diff(np.sort(coords)).max()
    else:
        dh = np.diff(np.sort(np.unique(grid[:, 1]))).max()

    extent = [grid[:, 0].min() - dh, grid[:, 0].max() + dh,
              grid[:, 1].min() - dh, grid[:, 1].max() + dh]
    ax = csep.utils.plots.plot_basemap(ax=ax,
                                       extent=extent,
                                       basemap='ESRI_terrain')
    ax.plot(*grid.T, '.', color='red')
    ax.gridlines(draw_labels=True)
    if fname:
        fig.savefig(fname)
    else:
        plt.show()



