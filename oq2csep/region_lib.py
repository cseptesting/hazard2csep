import geopandas as gpd
import pandas as pd
from openquake.hazardlib import nrml, sourceconverter
from openquake.hazardlib.source.complex_fault import ComplexFaultSource
from openquake.hazardlib.source.simple_fault import SimpleFaultSource
from openquake.hazardlib.source.area import AreaSource
from openquake.hazardlib.source.point import PointSource
from openquake.hazardlib.source.multi_point import MultiPointSource
from openquake.hazardlib.geo.surface import ComplexFaultSurface
from openquake.hazardlib.geo.line import Line
import numpy as np
import matplotlib.pyplot as plt
import shapely.geometry
import csep
import xml.etree.ElementTree as etree
import cartopy
from csep.core.regions import CartesianGrid2D
import logging


log = logging.getLogger('oq2csepLogger')


def cleaner_range(start, end, h):
    """ Returns array holding bin edges that doesn't contain floating point
    wander.

    Floating point wander can occur when repeatedly adding floating point
    numbers together. The errors propagate and become worse over the sum.
    This function generates the
    values on an integer grid and converts back to floating point numbers
     through multiplication.

     Args:
        start (float)
        end (float)
        h (float): magnitude spacing

    Returns:
        bin_edges (numpy.ndarray)
    """
    # convert to integers to prevent accumulating floating point errors
    const = 100000
    start = np.floor(const * start)
    end = np.floor(const * end)
    d = const * h
    return np.arange(start, end + d / 2, d) / const


def parse_source_model(fname):
    """
    Parse the OpenQuake Source Model file(s).
    Can be a single file pointint to a source model, or to a compound model
    consisting of two different files (e.g. point surfaces + subduction faults)

    In case fault conventions were deprecated, attemps to fix them.

    :param fname: path or list of paths (list, string)
    :return: returns the Openquake source model object
    """
    parser = sourceconverter.SourceConverter(
        area_source_discretization=10,
        width_of_mfd_bin=0.1
    )

    if not isinstance(fname, list):
        fname = [fname]
    try:
        src_model_nrml = nrml.read_source_models(fname, parser)
        src_models = list(src_model_nrml)

    except ValueError:
        src_models = []
        for file_ in fname:
            log.info('Fixing deprecated fault geometries')
            tree = etree.parse(file_)
            root = tree.getroot()
            xmlns = '{http://openquake.org/xmlns/nrml/0.4}'
            xmlns_gml = "{http://www.opengis.net/gml}"
            complexfaultsrcs = root[0].findall(f'{xmlns}complexFaultSource')
            for src in complexfaultsrcs:
                geom = src.find(f'{xmlns}complexFaultGeometry')
                top_edge = geom.find(f'{xmlns}faultTopEdge')
                top_ls = top_edge.find(f'{xmlns_gml}LineString')
                top_pos = top_ls.find(f'{xmlns_gml}posList')
                top_line = Line.from_coo(np.fromstring(top_pos.text,
                                                       sep=' ').reshape(-1, 3))

                bottom_edge = geom.find(f'{xmlns}faultBottomEdge')
                bottom_ls = bottom_edge.find(f'{xmlns_gml}LineString')
                bottom_pos = bottom_ls.find(f'{xmlns_gml}posList')
                bottom_line = Line.from_coo(np.fromstring(bottom_pos.text,
                                                          sep=' ').reshape(-1, 3))

                try:
                    ComplexFaultSurface.check_aki_richards_convention(
                        [top_line, bottom_line])

                except ValueError:
                    for edge_type in ['faultTopEdge', 'intermediateEdge',
                                      'faultBottomEdge']:
                        edge = geom.find(f'{xmlns}{edge_type}')
                        if edge:
                            ls = edge.find(f'{xmlns_gml}LineString')
                            pos = ls.find(f'{xmlns_gml}posList')
                            array = np.fromstring(pos.text, sep=' ').reshape(-1, 3)
                            array = np.flipud(array)
                            pos.text = ' '.join(array.ravel().astype(str))

            vparser = nrml.ValidatingXmlParser(nrml.validators, stop=None)
            src_model = vparser.parse_bytes(etree.tostring(root))
            xmlns = src_model.tag.split('}')[0][1:]
            src_model['xmlns'] = xmlns
            src_model['xmlns:gml'] = xmlns_gml
            output = nrml.get_source_model_04(src_model[0], fname,
                                              converter=parser)
            src_models.append(output)

    return src_models


def parse_srcs(source_model, trt=None):
    """
    From a source model, gets all the possible sources in order. If a tectonic
    region type is passed, return only sources pertaining to the trt.

    :param source_model: OpenQuake Source Model object
    :param trt: Tectonic region type  #todo pending
    :return:
    """

    sources = []
    if not isinstance(source_model, list):
        source_model = [source_model]

    for model in source_model:
        for grp in model.src_groups:
            for src in grp.sources:
                sources.append(src)

    src_types = set([i.__class__.__name__ for i in sources])
    log.info(f'Source types found: {src_types} ')
    return sources


def parse_region(sources, dh=0.1):
    """
    Creates a CSEP region from a Source Model. A Lat/Lon uniform grid is
     created from the boundaries of the SM. Active cells are determined if
     a source touches the cell.
     The final region is conformed from the Union of all sub-models.

    :param sources: list of sources
    :param dh: cell size
    :return:
    """

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

    # Create CSEP Grid from the Source Model's bounds and grid size
    x = cleaner_range(np.round(bounds[0], 1) - dh,
                      np.round(bounds[2], 1) + dh, h=dh)
    y = cleaner_range(np.round(bounds[1], 1) - dh,
                      np.round(bounds[3], 1) + dh, h=dh)
    grid = np.vstack([i.ravel() for i in np.meshgrid(x, y)]).T
    box_polygons = [shapely.geometry.Polygon(
        [(i[0], i[1]), (i[0] + dh, i[1]), (i[0] + dh, i[1] + dh),
         (i[0], i[1] + dh)]) for i in grid]
    box_geodf = gpd.GeoSeries(box_polygons)

    # Initialize array of active CSEP cells
    tag_cells = np.zeros(len(grid))

    # Check which cells are touched by Polygon-type Sources
    polygons = [shapely.geometry.Polygon(
        [(i, j) for i, j in zip(src.polygon.lons, src.polygon.lats)])
        for src in [*area_srcs, *sf_srcs, *cf_srcs]]
    for poly in polygons:
        tag_cells = np.logical_or(tag_cells, box_geodf.intersects(poly))

    # Check which cells are touched by Point-type Sources
    if pointsrc_coords.size > 0:
        df = pd.DataFrame(
            {'lon': pointsrc_coords[:, 0], 'lat': pointsrc_coords[:, 1]})
        df['coords'] = list(zip(df['lon'], df['lat']))
        df['coords'] = df['coords'].apply(shapely.geometry.Point)
        points = gpd.GeoDataFrame(df, geometry='coords')
        cells = gpd.GeoDataFrame(geometry=box_geodf)
        points_in_reg = gpd.tools.sjoin(points, cells, predicate='within',
                                        how='left')
        idxs = np.unique(points_in_reg.index_right.to_numpy())

        tag_from_points = np.zeros(len(grid))
        tag_from_points[idxs] = 1
        tag_cells = np.logical_or(tag_cells, tag_from_points)

    # Crops CSEP grid to those that were touched by sources
    if tag_cells.any():
        grid = grid[tag_cells]

    csep_region = CartesianGrid2D.from_origins(grid, dh=dh)

    return grid, csep_region


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
        region1 = CartesianGrid2D.from_origins(origins)
    for reg in regions:
        if isinstance(reg, str):
            origins = np.loadtxt(reg)
            reg = CartesianGrid2D.from_origins(origins)
        other_.append(reg)

    reg1_poly = [shapely.geometry.Polygon(
        [i.points[0], i.points[3], i.points[2], i.points[1]])
                    for i in region1.polygons]
    reg1_gdf = gpd.GeoSeries(reg1_poly)

    tag_cells = np.ones(len(reg1_poly))
    for reg2 in other_:
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
    dh = np.diff(np.unique(grid[:, 0])).max()
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



