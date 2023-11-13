import os
from os.path import join, abspath, normpath
import numpy
from openquake.hazardlib import nrml, sourceconverter
from openquake.hazardlib.source.complex_fault import ComplexFaultSource
from openquake.hazardlib.source.simple_fault import SimpleFaultSource
from openquake.hazardlib.source.area import AreaSource
from openquake.hazardlib.source.point import PointSource
from openquake.hazardlib.source.multi_point import MultiPointSource
from openquake.hazardlib.geo.surface import ComplexFaultSurface, \
    SimpleFaultSurface
import numpy as np
import matplotlib.pyplot as plt
import pickle
from os import path

import shapely.geometry
import oq2csep
import csep


def parse_elem_geom(src):
    if isinstance(src, AreaSource):
        print('hei')
        # poly = Polygon([(i,j) for i,j in zip(src.polygon.lons, src.polygon.lats)])
        # points = poly.intersection(MeshGrid)
        # # print(points)
        # try:
        #     xy = np.array([i.coords.xy for i in points]).squeeze()
        # except TypeError:
        #     xy = np.array([points.coords.xy]).squeeze().reshape((-1, 2))
        # counts = src.get_annual_occurrence_rates()
        # rates = []
        # magnitudes = []
        #
        # mws = []
        # rates_m = []
        # for i, j in counts:
        #     mws.append(i)
        #     rates_m.append(j)
        #
        # for k in xy:
        #     if len(k) > 0:
        #         cells.append([k[0], k[1]])
        #         magnitudes.append(mws)
        #         rates.append(np.array(rates_m)/xy.shape[0])
        #         src_type.append('area_source')


    elif isinstance(src, PointSource):
        print('ho')
    # rates_m = []
    # mws = []
    # for i, j in src.get_annual_occurrence_rates():
    #     mws.append(i)
    #     rates_m.append(j)
    # cells.append([src.location.longitude, src.location.latitude])
    # magnitudes.append(mws)
    # rates.append(np.array(rates_m))
    # src_type.append('point_source')


def parse_region(fnames):
    pass


if __name__ == '__main__':

    pkg_path = oq2csep.__path__[0]
    ssm_path = normpath(join(pkg_path, '..', 'examples', 'eshm13',
                             'SHARE_OQ_input_20140807', 'seifa_model.xml'))
    fsbg_path = normpath(join(pkg_path, '..', 'examples', 'eshm13',
                              'SHARE_OQ_input_20140807',
                              'faults_backg_source_model.xml'))
    parser = sourceconverter.SourceConverter(
        area_source_discretization=10,
        width_of_mfd_bin=0.1
    )

    # a = nrml.read_source_models([ssm_path], parser)
    a = nrml.read_source_models([fsbg_path], parser)
    b = list(a)
    # with open('test', 'wb') as file_:
    #     pickle.dump(b, file_ )
    # with open('test', 'rb') as file_:
    #     b = pickle.load(file_)
    # sources = []
    # for model in b:
    #     for grp in model.src_groups:
    #         for src in grp.sources:
    #             sources.append(src)

    points = []
    polygons = []

    points = shapely.geometry.MultiPoint(
        [(src.location.longitude, src.location.latitude)
         for src in sources if isinstance(src, PointSource)]
    )
    polygons = [shapely.geometry.Polygon(
        [(i, j) for i, j in zip(src.polygon.lons, src.polygon.lats)])
        for src in sources if isinstance(src, AreaSource)]


    # for src in sources:
    #     if isinstance(src, PointSource):
    #         point = shapely.geometry.Point(src.location.longitude,
    #                                        src.location.latitude)
    #         points.append(point)
    #     if isinstance(src, AreaSource):
    #         poly = shapely.geometry.Polygon(
    #             [(i, j) for i, j in zip(src.polygon.lons, src.polygon.lats)])
    #         polygons.append(poly)
    # points = shapely.geometry.MultiPoint(points)
