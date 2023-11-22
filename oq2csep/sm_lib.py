import os.path
from xml.etree import ElementTree as etree

import numpy as np
from openquake.hazardlib import sourceconverter, nrml
from openquake.hazardlib.geo import Line, ComplexFaultSurface
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

    log.info('Parsing source model: %s' % fname)
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
    log.info('Collecting sources from source model')
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


def parse_logictree_files(filename):
    " only linear filetress with 1 branchset"
    lt = nrml.read(filename)[0]
    lts = lt[0]

    if not ('logicTreeBranchSet' in lts.tag):
        lts = lts[0]

    branches = {}
    for ltbranch in lts:
        files = ltbranch[0].text.split(' ')
        branches[ltbranch.attrib['branchID']] = files

    return os.path.dirname(filename), branches
    # return os.path.dirname(filename), {i['branchID']: i[0].text.split(' ') for i in lts}

