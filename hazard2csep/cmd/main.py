import os
import logging
import argparse

import numpy
import matplotlib.pyplot as plt

import csep
from hazard2csep import __version__
from hazard2csep import region_lib
from hazard2csep import sm_lib
from hazard2csep import forecast_lib

log = logging.getLogger('hazard2csepLogger')


def region(files,
           dest=False,
           plot=False,
           fill=False, **_):

    log.info(f'CSEP: OpenQuake reader v{__version__} | Region parsing')
    if isinstance(files, list):
        log.info('Joining regions from models:')
        for file in files:
            log.info(f'\t> {file}')
    else:
        log.info('Creating region from model:')
        log.info(f'\t{files}')
    src_model = sm_lib.parse_source_model(files)
    srcs = sm_lib.parse_srcs(src_model)
    grid, csep_reg = region_lib.make_region(srcs, fill=fill)
    if not dest:
        dest = 'region.txt'
    numpy.savetxt(dest, grid, fmt='%.2f')
    log.info(f'Saved region to: {dest}')

    if plot:
        dest = dest.split('.')[0] + '.png'
        region_lib.plot_region(grid, dest)
        log.info(f'Saved figure to: {dest}')

    log.info('Finalized')


def intersect(files,
           dest=False,
           plot=False,
           fill=False, **_):

    log.info(f'CSEP: OpenQuake reader v{__version__} | Intersect Regions')
    log.info('Intersecting regions from models:')
    log.info(f'\t{files}')
    grid, csep_reg = region_lib.intersect_region(*files)
    if not dest:
        dest = 'region.txt'
    numpy.savetxt(dest, grid, fmt='%.2f')
    log.info(f'Saved region to: {dest}')


def project(files,
            projection_region=False,
            min_mag=4.5,
            max_mag=9.0,
            dm=0.2,
            max_depth=100,
            dest='forecast.csv',
            plot=False, **_):

    log.info(f'CSEP: OpenQuake reader v{__version__} | Rate projection')

    log.info('Projecting source models:')
    if isinstance(files, list):
        for file in files:
            log.info(f'\t{file}')
    else:
        log.info(f'\t{files}')

    if projection_region:
        log.info(f'Loading region: {projection_region}')
        csep_reg = csep.core.regions.CartesianGrid2D.from_origins(
            numpy.loadtxt(projection_region)
        )
    else:
        csep_reg = None

    src_model = sm_lib.parse_source_model(files)
    srcs = sm_lib.parse_srcs(src_model)
    forecast = forecast_lib.return_rates(srcs,
                                         region=csep_reg,
                                         min_mag=min_mag,
                                         max_mag=max_mag,
                                         dm=dm,
                                         max_depth=max_depth)

    destpath = os.path.split(dest)[0]
    if destpath:
        os.makedirs(destpath, exist_ok=True)

    forecast_lib.write_forecast(forecast, dest)

    if plot:
        log.info(f'Plotting forecast')
        dest = dest[::-1].split('.', maxsplit=1)[1][::-1] + '.png'
        forecast.plot(plot_args={'region_border': False})
        plt.savefig(dest)

    log.info('Finalized')
    return forecast


def hazard2csep():
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    parser.add_argument('func', type=str,
                        choices=['region', 'intersect', 'project'],
                        help='Source model parsing options')
    parser.add_argument('files', metavar='files', type=str,     # todo write arg docs
                        nargs='+',
                        help='path to source model files')
    parser.add_argument('-f', '--fill',  action='store_true',
                        default=False, help="Fill holes inside a region")
    parser.add_argument('-d', '--dest',
                        help="Destination to save the forecast/grid")
    parser.add_argument('-p', '--plot', action='store_true',
                        default=False, help="Plot results flag")
    args = parser.parse_args()
    try:
        func = globals()[args.func]
        args.__delattr__('func')

    except AttributeError:
        raise AttributeError('Function not implemented')
    func(**vars(args))
