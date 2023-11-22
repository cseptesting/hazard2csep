import os
import logging
import argparse

import numpy
import matplotlib.pyplot as plt

import csep
from oq2csep import __version__
from oq2csep import region_lib
from oq2csep import sm_lib
from oq2csep import forecast_lib

log = logging.getLogger('oq2csepLogger')


def region(files, intersect=False, dest=False, plot=False, fill=False, **_):

    log.info(f'OpenQuake to CSEP v{__version__} | Region parsing')

    if intersect:
        log.info('Intersecting regions from models:')
        log.info(f'\t{files}')
        grid, csep_reg = region_lib.intersect_region(*files)
        if not dest:
            dest = 'region.txt'
        numpy.savetxt(dest, grid, fmt='%.2f')
        log.info(f'Saved region to: {dest}')

    else:
        log.info('Joining regions from models:')
        if isinstance(files, list):
            for file in files:
                log.info(f'\t{file}')
        else:
            log.info(f'\t{files}')
        log.info('Parsing source models')
        src_model = sm_lib.parse_source_model(files)
        log.info('Getting source elements')
        srcs = sm_lib.parse_srcs(src_model)
        log.info('Computing region')
        grid, csep_reg = region_lib.parse_region(srcs, fill=fill)
        if not dest:
            dest = 'region.txt'
        numpy.savetxt(dest, grid, fmt='%.2f')
        log.info(f'Saved region to: {dest}')

    if plot:
        dest = dest.replace('.txt', '.png')
        region_lib.plot_region(grid, dest)
        log.info(f'Saved figure to: {dest}')
        log.info('Close figure to continue')

    log.info('Finalized')


def project(files, region=False, min_mag=4.7, max_mag=8, dm=0.2,
            dest=False, plot=False, **_):

    log.info(f'OpenQuake to CSEP v{__version__} | Rate projection')

    log.info('Projecting source models:')
    if isinstance(files, list):
        for file in files:
            log.info(f'\t{file}')
    else:
        log.info(f'\t{files}')
    if region:
        log.info(f'Loading region: {region}')
        csep_reg = csep.core.regions.CartesianGrid2D.from_origins(
            numpy.loadtxt(region)
        )
    else:
        csep_reg = None
    src_model = sm_lib.parse_source_model(files)
    srcs = sm_lib.parse_srcs(src_model)
    forecast = forecast_lib.return_rates(srcs,
                                         region=csep_reg,
                                         min_mag=min_mag,
                                         max_mag=max_mag,
                                         dm=dm)

    if dest:
        forecast_lib.write_forecast(forecast, dest)

    if plot:
        log.info(f'Plotting forecast')
        if dest is False:
            dest = 'projected_forecast.png'
        else:
            dest = dest.split('.')[0] + '.png'

        forecast.plot(plot_args={'region_border': False})
        plt.savefig(dest)

        data_m65 = forecast.data[:, forecast.magnitudes > 6.5].sum(axis=1)
        csep.utils.plots.plot_spatial_dataset(
            forecast.region.get_cartesian(numpy.log10(data_m65)),
            forecast.region,
            plot_args={'region_border': False})
        plt.savefig(dest.replace('.png', '_m65.png'))

    log.info('Finalized')
    return forecast


def oq2csep():
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    parser.add_argument('func', type=str,
                        choices=['region', 'project'],
                        help='Source model parsing options')
    parser.add_argument('files', metavar='files', type=str,     # todo write arg docs
                        nargs='+',
                        help='path to source model files')
    parser.add_argument('-i', '--intersect', action='store_true',
                        default=False, help="Timestamp results")
    parser.add_argument('-f', '--fill',
                        default=False, help="File to save the grid")
    parser.add_argument('-d', '--dest',
                        default=False, help="File to save the grid")
    parser.add_argument('-p', '--plot', action='store_true',
                        default=False, help="Plot esults")
    args = parser.parse_args()
    try:
        func = globals()[args.func]
        args.__delattr__('func')

    except AttributeError:
        raise AttributeError('Function not implemented')
    func(**vars(args))
