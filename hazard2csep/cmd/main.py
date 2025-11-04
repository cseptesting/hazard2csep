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
           fill=False,
           fault_buffer=0,
           shapefile=False,
           **_):

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
    grid, csep_reg = region_lib.make_region(srcs,
                                            fault_buffer=fault_buffer,
                                            fill=fill,
                                            shapefile=shapefile)
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
              shapefile=False,
              **_):

    log.info(f'CSEP: OpenQuake reader v{__version__} | Intersect Regions')
    log.info('Intersecting regions from models:')
    log.info(f'\t{files}')
    grid, csep_reg = region_lib.intersect_region(*files, shapefile=shapefile)
    if not dest:
        dest = 'region.txt'
    numpy.savetxt(dest, grid, fmt='%.2f')
    log.info(f'Saved region to: {dest}')
    if plot:
        dest = dest.split('.')[0] + '.png'
        region_lib.plot_region(grid, dest)
        log.info(f'Saved figure to: {dest}')


def project(files,
            projection_region=False,
            min_mag=4.5,
            max_mag=9.0,
            dm=0.2,
            dh=0.1,
            max_depth=100,
            buffer=0,
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
            numpy.loadtxt(projection_region), dh=dh
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
                                         buffer=buffer,
                                         max_depth=max_depth)

    destpath = os.path.split(dest)[0]
    if destpath:
        os.makedirs(destpath, exist_ok=True)

    forecast_lib.write_forecast(forecast, dest)

    if plot:
        log.info(f'Plotting forecast')
        dest = dest[::-1].split('.', maxsplit=1)[1][::-1] + '.png'
        forecast.plot(plot_args={'region_border': False, 'borders': True, "basemap": "ESRI_terrain"})
        plt.savefig(dest)

    log.info('Finalized')
    return forecast


def hazard2csep():
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    sub = parser.add_subparsers(dest="func", required=True)

    # project subcommand
    p = sub.add_parser("project", help="Project source model rates")
    p.add_argument("files", nargs="+", help="Path(s) to source model files")
    p.add_argument("-r", "--projection-region", dest="projection_region", help="Path to region origins file")
    p.add_argument("--min-mag", type=float, default=4.5)
    p.add_argument("--max-mag", type=float, default=9.0)
    p.add_argument("--dm", type=float, default=0.2)
    p.add_argument("--dh", type=float, default=0.1)
    p.add_argument("--max-depth", type=float, default=100)
    p.add_argument("--buffer", type=float, default=0)
    p.add_argument("-d", "--dest", default="forecast.csv")
    p.add_argument("--plot", action="store_true", default=False)

    # region subcommand (example — add the args that function needs)
    r = sub.add_parser("region", help="Region operations")
    r.add_argument("files", nargs="+")
    r.add_argument("-f", "--fill", action="store_true", default=False)
    r.add_argument("--dest")

    # intersect subcommand (example)
    i = sub.add_parser("intersect", help="Intersect operations")
    i.add_argument("files", nargs="+")
    i.add_argument("--dest")
    i.add_argument("--plot", action="store_true", default=False)

    args = parser.parse_args()
    func_name = args.func
    kwargs = vars(args)
    kwargs.pop("func")
    try:
        func = globals()[func_name]
    except KeyError:
        raise AttributeError("Function not implemented")
    func(**kwargs)
