import os

import numpy as np
import matplotlib.pyplot as plt

from openquake.hazardlib.mfd import TruncatedGRMFD, EvenlyDiscretizedMFD

from oq2csep import region_lib
from oq2csep import sm_lib
from oq2csep import forecast_lib
from oq2csep.cmd import main

_SRC_FOLDERS = [os.path.join(i) for i in [
                                  'area_src',
                                  'complexfault_src',
                                  'simplefault_src',
                                  'point_src',
                                  'multipoint_src',
                                  'simplefault_areasource_src',
                                  'simplefault_multipoint_src',
                                  'complexfault_crop_src'
            ]
                ]
_REGIONS = [os.path.join('regions', i) for i in [
    'test_simple.txt',
    'test_simple_fill.txt',
    'test_multiregion.txt',
    'test_multiregion_fill.txt']]


def test_region_fill():

    for reg in _REGIONS:
        coords = np.genfromtxt(reg)
        new_coords = region_lib.fill_holes(coords, dh=1)

        plt.figure(figsize=(5,5))
        plt.plot(*coords.T, 'bo', ms=7)
        plt.plot(*new_coords.T, 'r^', ms=4)
        filename = reg.replace('.txt', '.png')
        plt.savefig(filename)


def test_regions():

    for sm_dir in _SRC_FOLDERS:
        sm_file = os.path.join(sm_dir, 'test_src.xml')
        dest = os.path.join(sm_dir, 'region.txt')
        main.region(sm_file, plot=True, dest=dest)


def test_projection_rates():

    for _DIR in _SRC_FOLDERS[:-1]:  # Skip cropped ComplexFaultSource
        _PATH = os.path.join(_DIR, 'test_src.xml')
        src_model = sm_lib.parse_source_model(_PATH)

        for src in sm_lib.parse_srcs(src_model):

            min_mag = 5.0
            max_mag = 8.4
            dm = 0.1
            magnitudes = np.arange(min_mag, max_mag + dm / 2, dm)

            forecast = forecast_lib.return_rates([src],
                                                 min_mag=min_mag,
                                                 max_mag=max_mag,
                                                 dm=dm)
            mag_rates = src.get_annual_occurrence_rates()

            src_mags = np.array([i[0] for i in mag_rates])
            src_rates = np.array([i[1] for i in mag_rates])
            used_rates = src_rates[src_mags >= min_mag]

            np.testing.assert_allclose(forecast.magnitudes, magnitudes)
            np.testing.assert_allclose(used_rates.sum(), forecast.data.sum())

            if isinstance(src.mfd, TruncatedGRMFD):
                a_val = src.mfd.a_val
                b_val = src.mfd.b_val
                mbin = src.mfd.bin_width

                src_mmin = min([i[0] for i in mag_rates]) - mbin / 2.
                src_mmax = max([i[0] for i in mag_rates]) + mbin / 2.

                src_total_rate = (
                        10 ** (a_val - b_val * max(min_mag, src_mmin)) -
                        10 ** (a_val - b_val * min(max_mag, src_mmax)))

                np.testing.assert_allclose(src_total_rate, forecast.data.sum())


def test_model_projection():

    for sm_dir in _SRC_FOLDERS:

        forecast = main.project(os.path.join(sm_dir, 'test_src.xml'),
                                min_mag=4.7, max_mag=8.1, dm=0.2, plot=True,
                                dest=os.path.join(sm_dir, 'forecast.csv'))


if __name__ == '__main__':
    pass
    # test_regions()
    # test_rate_projection()

