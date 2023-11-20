import os

import numpy as np
import matplotlib.pyplot as plt

from oq2csep import region_lib
from oq2csep.cmd import main

_SRC_FOLDERS = [os.path.join(i) for i in ['area_src', 'complexfault_src',
                                          'simplefault_src', 'point_src',
                                          'multipoint_src',
                                          'simplefault_areasource_src',
                                          'simplefault_multipoint_src']]
_REGIONS = [os.path.join('regions', i) for i in [
    'test_simple.txt',
    'test_simple_fill.txt',
    'test_multiregion.txt',
    'test_multiregion_fill.txt']]


def test_region_utils():

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


def test_rate_projection():

    for sm_dir in _SRC_FOLDERS:

        forecast = main.project(os.path.join(sm_dir, 'test_src.xml'),
                                min_mag=4.7, max_mag=8.1, dm=0.2, plot=True,
                                dest=os.path.join(sm_dir, 'forecast.txt'))


if __name__ == '__main__':
    test_regions()
    test_rate_projection()

