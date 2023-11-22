import unittest
import numpy as np
import os
from hazard2csep.region_lib import fill_holes, make_region, intersect_region, plot_region
from hazard2csep.sm_lib import parse_source_model, parse_srcs
import geopandas as gpd
import shapely

_SRC_FOLDERS = [os.path.join(i) for i in [
    'area_src',
    'complexfault_src',
    'simplefault_src',
    'point_src',
    'multipoint_src',
    'simplefault_areasource_src',
    'simplefault_multipoint_src',
    'complexfault_crop_src']]


class TestRegionLib(unittest.TestCase):

    def setUp(self):

        # This method will be used to set up any objects that you will use in multiple tests.
        # For example, you might want to create some coordinates here that you will use in multiple tests.
        pass

    @staticmethod
    def get_srcs(path_):
        file_ = os.path.join(path_, 'test_src.xml')
        return parse_srcs(parse_source_model(file_))

    def test_fill_holes(self):

        # Test when `coords` is a non-empty array, but there are no holes
        # to fill
        coords = np.array([[0, 0.0], [0.1, 0.0], [0.2, 0],
                           [0, 0.1], [0.1, 0.1], [0.2, 0.1],
                           [0, 0.2], [0.1, 0.2], [0.2, 0.2]])
        result = fill_holes(coords)
        np.testing.assert_allclose(result, coords)

        # Test when `coords` is a non-empty array and there are holes to fill
        coords = np.array([[0, 0.0], [0.1, 0.0], [0.2, 0],
                           [0, 0.1], [0.2, 0.1],
                           [0, 0.2], [0.1, 0.2], [0.2, 0.2]])
        expected_result = np.array([[0, 0.0], [0.1, 0.0], [0.2, 0.0],
                                    [0, 0.1], [0.1, 0.1], [0.2, 0.1],
                                    [0, 0.2], [0.1, 0.2], [0.2, 0.2],
                                    ])
        result = fill_holes(coords)
        np.testing.assert_allclose(result, expected_result)

        # Test when `dh` is provided
        coords = np.array([[0, 0], [1, 0], [2, 0],
                           [0, 1], [2, 1],
                           [0, 2], [1, 2], [2, 2]])
        expected_result = np.array([[0, 0], [1, 0], [2, 0],
                                    [0, 1], [1, 1], [2, 1],
                                    [0, 2], [1, 2], [2, 2]])
        result = fill_holes(coords, dh=1)
        np.testing.assert_array_equal(result, expected_result)

    def test_make_region(self):

        for sm_dir in _SRC_FOLDERS:
            srcs = self.get_srcs(sm_dir)
            target = os.path.join(sm_dir, 'region.txt')

            coords, _ = make_region(srcs)
            region = np.genfromtxt(target)
            np.testing.assert_allclose(coords, region, atol=1e-3)


    def test_intersect_region(self):
        # Here you will test the intersect_region function.
        # You should test various cases, including edge cases.
        pass


if __name__ == '__main__':
    unittest.main()