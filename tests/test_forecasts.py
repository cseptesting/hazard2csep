import os
import tempfile

import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon, Point

from csep.core.regions import CartesianGrid2D
import csep.models
from csep.core.forecasts import GriddedForecast

from hazard2csep import sm_lib
from hazard2csep import forecast_lib


# -------------------------------------------------------------------------
# Fixtures / paths (mirror test_sources.py style)
# -------------------------------------------------------------------------
_SRC_FOLDERS = [os.path.join(i) for i in [
    'area_src',
    'complexfault_src',
    'simplefault_src',
    'point_src',
    'multipoint_src',
    'simplefault_areasource_src',
    'simplefault_multipoint_src',
    'complexfault_crop_src'
]]


def _test_src_xml_paths():
    for folder in _SRC_FOLDERS:
        yield os.path.join(folder, "test_src.xml")


# -------------------------------------------------------------------------
# project_mfd
# -------------------------------------------------------------------------
def test_project_mfd_center_simple():
    print("\n[forecast_lib.project_mfd] test_project_mfd_center_simple")

    # Three source bins, centers at 5.05, 5.15, 5.25 with rate 1 each
    src = np.array([
        [5.05, 1.0],
        [5.15, 1.0],
        [5.25, 1.0],
    ])
    # Target bins left edges: [5.0, 5.1, 5.2]
    mags = np.array([5.0, 5.1, 5.2])

    expected = np.array([1.0, 1.0, 1.0])  # one center per bin
    actual = forecast_lib.project_mfd(src, mags, bin_width=0.1, method="center")

    print("  expected:", expected)
    print("  actual:  ", actual)

    np.testing.assert_allclose(expected, actual)
    np.testing.assert_allclose(actual.sum(), src[:, 1].sum())


def test_project_mfd_overlap_rebin_preserves_total():
    print("\n[forecast_lib.project_mfd] test_project_mfd_overlap_rebin_preserves_total")

    # Source centers every 0.1, want to rebin to 0.2 with overlap
    centers = np.array([5.05, 5.15, 5.25, 5.35])
    values = np.array([1.0, 2.0, 3.0, 4.0])  # arbitrary
    src = np.column_stack([centers, values])

    # Target left edges with width 0.2: bins [5.0-5.2), [5.2-5.4)
    mags = np.array([5.0, 5.2])

    total_before = values.sum()
    actual = forecast_lib.project_mfd(src, mags, bin_width=0.1, method="overlap")
    total_after = actual.sum()

    print("  total_before:", total_before)
    print("  total_after: ", total_after)

    # Total rate must be preserved by construction
    np.testing.assert_allclose(total_before, total_after)


def test_project_mfd_overlap_single_bin_uses_bin_width():
    print("\n[forecast_lib.project_mfd] test_project_mfd_overlap_single_bin_uses_bin_width")

    # Single source bin, center 5.05, width = 0.1 (from bin_width argument)
    src = np.array([[5.05, 10.0]])
    mags = np.array([5.0, 5.1])  # two target bins of width 0.1

    # With perfect overlap, the 10.0 should fall entirely into first bin
    actual = forecast_lib.project_mfd(src, mags, bin_width=0.1, method="overlap")

    print("  actual:", actual)
    np.testing.assert_allclose(actual.sum(), 10.0)
    assert actual[0] > 0.0
def test_project_mfd_overlap_nonuniform_target_raises():
    print("\n[forecast_lib.project_mfd] test_project_mfd_overlap_nonuniform_target_raises")

    src = np.array([[5.05, 1.0], [5.15, 2.0]])
    # truly non-uniform spacing: diffs [0.15, 0.17]
    mags = np.array([5.0, 5.15, 5.32])

    try:
        forecast_lib.project_mfd(src, mags, bin_width=0.1, method="overlap")
    except ValueError as exc:
        print("  caught ValueError as expected:", exc)
    else:
        raise AssertionError("Expected ValueError for non-uniform target magnitudes")


def test_project_mfd_preserves_total_for_real_source():
    from hazard2csep import sm_lib, forecast_lib

    min_mag = 5.0
    max_mag = 8.4
    dm = 0.1

    sm_file = os.path.join("area_src", "test_src.xml")  # pick any of your dirs
    src_model = sm_lib.parse_source_model(sm_file)
    src = sm_lib.parse_srcs(src_model)[0]

    # OQ MFD
    mag_rates = src.get_annual_occurrence_rates()
    src_mags = np.array([m for (m, _) in mag_rates])
    src_rates = np.array([r for (_, r) in mag_rates])

    # CSEP target grid
    csep_mags = sm_lib.cleaner_range(min_mag, max_mag, dm).round(1)

    # Use your project_mfd
    rates_array = np.column_stack([src_mags, src_rates])
    rebinned = forecast_lib.project_mfd(rates_array, csep_mags, src.mfd.bin_width)

    print("sum src_rates[m>=min_mag]:", src_rates[src_mags >= min_mag].sum())
    print("sum rebinned:", rebinned.sum())

    np.testing.assert_allclose(
        src_rates[src_mags >= min_mag].sum(),
        rebinned.sum(),
        rtol=1e-5, atol=1e-8
    )

# -------------------------------------------------------------------------
# intersect_geom2grid
# -------------------------------------------------------------------------
def _make_simple_grid_gdf():
    # 2x2 grid: cells [0,1)x[0,1), [1,2)x[0,1), [0,1)x[1,2), [1,2)x[1,2)
    polys = [
        Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
        Polygon([(1, 0), (2, 0), (2, 1), (1, 1)]),
        Polygon([(0, 1), (1, 1), (1, 2), (0, 2)]),
        Polygon([(1, 1), (2, 1), (2, 2), (1, 2)]),
    ]
    return gpd.GeoDataFrame({"geometry": polys})

def test_intersect_geom2grid_polygons_multiple_cells():
    print("\n[forecast_lib.intersect_geom2grid] test_intersect_geom2grid_polygons_multiple_cells")

    csep_grid = _make_simple_grid_gdf()
    # Polygon covering bottom row (y in [0,1])
    big_poly = Polygon([(0, 0), (2, 0), (2, 1), (0, 1)])

    maps = forecast_lib.intersect_geom2grid([big_poly], csep_grid)
    idxs, frac = maps[0]

    print("  indices:", idxs)
    print("  frac:   ", frac)
    print("  sum(frac):", frac.sum())

    idxs_list = idxs.tolist() if isinstance(idxs, np.ndarray) else [idxs]
    idxs_set = set(idxs_list)

    # We *at least* want the two bottom cells in the intersection
    assert {0, 1}.issubset(idxs_set)

    # Fractions should sum to ~1 because grid fully covers the source polygon
    np.testing.assert_allclose(frac.sum(), 1.0, rtol=1e-6, atol=1e-8)

    # Cells 2 and 3 only touch along the boundary -> intersection area should be 0
    # (or numerically extremely close to 0)
    for i, idx in enumerate(idxs_list):
        if idx in (2, 3):
            print(f"  cell {idx} frac (should be ~0):", frac[i])
            assert abs(frac[i]) < 1e-10




def test_intersect_geom2grid_points_basic():
    print("\n[forecast_lib.intersect_geom2grid] test_intersect_geom2grid_points_basic")

    csep_grid = _make_simple_grid_gdf()
    # Point in cell 0 and point outside any cell
    p_inside = Point(0.5, 0.5)
    p_outside = Point(10.0, 10.0)

    maps = forecast_lib.intersect_geom2grid([p_inside, p_outside], csep_grid)

    for i, (idxs, frac) in enumerate(maps):
        print(f"  point {i} -> idxs: {idxs}, frac: {frac}")

    # First point should map to exactly one cell with frac = 1
    assert isinstance(maps[0][0], np.ndarray)
    assert maps[0][0].size == 1
    assert maps[0][1] == 1

    # Second point should map to empty (outside grid)
    assert maps[1][0].size == 0
    assert maps[1][1].size == 0


# -------------------------------------------------------------------------
# read_forecast + write_forecast (round-trip)
# -------------------------------------------------------------------------
def test_write_read_forecast_roundtrip(tmp_path):
    print("\n[forecast_lib.write_forecast + read_forecast] test_write_read_forecast_roundtrip")

    # Build a tiny 1x2 grid (two cells) with dh = 0.1
    bbox1 = ((0.0, 0.0), (0.1, 0.0), (0.1, 0.1), (0.0, 0.1))
    bbox2 = ((0.1, 0.0), (0.2, 0.0), (0.2, 0.1), (0.1, 0.1))
    polys = [csep.models.Polygon(bbox1), csep.models.Polygon(bbox2)]
    dh = 0.1

    region = CartesianGrid2D(polys, dh)
    mags = np.array([5.0, 5.1])
    data = np.array([[0.1, 0.2],
                     [0.3, 0.4]])

    original = GriddedForecast(data=data, region=region, magnitudes=mags)

    dest = tmp_path / "forecast.csv"
    forecast_lib.write_forecast(original, dest=str(dest), fmt="csep", depths=(0, 30), dh=dh)

    # Now read it back
    reloaded = forecast_lib.read_forecast(str(dest))

    print("  original magnitudes:", original.magnitudes)
    print("  reloaded magnitudes:", reloaded.magnitudes)
    print("  original data:\n", original.data)
    print("  reloaded data:\n", reloaded.data)

    np.testing.assert_allclose(original.magnitudes, reloaded.magnitudes)
    np.testing.assert_allclose(original.data, reloaded.data)


# -------------------------------------------------------------------------
# return_rates (similar to test_sources but more debug; real XML)
# -------------------------------------------------------------------------
def test_return_rates_per_source_matches_total_rate():
    print("\n[forecast_lib.return_rates] test_return_rates_per_source_matches_total_rate")

    min_mag = 5.0
    max_mag = 8.4
    dm = 0.1
    expected_mags = np.arange(min_mag, max_mag + dm / 2, dm)

    for sm_dir in _SRC_FOLDERS[:-1]:  # skip cropped complexfault as test_sources does
        sm_file = os.path.join(sm_dir, "test_src.xml")
        print(f"\n  -> source model file: {sm_file}")

        src_model = sm_lib.parse_source_model(sm_file)
        for src in sm_lib.parse_srcs(src_model):
            print(f"    source_id: {getattr(src, 'source_id', 'N/A')}, type: {src.__class__.__name__}")

            forecast = forecast_lib.return_rates([src],
                                                 min_mag=min_mag,
                                                 max_mag=max_mag,
                                                 dm=dm)

            mag_rates = src.get_annual_occurrence_rates()
            src_mags = np.array([m for (m, _) in mag_rates])
            src_rates = np.array([r for (_, r) in mag_rates])

            used_rates = src_rates[src_mags >= min_mag]

            print("      expected magnitudes:", expected_mags)
            print("      forecast magnitudes:", forecast.magnitudes)
            print("      sum used_rates:", used_rates.sum())
            print("      sum forecast.data:", forecast.data.sum())

            np.testing.assert_allclose(forecast.magnitudes, expected_mags)
            np.testing.assert_allclose(used_rates.sum(), forecast.data.sum(), rtol=1e-5, atol=1e-8)


def test_return_rates_shape_and_nonzero_for_simple_case():
    print("\n[forecast_lib.return_rates] test_return_rates_shape_and_nonzero_for_simple_case")

    # Use the first source of the simplefault_src model as a representative
    sm_file = os.path.join("simplefault_src", "test_src.xml")
    src_model = sm_lib.parse_source_model(sm_file)
    srcs = sm_lib.parse_srcs(src_model)

    assert len(srcs) > 0
    src = srcs[0]

    min_mag = 4.7
    max_mag = 8.1
    dm = 0.2

    forecast = forecast_lib.return_rates([src],
                                         min_mag=min_mag,
                                         max_mag=max_mag,
                                         dm=dm)

    print("  forecast magnitudes:", forecast.magnitudes)
    print("  forecast data shape:", forecast.data.shape)
    print("  total events:", forecast.data.sum())

    # Shape: (n_cells, n_mags)
    assert forecast.data.ndim == 2
    assert forecast.data.shape[1] == len(forecast.magnitudes)
    # Some non-zero content expected
    assert forecast.data.sum() > 0.0

import os
import numpy as np

from openquake.hazardlib.source.area import AreaSource
from openquake.hazardlib.source.simple_fault import SimpleFaultSource
from openquake.hazardlib.source.point import PointSource
from openquake.hazardlib.source.multi_point import MultiPointSource

from hazard2csep import sm_lib
from hazard2csep import forecast_lib


# Reuse same structure as test_sources.py
_SRC_FOLDERS = [
    "area_src",
    "complexfault_src",
    "simplefault_src",
    "point_src",
    "multipoint_src",
    "simplefault_areasource_src",
    "simplefault_multipoint_src",
    "complexfault_crop_src",  # we'll skip this one where needed
]


def _load_first_src_of_type(folder, cls):
    """Helper: load first source of given class from test_src.xml in folder."""
    sm_file = os.path.join(folder, "test_src.xml")
    src_models = sm_lib.parse_source_model(sm_file)
    srcs = sm_lib.parse_srcs(src_models)
    for src in srcs:
        if isinstance(src, cls):
            return src
    raise RuntimeError(f"No {cls.__name__} found in {sm_file}")


def _compute_expected_truncated_rate(src, min_mag):
    """Sum of src rates for mags >= min_mag."""
    mag_rates = src.get_annual_occurrence_rates()
    mags = np.array([m for (m, _) in mag_rates])
    vals = np.array([v for (_, v) in mag_rates])
    used = vals[mags >= min_mag]
    return used.sum()


def test_return_rates_area_source_from_xml():
    print("\n[forecast_lib.return_rates] test_return_rates_area_source_from_xml")
    folder = "area_src"
    src = _load_first_src_of_type(folder, AreaSource)

    min_mag = 5.0
    max_mag = 8.4
    dm = 0.1

    expected_total = _compute_expected_truncated_rate(src, min_mag)
    forecast = forecast_lib.return_rates([src], min_mag=min_mag, max_mag=max_mag, dm=dm)

    actual_total = forecast.data.sum()
    print("  folder:", folder)
    print("  expected_total_rate (m>=5.0):", expected_total)
    print("  actual_total_rate (forecast):", actual_total)

    # Magnitude grid should match np.arange definition
    mags_ref = np.arange(min_mag, max_mag + dm / 2, dm)
    print("  forecast.magnitudes:", forecast.magnitudes)
    print("  mags_ref:", mags_ref)
    np.testing.assert_allclose(forecast.magnitudes, mags_ref)

    np.testing.assert_allclose(expected_total, actual_total, rtol=1e-6)


def test_return_rates_simple_fault_from_xml():
    print("\n[forecast_lib.return_rates] test_return_rates_simple_fault_from_xml")
    folder = "simplefault_src"
    src = _load_first_src_of_type(folder, SimpleFaultSource)

    min_mag = 5.0
    max_mag = 8.4
    dm = 0.1

    expected_total = _compute_expected_truncated_rate(src, min_mag)
    forecast = forecast_lib.return_rates([src], min_mag=min_mag, max_mag=max_mag, dm=dm)

    actual_total = forecast.data.sum()
    print("  folder:", folder)
    print("  expected_total_rate (m>=5.0):", expected_total)
    print("  actual_total_rate (forecast):", actual_total)

    np.testing.assert_allclose(expected_total, actual_total, rtol=1e-6)


def test_return_rates_point_source_from_xml():
    print("\n[forecast_lib.return_rates] test_return_rates_point_source_from_xml")
    folder = "point_src"
    src = _load_first_src_of_type(folder, PointSource)

    min_mag = 5.0
    max_mag = 8.4
    dm = 0.1

    expected_total = _compute_expected_truncated_rate(src, min_mag)
    forecast = forecast_lib.return_rates([src], min_mag=min_mag, max_mag=max_mag, dm=dm)

    actual_total = forecast.data.sum()
    print("  folder:", folder)
    print("  expected_total_rate (m>=5.0):", expected_total)
    print("  actual_total_rate (forecast):", actual_total)

    np.testing.assert_allclose(expected_total, actual_total, rtol=1e-6)


def test_return_rates_multipoint_source_from_xml():
    print("\n[forecast_lib.return_rates] test_return_rates_multipoint_source_from_xml")
    folder = "multipoint_src"
    src = _load_first_src_of_type(folder, MultiPointSource)

    min_mag = 5.0
    max_mag = 8.4
    dm = 0.1

    expected_total = _compute_expected_truncated_rate(src, min_mag)
    forecast = forecast_lib.return_rates([src], min_mag=min_mag, max_mag=max_mag, dm=dm)

    actual_total = forecast.data.sum()
    print("  folder:", folder)
    print("  expected_total_rate (m>=5.0):", expected_total)
    print("  actual_total_rate (forecast):", actual_total)

    np.testing.assert_allclose(expected_total, actual_total, rtol=1e-6)
def test_multipoint_rate_conservation_above_min_mag():
    src_model = sm_lib.parse_source_model("multipoint_src/test_src.xml")
    mps = sm_lib.parse_srcs(src_model)[0]  # MultiPointSource

    # Extract exactly as get_rate_mpoint_source does
    points, rates_list, bin_widths = forecast_lib.get_rate_mpoint_source([mps])

    min_mag = 5.0
    max_mag = 8.4
    dm      = 0.1
    mags    = np.arange(min_mag, max_mag + dm / 2, dm)

    total_true_above_min = 0.0
    total_rebinned = 0.0

    for rate, bw in zip(rates_list, bin_widths):
        centers = np.array([m for m, _ in rate])
        values  = np.array([v for _, v in rate])

        total_true_above_min += values[centers >= min_mag].sum()

        reb = forecast_lib.project_mfd(rate, mags, bin_width=bw)
        total_rebinned += reb.sum()

    assert np.isclose(total_rebinned, total_true_above_min)
def test_complex_fault_rate_conservation():
    src_model = sm_lib.parse_source_model("complexfault_src/test_src.xml")
    src = sm_lib.parse_srcs(src_model)[0]  # ComplexFaultSource

    polygons, rates_list, bin_widths = forecast_lib.get_rate_complex_fault(
        [src],
        max_depth=200
    )

    min_mag = 5.0
    max_mag = 8.4
    dm      = 0.1
    mags    = np.arange(min_mag, max_mag + dm / 2, dm)

    total_true = 0.0
    total_rebinned = 0.0

    for rate, bw in zip(rates_list, bin_widths):
        centers = np.array([m for m, _ in rate])
        values  = np.array([v for _, v in rate])
        total_true += values.sum()
        total_rebinned += forecast_lib.project_mfd(rate, mags, bw).sum()

    assert np.isclose(total_rebinned, total_true)


def test_projection_rates_multi_source():
    for _DIR in _SRC_FOLDERS:
        _PATH = os.path.join(_DIR, 'test_src.xml')
        src_model = sm_lib.parse_source_model(_PATH)
        srcs = list(sm_lib.parse_srcs(src_model))

        min_mag = 5.0
        max_mag = 8.4
        dm = 0.1

        forecast_all = forecast_lib.return_rates(
            srcs, min_mag=min_mag, max_mag=max_mag, dm=dm
        )

        # Sum of src_mfd rates above min_mag vs forecast_total
        src_total = 0.0
        for src in srcs:
            mag_rates = src.get_annual_occurrence_rates()
            src_mags = np.array([m for m, _ in mag_rates])
            src_rates = np.array([r for _, r in mag_rates])
            src_total += src_rates[src_mags >= min_mag].sum()

        np.testing.assert_allclose(src_total, forecast_all.data.sum(),
                                   atol=1e-3,
                                   rtol=5e-3)

        # Order invariance
        shuffled = srcs[::-1]  # or random.shuffle copy
        forecast_shuffled = forecast_lib.return_rates(
            shuffled, min_mag=min_mag, max_mag=max_mag, dm=dm
        )
        np.testing.assert_allclose(
            forecast_all.data, forecast_shuffled.data
        )
        forecast, diagnostics = forecast_lib.return_rates(
            srcs,
            min_mag=min_mag,
            max_mag=max_mag,
            dm=dm,
            return_diag=True,
        )

        group_rebinned_total = sum(g["rebinned_total_rate"] for g in diagnostics["groups"])

        np.testing.assert_allclose(group_rebinned_total, forecast.data.sum())


def test_projection_rates_mag_window_consistency():
    for _DIR in _SRC_FOLDERS:
        print(_DIR)
        _PATH = os.path.join(_DIR, 'test_src.xml')
        src_model = sm_lib.parse_source_model(_PATH)
        srcs = list(sm_lib.parse_srcs(src_model))

        min_mag = 5.0
        max_mag = 8.4
        dm = 0.1

        # IMPORTANT: use keyword args
        fore_full = forecast_lib.return_rates(
            srcs,
            min_mag=min_mag,
            max_mag=max_mag,
            dm=dm,
        )

        # A narrower magnitude window
        min2, max2 = 5.4, 8.0
        fore_sub = forecast_lib.return_rates(
            srcs,
            min_mag=min2,
            max_mag=max2,
            dm=dm,
        )

        # Check that magnitudes line up
        mask = (fore_full.magnitudes >= min2) & (fore_full.magnitudes <= max2)
        np.testing.assert_allclose(
            fore_full.magnitudes[mask],
            fore_sub.magnitudes,
        )

        # And that data are identical in the overlapping window
        np.testing.assert_allclose(
            fore_full.data[:, mask],
            fore_sub.data,
            rtol=1e-6,
            atol=1e-9,

        )


from openquake.hazardlib.source import PointSource
from openquake.hazardlib.mfd import TruncatedGRMFD
from openquake.hazardlib.geo.point import Point

def make_toy_point_source():
    mfd = TruncatedGRMFD(
        a_val=3.0,
        b_val=1.0,
        min_mag=5.0,
        max_mag=6.0,
        bin_width=0.1,
    )
    src = PointSource(
        source_id="toy",
        name="toy",
        tectonic_region_type="Active Shallow Crust",
        mfd=mfd,
        rupture_mesh_spacing=1.0,
        magnitude_scaling_relationship=None,  # whatever you use
        rupture_aspect_ratio=1.0,
        upper_seismogenic_depth=0.0,
        lower_seismogenic_depth=15.0,
        location=Point(13.0, 43.0),
        nodal_plane_distribution=None,
        hypocenter_distribution=None,
    )
    return src

