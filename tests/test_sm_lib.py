import os
import tempfile

import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon

from hazard2csep import sm_lib


# -------------------------------------------------------------------------
# Fixtures / paths (mirror test_sources.py style)
# -------------------------------------------------------------------------
_SRC_FOLDERS = [
    'area_src',
    'complexfault_src',
    'simplefault_src',
    'point_src',
    'multipoint_src',
    'simplefault_areasource_src',
    'simplefault_multipoint_src',
    'complexfault_crop_src',
]


def _test_src_xml_paths():
    """Yield full paths to test_src.xml in each fixture folder."""
    for folder in _SRC_FOLDERS:
        yield os.path.join(folder, "test_src.xml")


# -------------------------------------------------------------------------
# cleaner_range
# -------------------------------------------------------------------------
def test_cleaner_range_basic():
    print("\n[sm_lib.cleaner_range] test_cleaner_range_basic")
    start, end, h = 4.7, 8.1, 0.2
    expected = np.arange(start, end + h / 2, h)
    actual = sm_lib.cleaner_range(start, end, h)

    print("  expected:", expected)
    print("  actual:  ", actual)

    np.testing.assert_allclose(expected, actual)


def test_cleaner_range_single_value():
    print("\n[sm_lib.cleaner_range] test_cleaner_range_single_value")
    start = end = 5.0
    expected = np.array([5.0])
    actual = sm_lib.cleaner_range(start, end, 0.2)

    print("  expected:", expected)
    print("  actual:  ", actual)

    np.testing.assert_allclose(expected, actual)


def test_cleaner_range_non_integer_spacing():
    print("\n[sm_lib.cleaner_range] test_cleaner_range_non_integer_spacing")
    start, end, h = 0.0, 1.0, 0.3
    actual = sm_lib.cleaner_range(start, end, h)

    print("  actual:", actual)
    # 0.0, 0.3, 0.6, 0.9
    assert len(actual) == 4
    assert actual[0] == 0.0
    assert actual[-1] < end  # strictly below 1.0


def test_cleaner_range_reversed():
    print("\n[sm_lib.cleaner_range] test_cleaner_range_reversed")
    actual = sm_lib.cleaner_range(5.0, 4.0, 0.1)
    print("  actual:", actual)
    # numpy.arange with start > stop and positive step → empty array
    assert actual.size == 0


# -------------------------------------------------------------------------
# parse_source_model + parse_srcs on REAL XMLs
# -------------------------------------------------------------------------
def test_parse_source_model_and_srcs_real_xml():
    print("\n[sm_lib.parse_source_model + sm_lib.parse_srcs] test_parse_source_model_and_srcs_real_xml")

    for path in _test_src_xml_paths():
        print(f"  -> source model file: {path}")
        models = sm_lib.parse_source_model(path)
        print("    models returned:", len(models))
        assert len(models) >= 1

        sources = sm_lib.parse_srcs(models)
        print("    sources returned:", len(sources))
        assert len(sources) >= 1

        classes = sorted({s.__class__.__name__ for s in sources})
        print("    source classes:", classes)
        assert len(classes) > 0


def test_parse_srcs_empty_list():
    print("\n[sm_lib.parse_srcs] test_parse_srcs_empty_list")
    result = sm_lib.parse_srcs([])
    print("  expected: []")
    print("  actual:  ", result)
    assert result == []


# -------------------------------------------------------------------------
# parse_fault_source_buffer (slightly higher-level but controlled)
# -------------------------------------------------------------------------
class DummySrc:
    def __init__(self, source_id, fallback_poly):
        self.source_id = source_id
        self.polygon = type("PolyWrapper", (), {"_polygon2d": fallback_poly})


def test_parse_fault_source_buffer_match_and_fallback(tmp_path):
    print("\n[sm_lib.parse_fault_source_buffer] test_parse_fault_source_buffer_match_and_fallback")

    # Build a small GeoDataFrame with one matching and one non-matching ID
    polys = [
        Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
        Polygon([(2, 2), (3, 2), (3, 3), (2, 3)]),
    ]
    gdf = gpd.GeoDataFrame({"IDFS": ["SRC_MATCH", "SRC_OTHER"], "geometry": polys})
    shp_path = tmp_path / "buffers.shp"
    gdf.to_file(shp_path)

    # Dummy sources: one with matching IDFS, one without
    fallback1 = Polygon([(9, 9), (10, 9), (10, 10), (9, 10)])
    fallback2 = Polygon([(19, 19), (20, 19), (20, 20), (19, 20)])
    src_match = DummySrc("SRC_MATCH", fallback1)
    src_nomatch = DummySrc("SRC_NOMATCH", fallback2)

    updated = sm_lib.parse_fault_source_buffer(
        [src_match, src_nomatch],
        buffer_shp=str(shp_path),
    )

    print("  expected match polygon:", polys[0])
    print("  actual match polygon:  ", updated[0].buffer_polygon)
    print("  expected fallback polygon:", fallback2)
    print("  actual fallback polygon:  ", updated[1].buffer_polygon)

    # First source should use shapefile geometry
    assert updated[0].buffer_polygon.equals(polys[0])
    # Second source should use its own polygon._polygon2d
    assert updated[1].buffer_polygon.equals(fallback2)


# -------------------------------------------------------------------------
# parse_logictree_files using realistic-style logic tree XML
# -------------------------------------------------------------------------
REALISTIC_LT_XML = """<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns:gml="http://www.opengis.net/gml"
      xmlns="http://openquake.org/xmlns/nrml/0.4">
    <logicTree logicTreeID="lt1">
        <logicTreeBranchSet branchSetID="bs1" uncertaintyType="sourceModel">
            <logicTreeBranch branchID="ABL_ML">
                <uncertaintyModel>

                     source_models/asm_v12e/asm_ver12e_winGT_fs017_lo_abgrs_maxmag_low.xml
                     source_models/asm_v12e/asm_ver12e_winGT_fs017_twingr.xml

                     source_models/deep_v12e/asm_deep_ver12e_winGT_fs017_lo_abgrs_maxmag_low.xml
                     source_models/volcanic_v12e/asm_volcanic_ver12e_winGT_fs017_lo_abgrs_maxmag_low.xml

                     source_models/interface_v12b/CaA_IF2222222_M40.xml
                     source_models/interface_v12b/CyA_IF2222222_M40.xml

                </uncertaintyModel>
                <uncertaintyWeight>0.030</uncertaintyWeight>
            </logicTreeBranch>
            <logicTreeBranch branchID="ABL_MC">
                <uncertaintyModel>

                     source_models/asm_v12e/asm_ver12e_winGT_fs017_lo_abgrs_maxmag_mid.xml
                     source_models/asm_v12e/asm_ver12e_winGT_fs017_twingr.xml

                     source_models/deep_v12e/asm_deep_ver12e_winGT_fs017_lo_abgrs_maxmag_mid.xml
                     source_models/volcanic_v12e/asm_volcanic_ver12e_winGT_fs017_lo_abgrs_maxmag_mid.xml

                     source_models/interface_v12b/CaA_IF2222222_M40.xml
                     source_models/interface_v12b/CyA_IF2222222_M40.xml

                </uncertaintyModel>
                <uncertaintyWeight>0.024</uncertaintyWeight>
            </logicTreeBranch>
        </logicTreeBranchSet>
    </logicTree>
</nrml>
"""


def test_parse_logictree_files_realistic(tmp_path):
    print("\n[sm_lib.parse_logictree_files] test_parse_logictree_files_realistic")

    lt_path = tmp_path / "logic_tree.xml"
    lt_path.write_text(REALISTIC_LT_XML, encoding="utf-8")

    root_dir, branches = sm_lib.parse_logictree_files(str(lt_path))

    print("  root_dir:", root_dir)
    print("  branches keys:", list(branches.keys()))
    print("  raw branches dict:", branches)

    assert root_dir == str(tmp_path)
    assert set(branches.keys()) == {"ABL_ML", "ABL_MC"}

    # Clean up whitespace and empty tokens
    clean_branches = {
        bid: [f.strip() for f in files if f.strip()]
        for bid, files in branches.items()
    }

    print("  cleaned branches dict:", clean_branches)

    # Basic sanity: each branch should have the interface file in the list
    assert "source_models/interface_v12b/CaA_IF2222222_M40.xml" in clean_branches["ABL_ML"]
    assert "source_models/interface_v12b/CyA_IF2222222_M40.xml" in clean_branches["ABL_ML"]
    assert "source_models/interface_v12b/CaA_IF2222222_M40.xml" in clean_branches["ABL_MC"]
    assert "source_models/interface_v12b/CyA_IF2222222_M40.xml" in clean_branches["ABL_MC"]
