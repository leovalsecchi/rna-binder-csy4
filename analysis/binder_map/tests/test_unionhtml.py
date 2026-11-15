import base64
import gzip

import numpy as np
import pandas as pd
import pytest

from binder_map.htmlbuild import encode_floats, encode_ints
from binder_map.unionhtml import (bake, display_rows, encode_stacks,
                                  join_layers, layer_mask, pack_medoids,
                                  stack_of, union_payload)


# --- the encoders the union page reuses from htmlbuild ---------------------
# unionhtml does NOT reimplement these; these tests pin the behaviour the
# page's Float32Array / Uint8Array views depend on.

def test_encode_floats_is_float32_not_float64():
    # 4 bytes per value: a float64 payload would be twice the size
    assert len(base64.b64decode(encode_floats([1.0, 2.0]))) == 8


def test_encode_floats_roundtrips_little_endian():
    raw = np.frombuffer(base64.b64decode(encode_floats([1.5, -2.25])), dtype="<f4")
    assert raw.tolist() == [1.5, -2.25]


def test_encode_ints_uint8_roundtrips():
    raw = np.frombuffer(base64.b64decode(encode_ints([0, 7, 255], "uint8")),
                        dtype=np.uint8)
    assert raw.tolist() == [0, 7, 255]


# --- the layer bitmask, which is what is new here -------------------------

def test_layer_mask_sets_one_bit_per_layer():
    flags = pd.DataFrame({"in_sample500k": [True, False, True],
                          "in_passing60k": [False, True, True],
                          "in_12k": [False, False, True]})
    assert layer_mask(flags).tolist() == [1, 2, 7]


def test_layer_mask_conserved_nodes_have_multiple_bits():
    flags = pd.DataFrame({"in_sample500k": [True], "in_passing60k": [True], "in_12k": [False]})
    mask = layer_mask(flags)
    assert bin(int(mask[0])).count("1") == 2


# --- the join: three tables, three different keys -------------------------
# meta/layers.parquet keys on coord_hash, everything else on node_id. An inner
# join that silently drops rows is the failure mode that makes the page wrong
# without making it look wrong, so join_layers refuses to lose one.

def _tables():
    nodes = pd.DataFrame({"node_id": [0, 1], "struct_name": ["d1", "d2"],
                          "n_designs": [1, 3]})
    canonical = pd.DataFrame({"coord_hash": ["h1", "h2"], "struct_name": ["d1", "d2"]})
    layers = pd.DataFrame({"coord_hash": ["h2", "h1"],
                           "in_sample500k": [True, True],
                           "in_passing60k": [True, False],
                           "in_12k": [False, False]})
    return nodes, canonical, layers


def test_join_layers_keeps_every_node_and_its_own_flags():
    joined = join_layers(*_tables())
    assert len(joined) == 2
    assert joined["node_id"].tolist() == [0, 1]
    # layers/ arrived in the other order: the flags must follow coord_hash,
    # not row position.
    assert joined["in_passing60k"].tolist() == [False, True]


def test_join_layers_raises_when_a_node_has_no_layer_row():
    nodes, canonical, layers = _tables()
    with pytest.raises(ValueError, match="no layer row"):
        join_layers(nodes, canonical, layers.iloc[:1])


def test_join_layers_raises_when_a_node_has_no_coord_hash():
    nodes, canonical, layers = _tables()
    with pytest.raises(ValueError, match="no coord_hash"):
        join_layers(nodes, canonical.iloc[:1], layers)


def test_join_layers_raises_when_a_coord_hash_is_duplicated():
    # a duplicated key multiplies rows instead of dropping them; the node
    # count has to be defended in both directions.
    nodes, canonical, layers = _tables()
    doubled = pd.concat([layers, layers.iloc[:1]], ignore_index=True)
    with pytest.raises(ValueError, match="row count"):
        join_layers(nodes, canonical, doubled)


# --- design id -> stack ---------------------------------------------------

def test_stack_of_splits_the_design_id_prefix():
    ids = pd.Series(["csy4free_HD_LiMPNN_len80_d045f606_recycle7_sample7_pack1",
                     "csy4free_BG_BIF_len80_c55bf73f_design16_seq8"])
    assert stack_of(ids).tolist() == ["HD+LiMPNN", "BG+BIF"]


def test_stack_of_rejects_an_unparseable_id():
    with pytest.raises(ValueError, match="design id"):
        stack_of(pd.Series(["not_a_design_id"]))


# --- payload assembly -----------------------------------------------------

def _frame():
    """Four nodes: campaign-only, campaign-only, passing-only, all three."""
    return pd.DataFrame({
        "node_id": [0, 1, 2, 3],
        "struct_name": [f"csy4free_HD_LiMPNN_len80_000000{i}_recycle1_sample1_pack1"
                        for i in range(4)],
        "n_designs": [1, 1, 2, 5],
        "x": [1.5, -2.25, 0.0, 4.0], "y": [0.0, 3.75, -1.0, 2.0],
        "cluster_r80": [0, 1, 1, 2],
        "in_sample500k": [True, True, False, True],
        "in_passing60k": [False, False, True, True],
        "in_12k": [False, False, False, True],
    })


def test_union_payload_interleaves_x_and_y_as_float32():
    payload = union_payload(_frame(), {"any": 1})
    raw = np.frombuffer(base64.b64decode(payload["coords"]), dtype="<f4")
    assert raw.tolist() == [1.5, 0.0, -2.25, 3.75, 0.0, -1.0, 4.0, 2.0]


def test_union_payload_mask_is_one_uint8_per_node():
    payload = union_payload(_frame(), {"any": 1})
    raw = np.frombuffer(base64.b64decode(payload["mask"]), dtype=np.uint8)
    assert raw.tolist() == [1, 1, 2, 7]
    assert payload["n"] == 4


def test_union_payload_counts_every_layer_it_drew():
    payload = union_payload(_frame(), {"any": 1})
    assert payload["stats"]["layers"] == {"sample": 3, "passing": 2, "library": 1}
    assert payload["stats"]["designs"] == 9
    assert payload["stats"]["decimation"] == 1


# --- decimation: thin the campaign, never the populations being compared ---

def test_display_rows_keeps_everything_at_factor_one():
    assert display_rows(_frame(), 1).all()


def test_display_rows_thins_only_campaign_only_nodes():
    keep = display_rows(_frame(), 2)
    # nodes 0,1 are campaign-only -> one of the two survives;
    # node 2 (passing-only) and node 3 (conserved) must both survive.
    assert keep.tolist() == [True, False, True, True]


def test_display_rows_never_drops_a_conserved_node():
    frame = _frame()
    for factor in (2, 3, 8, 32):
        keep = display_rows(frame, factor)
        conserved = layer_mask(frame) == 7
        assert keep[conserved].all(), factor
        assert keep[layer_mask(frame) == 2].all(), factor


# --- stack codes: a five-value category does not need two bytes ------------
# The union map draws 526,467 nodes; htmlbuild's int16 codes would spend half
# a megabyte of the page on a high byte that is always zero.

def test_encode_stacks_is_one_byte_per_node_and_sorted():
    lookup, codes = encode_stacks(pd.Series(["B+X", "A+Y", "B+X"]))
    assert lookup == ["A+Y", "B+X"]
    raw = np.frombuffer(base64.b64decode(codes), dtype=np.uint8)
    assert raw.tolist() == [1, 0, 1]


def test_union_payload_stack_codes_index_the_lookup():
    payload = union_payload(_frame(), {"any": 1})
    raw = np.frombuffer(base64.b64decode(payload["stack_codes"]), dtype=np.uint8)
    assert len(raw) == 4
    assert [payload["stacks"][code] for code in raw] == ["HD+LiMPNN"] * 4


# --- medoid CA traces -----------------------------------------------------
# _frame() has clusters [0, 1, 1, 2] on nodes [0, 1, 2, 3]. One medoid per
# cluster, plus a coarser resolution that must be ignored entirely.

def _name(i):
    return f"csy4free_HD_LiMPNN_len80_000000{i}_recycle1_sample1_pack1"


def _line(axis, nres=5, step=3.8):
    """A straight CA trace along one axis: tells the three medoids apart."""
    out = []
    for r in range(nres):
        point = [0.0, 0.0, 0.0]
        point[axis] = r * step
        out.append(point)
    return out


def _medoid_table():
    return pd.DataFrame({
        "resolution": ["cluster_r80"] * 3 + ["cluster_r10"] * 2,
        # deliberately out of order: pack_medoids sorts, it does not trust
        "cluster": [2, 0, 1, 0, 1],
        "node_id": [3, 0, 2, 1, 3],
        "struct_name": [_name(3), _name(0), _name(2), _name(1), _name(3)],
        "size": [1, 1, 2, 3, 1],
    })


def _trace_bundle():
    return {_name(0): _line(0), _name(1): _line(0),
            _name(2): _line(1), _name(3): _line(2)}


def _decode(payload):
    """Undo the page's own decode: gunzip, reshape, cumulative sum, rescale."""
    raw = np.frombuffer(gzip.decompress(base64.b64decode(payload["medoid_ca_gz"])),
                        dtype="<i2")
    delta = raw.reshape(payload["medoid_n"], payload["medoid_nres"], 3).astype(np.int64)
    return np.cumsum(delta, axis=1) / payload["medoid_scale"]


def _payload_with_medoids():
    return union_payload(_frame(), {"any": 1}, 1, _medoid_table(), _trace_bundle())


def test_pack_medoids_roundtrips_every_trace_onto_the_grid():
    payload = _payload_with_medoids()
    got = _decode(payload)
    assert got.shape == (3, 5, 3)
    for cluster, axis in enumerate((0, 1, 2)):
        want = np.asarray(_line(axis))
        # stored about each trace's own centroid, so compare centred
        centred = want - want.mean(axis=0)
        assert np.abs(got[cluster] - centred).max() <= 1 / payload["medoid_scale"]


def test_pack_medoids_indexes_by_cluster_id_not_row_order():
    # the table lists cluster 2 first; if pack_medoids trusted row order,
    # cluster 0 would get the z-axis trace instead of the x-axis one.
    got = _decode(_payload_with_medoids())
    spread = np.abs(got - got.mean(axis=1, keepdims=True)).sum(axis=1)
    assert spread[0].argmax() == 0
    assert spread[1].argmax() == 1
    assert spread[2].argmax() == 2


def test_pack_medoids_row_points_at_the_drawn_node():
    payload = _payload_with_medoids()
    rows = np.frombuffer(base64.b64decode(payload["medoid_row"]), dtype="<i4")
    # cluster 0 -> node 0 -> view row 0; cluster 1 -> node 2 -> row 2
    assert rows.tolist() == [0, 2, 3]
    sizes = np.frombuffer(base64.b64decode(payload["medoid_size"]), dtype="<i4")
    assert sizes.tolist() == [1, 2, 1]


def test_pack_medoids_row_is_minus_one_when_the_medoid_is_not_drawn():
    # decimation 2 drops node 1 (campaign-only). Rebuild the medoid table so
    # cluster 0's medoid IS that dropped node: the page must be told there is
    # no drawn row rather than handed a stale index into a shorter view.
    table = _medoid_table()
    table.loc[table.index[1], "node_id"] = 1            # cluster 0 -> node 1
    payload = union_payload(_frame(), {"any": 1}, 2, table, _trace_bundle())
    rows = np.frombuffer(base64.b64decode(payload["medoid_row"]), dtype="<i4")
    assert rows[0] == -1


def test_pack_medoids_rejects_cluster_ids_that_are_not_dense():
    table = _medoid_table()
    table.loc[table.index[0], "cluster"] = 7            # 0, 1, 7
    with pytest.raises(ValueError, match="dense"):
        pack_medoids(table, _trace_bundle(), {})


def test_pack_medoids_rejects_a_medoid_with_no_trace():
    bundle = _trace_bundle()
    del bundle[_name(2)]
    with pytest.raises(ValueError, match="no CA trace"):
        pack_medoids(_medoid_table(), bundle, {})


def test_pack_medoids_rejects_ragged_trace_lengths():
    bundle = _trace_bundle()
    bundle[_name(2)] = _line(1, nres=4)
    with pytest.raises(ValueError, match="same length"):
        pack_medoids(_medoid_table(), bundle, {})


def test_union_payload_without_medoids_carries_no_trace_keys():
    payload = union_payload(_frame(), {"any": 1})
    assert "medoid_ca_gz" not in payload
    assert "medoids" not in payload["stats"]


def test_union_payload_counts_the_medoids_it_packed():
    payload = _payload_with_medoids()
    assert payload["stats"]["medoids"] == 3
    assert payload["medoid_nres"] == 5


def test_union_payload_refuses_a_drawn_cluster_with_no_medoid():
    # node 3 sits in cluster 2; drop cluster 2's medoid and the page would
    # index past the end of the trace array and draw nothing, silently.
    table = _medoid_table()
    table = table[~((table["resolution"] == "cluster_r80") & (table["cluster"] == 2))]
    with pytest.raises(ValueError, match="only 2 medoids"):
        union_payload(_frame(), {"any": 1}, 1, table, _trace_bundle())


# --- baking: 3Dmol is optional -------------------------------------------
# Without medoids there is nothing for the viewer to draw, and half a megabyte
# of 3Dmol would be dead weight on a page that is already at the size limit.

def _template(tmp_path):
    path = tmp_path / "t.html"
    path.write_text("<!doctype html><script>/*__FFLATE__*/</script>"
                    "<script>/*__3DMOL__*/</script>"
                    "<script>window.UNION = __PAYLOAD__;</script></html>")
    return path


def test_bake_injects_3dmol_when_it_is_given(tmp_path):
    fflate = tmp_path / "f.js"; fflate.write_text("var FFLATE=1;")
    threedmol = tmp_path / "m.js"; threedmol.write_text("var $3Dmol=1;")
    out = tmp_path / "out.html"
    bake({"n": 0}, _template(tmp_path), fflate, out, threedmol)
    html = out.read_text()
    assert "var $3Dmol=1;" in html and "var FFLATE=1;" in html
    assert "__PAYLOAD__" not in html and "/*__3DMOL__*/" not in html


def test_bake_empties_the_3dmol_placeholder_when_it_is_not(tmp_path):
    fflate = tmp_path / "f.js"; fflate.write_text("var FFLATE=1;")
    out = tmp_path / "out.html"
    bake({"n": 0}, _template(tmp_path), fflate, out)
    html = out.read_text()
    assert "/*__3DMOL__*/" not in html
    assert "<script></script>" in html
