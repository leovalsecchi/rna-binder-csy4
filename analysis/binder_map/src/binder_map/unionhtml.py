"""Assemble the union fold map's payload for the browser.

526,467 nodes cannot travel as JSON objects: one dict per node would put the
page in the hundreds of MB. Everything the page needs per node is a fixed-width
column, so it ships as the base64 typed arrays `htmlbuild` already encodes --
this module holds only what is specific to the union map and imports the rest.

What is genuinely new here is the layer bitmask: three populations were
co-embedded into one map, and a node belonging to more than one of them is
"conserved", which is the whole point of the exercise.

The page also carries ONE CA trace per union cluster -- its medoid -- so a
click can show what the fold under the cursor actually looks like. The v6
explorer shipped a trace for every node; that is what made it 24 MB, and it is
deliberately not repeated here. Medoids only, at the operative resolution only.
"""
from __future__ import annotations

import base64
import gzip
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from binder_map.htmlbuild import (INT16_MAX, encode_floats, encode_ints,
                                  encode_strings)

#: Flag columns of meta/layers.parquet, low bit first.
LAYER_COLUMNS = ("in_sample500k", "in_passing60k", "in_12k")

#: `csy4free_<gen>_<invfold>_len<N>_...`. The stack is the two tokens between
#: the target name and the length, written with a `+` everywhere else in this
#: project ("RFD3+LiMPNN"), so it is normalised here rather than at each call.
_STACK_RE = r"^csy4free_([A-Za-z0-9]+)_([A-Za-z0-9]+)_len\d+_"

#: The operative resolution. The coarser columns ride along in clusters.parquet
#: but r80 is the one the campaign reads.
RESOLUTION = "cluster_r80"

#: Grid for the medoid CA traces, in steps per Angstrom: 4 = a 0.25 A grid.
#:
#: The traces are the single biggest thing on the page, so the grid is a size
#: decision, not a precision one. 32,394 medoids x 80 residues at stage 7's
#: native 0.01 A is 18 MB gzipped -- alone more than the rest of the page. On a
#: 0.25 A grid, delta-encoded along the chain, the same traces are 6.2 MB
#: gzipped. A quarter of an Angstrom is an order of magnitude below what a
#: ~250 px CA ribbon can show, so nothing visible is lost; these are for
#: recognising a fold, and anything measured comes from the real structure.
TRACE_SCALE = 4

#: From the stage 6 run log (logs/union_embed_fold_15438409.log):
#:   "design/fold: embedding 526,467 nodes from 12,027,274 layout edges"
#:   "layout graph has 96 components -> init=random"
#: The page states both: with 96 components UMAP could not use spectral init,
#: so where one component sits relative to another is arbitrary.
LAYOUT_COMPONENTS = 96
LAYOUT_EDGES = 12_027_274


#: bit 0 = campaign sample, bit 1 = passing pool, bit 2 = ordered library.
#: A node with more than one bit set is conserved -- that is the whole signal.
def layer_mask(flags: pd.DataFrame) -> np.ndarray:
    return (flags["in_sample500k"].to_numpy(bool).astype(np.uint8)
            | (flags["in_passing60k"].to_numpy(bool).astype(np.uint8) << 1)
            | (flags["in_12k"].to_numpy(bool).astype(np.uint8) << 2))


def stack_of(struct_names: pd.Series) -> pd.Series:
    """Generation stack for each design id, e.g. `HD+LiMPNN`.

    Parsed from the id rather than looked up in a manifest: no single manifest
    covers all three source populations, and the id is what stage 2 elected as
    canonical, so it is the only key guaranteed present for every union node.
    An id that does not parse raises -- a silently-blank stack would show up as
    an unlabelled colour in the legend and nowhere else.
    """
    parts = struct_names.str.extract(_STACK_RE)
    bad = parts[0].isna()
    if bool(bad.any()):
        raise ValueError(
            f"{int(bad.sum())} design id(s) do not match "
            f"csy4free_<gen>_<invfold>_len<N>_, first: {struct_names[bad].iloc[0]!r}")
    return parts[0] + "+" + parts[1]


def encode_stacks(stacks: pd.Series) -> tuple[list[str], str]:
    """(sorted stack names, base64 uint8 codes).

    `htmlbuild.encode_categorical` is not reused: it emits int16 codes, and at
    526,467 nodes the high byte of a five-value category is half a megabyte of
    zeroes on the wire. Five stacks fit a byte with room to spare.
    """
    lookup = sorted(stacks.unique())
    if len(lookup) > 255:
        raise ValueError(f"too many stacks for uint8 codes: {len(lookup)}")
    index = {name: i for i, name in enumerate(lookup)}
    return lookup, encode_ints(stacks.map(index).to_numpy(), "uint8")


def join_layers(nodes: pd.DataFrame, canonical: pd.DataFrame,
                layers: pd.DataFrame) -> pd.DataFrame:
    """Attach each node's layer flags: node_id -> struct_name -> coord_hash.

    `meta/layers.parquet` keys on `coord_hash` while `dedup/nodes.tsv` keys on
    `node_id`, and `dedup/canonical.tsv` is the only bridge between them. The
    joins are LEFT, not inner, on purpose: an inner join answers a missing key
    by dropping the row, which would give a page that renders happily with
    fewer nodes than it claims. Left joins turn the same event into a NaN this
    function can refuse.
    """
    expected = len(nodes)
    bridged = nodes.merge(canonical[["struct_name", "coord_hash"]],
                          on="struct_name", how="left", validate="one_to_one")
    orphan = bridged["coord_hash"].isna()
    if bool(orphan.any()):
        raise ValueError(
            f"{int(orphan.sum())} node(s) have no coord_hash in canonical.tsv, "
            f"first: {bridged.loc[orphan, 'struct_name'].iloc[0]!r}")

    joined = bridged.merge(layers[["coord_hash", *LAYER_COLUMNS]],
                           on="coord_hash", how="left")
    if len(joined) != expected:
        raise ValueError(
            f"join changed the row count: {expected} nodes in, {len(joined)} out "
            "(a duplicated coord_hash in meta/layers.parquet multiplies rows)")
    missing = joined[list(LAYER_COLUMNS)].isna().any(axis=1)
    if bool(missing.any()):
        raise ValueError(
            f"{int(missing.sum())} node(s) have no layer row in layers.parquet, "
            f"first coord_hash: {joined.loc[missing, 'coord_hash'].iloc[0]!r}")
    return joined


def display_rows(frame: pd.DataFrame, decimation: int) -> np.ndarray:
    """Boolean keep-mask for `decimation`. Thins campaign-ONLY nodes and
    nothing else: the passing pool, the library and every conserved backbone
    stay whole, because they are the populations the page exists to compare.
    """
    keep = np.ones(len(frame), dtype=bool)
    if decimation <= 1:
        return keep
    campaign_only = layer_mask(frame) == 1
    where = np.flatnonzero(campaign_only)
    keep[where] = False
    keep[where[::decimation]] = True
    return keep


def pack_medoids(medoids: pd.DataFrame, traces: dict[str, list],
                 node_row: dict[int, int]) -> dict:
    """One CA trace per union cluster, indexed BY cluster id.

    The traces are stored in cluster order and nothing else, so the page looks
    one up with a multiply instead of carrying a 32k-entry cluster->name map it
    would otherwise need. That only works if the cluster ids really are a dense
    0..K-1 range with exactly one medoid each, which is asserted rather than
    assumed: an off-by-one here shows the WRONG backbone for a cluster, and a
    wrong backbone looks exactly as convincing as a right one.

    Coordinates are quantised onto the `TRACE_SCALE` grid and then delta-encoded
    along the chain. Consecutive CA atoms sit ~3.8 A apart, so the steps are
    small numbers that deflate compresses far better than absolute coordinates
    (6.2 MB against 9.4 MB, measured on this tree).
    """
    frame = medoids[medoids["resolution"] == RESOLUTION].sort_values("cluster")
    ids = frame["cluster"].to_numpy()
    if len(ids) == 0:
        raise ValueError(f"medoids.parquet has no {RESOLUTION} rows")
    if ids[0] != 0 or not bool((np.diff(ids) == 1).all()):
        raise ValueError(
            f"{RESOLUTION} medoid cluster ids are not a dense 0..{len(ids) - 1} "
            f"range (min {ids.min()}, max {ids.max()}, {len(ids)} rows): the "
            "page indexes traces by cluster id and would show the wrong one")

    names = frame["struct_name"].tolist()
    absent = [name for name in names if name not in traces]
    if absent:
        raise ValueError(f"{len(absent)} medoid(s) have no CA trace in the "
                         f"stage 7 bundle, first: {absent[0]!r}")
    lengths = {len(traces[name]) for name in names}
    if len(lengths) != 1:
        raise ValueError(
            "medoid CA traces are not all the same length "
            f"({sorted(lengths)[:4]}...); the page slices one flat array at a "
            "fixed stride")
    nres = lengths.pop()

    grid = np.rint(np.asarray([traces[name] for name in names],
                              dtype=np.float64) * TRACE_SCALE).astype(np.int32)
    grid -= np.rint(grid.mean(axis=1, keepdims=True)).astype(np.int32)
    delta = grid.copy()
    delta[:, 1:, :] = grid[:, 1:, :] - grid[:, :-1, :]
    if delta.min() < -INT16_MAX - 1 or delta.max() > INT16_MAX:
        raise ValueError(f"a CA step overflows int16 at scale {TRACE_SCALE}: "
                         f"range {delta.min()}..{delta.max()}")

    rows = np.array([node_row.get(int(node), -1) for node in frame["node_id"]],
                    dtype=np.int32)
    blob = gzip.compress(np.ascontiguousarray(delta, dtype="<i2").tobytes(), 9)
    return {
        "medoid_ca_gz": base64.b64encode(blob).decode("ascii"),
        "medoid_nres": int(nres),
        "medoid_scale": int(TRACE_SCALE),
        "medoid_row": encode_ints(rows, "int32"),
        "medoid_size": encode_ints(frame["size"].to_numpy(), "int32"),
        "medoid_n": int(len(ids)),
    }


def union_payload(frame: pd.DataFrame, conserved: dict, decimation: int = 1,
                  medoids: pd.DataFrame | None = None,
                  traces: dict[str, list] | None = None) -> dict:
    """Pack one row per drawn backbone into the page's payload."""
    view = frame.loc[display_rows(frame, decimation)]

    # Interleaved [x0,y0,x1,y1,...]: the renderer walks nodes in order, and one
    # strided read beats two. encode_floats does the float32 + base64.
    xy = np.empty(2 * len(view), dtype="<f4")
    xy[0::2] = view["x"].to_numpy()
    xy[1::2] = view["y"].to_numpy()

    clusters = view[RESOLUTION].to_numpy()
    if clusters.max() > INT16_MAX:
        raise ValueError(f"cluster id {clusters.max()} overflows the int16 codes")
    stack_lookup, stack_codes = encode_stacks(stack_of(view["struct_name"]))

    payload = {
        "n": int(len(view)),
        "coords": encode_floats(xy),
        "mask": encode_ints(layer_mask(view), "uint8"),
        "cluster": encode_ints(clusters, "int16"),
        "stack_codes": stack_codes,
        "stacks": stack_lookup,
        "ids_gz": encode_strings(view["struct_name"].tolist()),
        "resolution": RESOLUTION,
        "built": date.today().isoformat(),
        "stats": {
            "nodes": int(len(view)),
            "designs": int(view["n_designs"].sum()),
            "layers": {"sample": int(view["in_sample500k"].sum()),
                       "passing": int(view["in_passing60k"].sum()),
                       "library": int(view["in_12k"].sum())},
            "conserved": conserved,
            "clusters": int(view[RESOLUTION].nunique()),
            "singletons": int((view[RESOLUTION].value_counts() == 1).sum()),
            "components": LAYOUT_COMPONENTS,
            "edges": LAYOUT_EDGES,
            "decimation": int(decimation),
        },
    }

    if medoids is not None and traces is not None:
        node_row = {int(node): row
                    for row, node in enumerate(view["node_id"].to_numpy())}
        packed = pack_medoids(medoids, traces, node_row)
        if clusters.max() >= packed["medoid_n"]:
            raise ValueError(
                f"cluster {clusters.max()} is drawn but only "
                f"{packed['medoid_n']} medoids were packed")
        payload.update(packed)
        payload["stats"]["medoids"] = packed["medoid_n"]
    return payload


def bake(payload: dict, template: Path, fflate: Path, out: Path,
         threedmol: Path | None = None) -> int:
    """Substitute the payload and the vendored libraries into the template.

    `htmlbuild.bake` is not reused here: it injects a separate app.js and
    style.css the union page does not have, and hard-wires 3Dmol in rather than
    leaving it optional. 3Dmol IS optional here -- without the medoid bundle
    there is nothing for it to draw, and half a megabyte of viewer would be
    dead weight -- so the placeholder is emptied rather than filled.
    """
    import json

    html = template.read_text()
    html = html.replace("/*__FFLATE__*/", fflate.read_text())
    html = html.replace("/*__3DMOL__*/",
                        threedmol.read_text() if threedmol is not None else "")
    html = html.replace("__PAYLOAD__",
                        json.dumps(payload, separators=(",", ":"), allow_nan=False))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html)
    return out.stat().st_size
