#!/usr/bin/env python3
"""
harvest_designs.py  (runs on Euler)
===================================

Pulls a sampled, equal-N set of DESIGN CIFs out of the campaign archives, one
group per backbone design method, ready for extract_binders.py.

Four groups, because the backbone is what the density is about and the inverse
folder moves no atoms: RFD3, BoltzGen, BoltzDesign, HalluDesign.  BG+BIF is left
out on purpose - it shares BoltzGen's backbones, and its stored design file is a
Boltz re-prediction whose pose is wrong by ~4.4 A median.

Archive shapes (both handled):
    big    <dir>/midrun_<c>.tar.zst        -> midrun_<c>/designs/<uuid>.tar.zst -> *.cif
    loose  <dir>/midrun_<c>_<uuid>.tar.zst -> midrun_<c>_<uuid>/designs.tar.zst -> *.cif
A loose per-uuid tar is preferred when it exists: a few MB instead of streaming
30 GB.  Nothing is written inside the campaign tree.

    python3 harvest_designs.py --out <dir> --n 800 --dry-run   # plan only
    python3 harvest_designs.py --out <dir> --n 800
"""

import argparse
import collections
import csv
import io
import os
import random
import sys
import tarfile

import zstandard as zstd

RB = "/cluster/work/wicky/lvalsecchi/rna_binder"
CAMP = RB + "/campaigns/csy4_crrna_free"
PASSING = RB + "/analyses/foldseek_map/passing_v6_final.csv"

# group -> (stack label in the table, archive directory, campaigns or None=all)
# One campaign per group where the group is large enough, so that each group is
# one generation batch rather than a blend; HalluDesign has to pool to reach N,
# and its designs are loose per-uuid tars anyway, so pooling costs nothing.
GROUPS = {
    "RFD3":        ("RFD3+LiMPNN", "rfd3__ligandmpnn", ["v7"]),
    "BoltzGen":    ("BG+LiMPNN", "boltzgen__ligandmpnn", ["v8"]),
    "BoltzDesign": ("BD1+LiMPNN", "boltzdesign__ligandmpnn", ["v8"]),
    "HalluDesign": ("HD+LiMPNN", "halludesign__ligandmpnn", None),
}


def read_passing(path):
    """(stack, campaign, uuid8, design_id) for every row; '#' header is skipped."""
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                head = next(csv.reader([line]))
                break
        cols = {c: i for i, c in enumerate(head)}
        for row in csv.reader(fh):
            if not row:
                continue
            yield (row[cols["stack"]], row[cols["campaign"]],
                   row[cols["uuid8"]], row[cols["design_id"]])


def loose_has_designs(path):
    """Does this loose per-uuid tar actually hold designs?

    Many HalluDesign ones are 2 KB stubs carrying only
    designs/processing_results.csv, while the designs themselves are in the big
    campaign tar.  Treating "a loose tar exists" as "the designs are here" let a
    stub shadow the real archive and silently lost 341 of 800 designs.  The
    .list sidecar answers it without opening anything.
    """
    lst = path + ".list"
    if not os.path.exists(lst):
        return True                      # no index: try it and let the retry sort it out
    with open(lst) as fh:
        for line in fh:
            line = line.rstrip()
            if line.endswith("designs.tar.zst") or line.endswith(".cif"):
                return True
    return False


def big_path(archive_dir, campaign):
    return "%s/%s/midrun_%s.tar.zst" % (CAMP, archive_dir, campaign)


def sources(archive_dir, campaign, uuids):
    """[(kind, path, uuid(s))] - a loose tar per uuid where one exists AND holds
    designs, else one big tar carrying all the remaining uuids."""
    out, rest = [], []
    for u in sorted(uuids):
        loose = "%s/%s/midrun_%s_%s.tar.zst" % (CAMP, archive_dir, campaign, u)
        if os.path.exists(loose) and loose_has_designs(loose):
            out.append(("loose", loose, {u}))
        else:
            rest.append(u)
    if rest:
        out.append(("big", big_path(archive_dir, campaign), set(rest)))
    return out


def nested_cifs(blob, wanted, write):
    """Extract the wanted <design_id>.cif members from a nested .tar.zst blob."""
    got = set()
    with zstd.ZstdDecompressor().stream_reader(io.BytesIO(blob)) as r, \
            tarfile.open(fileobj=r, mode="r|") as tf:
        for m in tf:
            if not m.isfile():
                continue
            base = os.path.basename(m.name)
            if base.endswith(".cif") and base[:-4] in wanted:
                write(base[:-4], tf.extractfile(m).read())
                got.add(base[:-4])
                if got == wanted:
                    break
    return got


def harvest(kind, path, uuids, want_by_uuid, campaign, write, log):
    """One pass over one archive.  Returns the set of design_ids written."""
    got = set()
    if not os.path.exists(path):
        log("  MISSING %s" % path)
        return got
    wanted_members = {}
    if kind == "big":
        for u in uuids:
            wanted_members["midrun_%s/designs/%s.tar.zst" % (campaign, u)] = u
    else:
        u = next(iter(uuids))
        wanted_members["midrun_%s_%s/designs.tar.zst" % (campaign, u)] = u
    seen = set()
    with open(path, "rb") as fh, \
            zstd.ZstdDecompressor().stream_reader(fh) as rd, \
            tarfile.open(fileobj=rd, mode="r|") as tf:
        for m in tf:
            u = wanted_members.get(m.name)
            if u is None or not m.isfile():
                continue
            got |= nested_cifs(tf.extractfile(m).read(), want_by_uuid[u], write)
            seen.add(u)
            if seen == set(uuids):
                break
    return got


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=800,
                    help="designs per group; equal N keeps the percentage "
                         "scales comparable between groups")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--passing", default=PASSING)
    ap.add_argument("--dry-run", action="store_true",
                    help="resolve and print the archives, read none of them")
    ap.add_argument("--only", default=None,
                    help="one group name, so the four can run as an array task "
                         "each (the sample is seeded, so it does not shift)")
    args = ap.parse_args()

    rows = collections.defaultdict(list)
    for stack, campaign, uuid8, design_id in read_passing(args.passing):
        rows[stack].append((campaign, uuid8, design_id))

    rng = random.Random(args.seed)
    plan, picked = {}, {}
    for group, (stack, archive_dir, campaigns) in sorted(GROUPS.items()):
        pool = [r for r in rows.get(stack, ())
                if campaigns is None or r[0] in campaigns]
        if not pool:
            print("%-12s NO ROWS for stack %s" % (group, stack))
            continue
        take = pool if len(pool) <= args.n else rng.sample(pool, args.n)
        picked[group] = take
        by_src = collections.defaultdict(lambda: collections.defaultdict(set))
        for campaign, uuid8, design_id in take:
            by_src[campaign][uuid8].add(design_id)
        plan[group] = (archive_dir, by_src)
        print("%-12s %s: %d of %d available, %d uuids, campaigns %s"
              % (group, stack, len(take), len(pool),
                 sum(len(v) for v in by_src.values()),
                 ",".join(sorted(by_src))))

    # --only is applied AFTER sampling, so each array task draws the same sample
    # it would have drawn in a single run
    if args.only:
        if args.only not in GROUPS:
            ap.error("--only must be one of %s" % ", ".join(sorted(GROUPS)))
        plan = {k: v for k, v in plan.items() if k == args.only}

    total_big = 0
    for group, (archive_dir, by_src) in sorted(plan.items()):
        for campaign, want_by_uuid in sorted(by_src.items()):
            for kind, path, uuids in sources(archive_dir, campaign, want_by_uuid):
                size = os.path.getsize(path) / 1e9 if os.path.exists(path) else -1
                if kind == "big":
                    total_big += size
                    print("  %-12s %-5s %-5s %5.1f GB  %3d uuids  %s"
                          % (group, campaign, kind, size, len(uuids),
                             os.path.basename(path)))
        n_loose = sum(1 for c, w in by_src.items()
                      for k, _, _ in sources(archive_dir, c, w) if k == "loose")
        if n_loose:
            print("  %-12s %d loose per-uuid tars" % (group, n_loose))
    print("big-tar bytes to stream: %.1f GB" % total_big)
    if args.dry_run:
        return 0

    os.makedirs(args.out, exist_ok=True)
    manifest = []
    for group, (archive_dir, by_src) in sorted(plan.items()):
        gdir = os.path.join(args.out, group)
        os.makedirs(gdir, exist_ok=True)
        written = set()

        def write(design_id, data, gdir=gdir, written=written):
            with open(os.path.join(gdir, design_id + ".cif"), "wb") as fh:
                fh.write(data)
            written.add(design_id)

        for campaign, want_by_uuid in sorted(by_src.items()):
            srcs = sources(archive_dir, campaign, want_by_uuid)
            for kind, path, uuids in srcs:
                harvest(kind, path, uuids, want_by_uuid, campaign, write,
                        lambda m: print(m, flush=True))
            # a loose tar that existed but did not hold what we wanted: fall back
            # to this campaign's big tar for those uuids, in one extra pass
            tried_big = {u for k, _, us in srcs if k == "big" for u in us}
            retry = {u for u, ids in want_by_uuid.items()
                     if ids - written and u not in tried_big}
            if retry:
                print("  %-12s %-4s retry %d uuids in the big tar"
                      % (group, campaign, len(retry)), flush=True)
                harvest("big", big_path(archive_dir, campaign), retry,
                        want_by_uuid, campaign, write,
                        lambda m: print(m, flush=True))
            print("  %-12s %-4s -> %d files so far"
                  % (group, campaign, len(written)), flush=True)
        wanted = {d for w in by_src.values() for s in w.values() for d in s}
        print("%-12s %d/%d written%s"
              % (group, len(written), len(wanted),
                 "" if written == wanted else "  MISSING %d" % len(wanted - written)),
              flush=True)
        for design_id in sorted(written):
            manifest.append((os.path.join(gdir, design_id + ".cif"), group))

    name = "designs_%s.csv" % args.only if args.only else "designs.csv"
    with open(os.path.join(args.out, name), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["design_cif", "group"])
        w.writerows(manifest)
    print("wrote %d designs + %s to %s" % (len(manifest), name, args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
