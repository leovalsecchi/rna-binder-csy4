"""Derive inputs/screen_points.csv from the splicing report.

WHY THIS EXISTS: the per-measurement source table `matched_comparisons.csv`
(in the screen's 12_analysis folder) is owner-only (chmod 700, written by the
ONT pipeline) and could not be read. The same data is embedded in the
interactive report `splicing_report.html` as the JS objects ALL_POINTS,
NEG_CTRL_MAP and TRANSCRIPT_META, so we reconstruct the points from there using
the report's own fold-change logic (_computeFC, FC_CAP = 50).

VALIDATION: the 32 rows of gated_Q1_75Reads_OrderedHits.xlsx were re-derived and
matched exactly on reads (`total`), splicing ratio (`removed_fraction`) and
`matching_fold_change`. See MANIFEST.

USAGE (not part of run_all.sh; needs the ~48 MB report, which is not in the repo):
    python3 extract_from_report.py /path/to/splicing_report.html screen_points.csv

Columns written: transcript, comparison, reads, splicing_ratio,
matching_fold_change, log2fc, fold_change_raw, neg_reads.
neg_reads = read count of the matched NEG control (neg<k>_total).
matching_fold_change is capped at 50 like the report's export (so it matches the
xlsx); fold_change_raw is the same ratio uncapped (POS splicing / matched NEG
splicing, NEG = 1/NEG reads when no NEG read was spliced; checked equal to
on/neg for every row with NEG > 0).
"""
import sys, re, json, csv, math

FC_CAP = 50.0


def extract_obj(s, varname):
    m = re.search(varname + r'\s*=\s*', s)
    i = s.index('{', m.end())
    depth = 0; instr = False; esc = False; q = ''
    for j in range(i, len(s)):
        c = s[j]
        if instr:
            if esc: esc = False
            elif c == '\\': esc = True
            elif c == q: instr = False
        else:
            if c in '"\'': instr = True; q = c
            elif c == '{': depth += 1
            elif c == '}':
                depth -= 1
                if depth == 0:
                    return json.loads(s[i:j+1])
    raise ValueError("unterminated object for " + varname)


def main(html_path, out_csv):
    s = open(html_path, encoding="utf-8", errors="replace").read()
    ap = extract_obj(s, "var ALL_POINTS")["scatter"]
    ncm = extract_obj(s, "var NEG_CTRL_MAP")["scatter"]
    tm = extract_obj(s, "var TRANSCRIPT_META")
    pt = {p["t"]: (p["x"], p["y"]) for p in ap}   # transcript -> (reads, splicing)

    def is_matching(pos, neg, entry, k):
        m = tm.get(pos)
        if m and m.get("matched_neg") is not None:
            return m["matched_neg"] == neg
        return entry.get(f"neg{k}_matching") is True

    def matching_ratio(pos, entry):
        for k in range(1, 10):
            neg = entry.get(f"neg{k}")
            if neg and is_matching(pos, neg, entry, k):
                return entry.get(f"ratio{k}"), entry.get(f"neg{k}_total")
        return None, None

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["transcript", "comparison", "reads", "splicing_ratio",
                    "matching_fold_change", "log2fc", "fold_change_raw", "neg_reads"])
        n = 0
        for pos, entry in ncm.items():
            if pos not in pt: continue
            reads, spl = pt[pos]
            r, neg_reads = matching_ratio(pos, entry)
            fc = min(r, FC_CAP) if r is not None else None
            log2fc = math.log2(fc) if fc and fc > 0 else ""
            w.writerow([pos, pos.split("::")[0], reads, spl,
                        fc if fc is not None else "", log2fc,
                        r if r is not None else "",
                        neg_reads if neg_reads is not None else ""])
            n += 1
    print(f"wrote {n} points to {out_csv}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
