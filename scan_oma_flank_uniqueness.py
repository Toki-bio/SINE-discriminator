#!/usr/bin/env python3
"""Scan all oma top100/rand100 alignments; write TSV summary."""
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flank_uniqueness import scan

ALN = os.path.join(os.path.dirname(__file__), "alignments", "oma")


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "oma_flank_uniqueness.tsv"
    paths = sorted(glob.glob(os.path.join(ALN, "oma_*_top100.aln.fa")))
    paths += sorted(glob.glob(os.path.join(ALN, "oma_*_rand100.aln.fa")))
    rows = []
    for p in paths:
        r = scan(p)
        base = os.path.basename(p).replace(".aln.fa", "")
        for side in ("left", "right"):
            s = r.get(side, {})
            rows.append({
                "set": base,
                "tier": r.get("tier"),
                "side": side,
                "measured": s.get("measured"),
                "n_measured": s.get("n_measured"),
                "unique_frac": s.get("unique_frac"),
                "shared_copy_frac": s.get("shared_copy_frac"),
                "largest_cluster_frac": s.get("largest_cluster_frac"),
                "largest_cluster": s.get("largest_cluster"),
                "worst_flag": r.get("worst_flag"),
            })
    fields = list(rows[0].keys()) if rows else []
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    flagged = sum(1 for p in paths if scan(p).get("worst_flag"))
    print("wrote %s (%d side-rows, %d/%d sets flagged)" % (out, len(rows), flagged, len(paths)))


if __name__ == "__main__":
    main()
