#!/usr/bin/env python3
"""Where does each COPY stop being the element, and do the copies agree?

Every boundary measure so far has been anchored on the consensus: it asks how
identity behaves at the column where the consensus happens to end. That assumes
the consensus edge is the element edge, which is exactly what is in doubt when
Sergei says "not firm edges" (hyd_SINE_6) or "weak right end" (hyd_SINE_5).

This asks the copies instead. For each copy, slide a 15 bp window across the
element and 150 bp of flank either side, and find the outermost position on each
side where the copy still matches the consensus above chance. That position is
where THAT copy stops being the element.

Then the two numbers that matter:

  spread5 / spread3   how much the copies disagree about where the element ends
                      (interquartile range of those positions, in bp)
  offset5 / offset3   median distance from the consensus edge to where the
                      copies actually stop - negative means the copies stop
                      short, positive means they carry on past it

A real family has copies that end together, at the consensus edge: small spread,
offset near zero. A ragged edge shows as a large spread whatever the mean does,
and that is a different fault from the consensus simply being the wrong length,
which shows as a large offset with a small spread.
"""
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np

GAP = "-"
WIN = 15
OUT = 150
HIT = 0.60          # a window this concordant still counts as element
MIN_COPIES = 12


def read_fa(p):
    names, seqs, cur, buf = [], [], None, []
    for line in open(p):
        line = line.rstrip("\n\r")
        if line.startswith(">"):
            if cur is not None:
                seqs.append("".join(buf))
            cur = line[1:]
            names.append(cur)
            buf = []
        else:
            buf.append(line.strip())
    if cur is not None:
        seqs.append("".join(buf))
    return names, seqs


def measure(path):
    names, seqs = read_fa(path)
    ci = [i for i, n in enumerate(names) if "CONSENSUS_" in n.upper()]
    if not ci or len(seqs) < MIN_COPIES:
        return None
    k = ci[0]
    cons = seqs[k].upper()
    nz = [i for i, c in enumerate(cons) if c != GAP]
    if len(nz) < 60:
        return None
    lo, hi = nz[0], nz[-1]

    # the consensus has no sequence outside the element, so a copy's flank can
    # only be compared to other copies - use the column-majority base there
    rows = [s.upper() for i, s in enumerate(seqs) if i != k]
    A = np.array([list(r) for r in rows])
    ref = np.array(list(cons))
    a, b = max(0, lo - OUT), min(A.shape[1], hi + 1 + OUT)
    for j in list(range(a, lo)) + list(range(hi + 1, b)):
        col = A[:, j]
        col = col[np.isin(col, list("ACGT"))]
        if len(col) >= 8:
            v, c = np.unique(col, return_counts=True)
            ref[j] = v[np.argmax(c)]

    sub = A[:, a:b]
    rsub = ref[a:b]
    el_lo, el_hi = lo - a, hi - a
    W = sub.shape[1]
    if W < 3 * WIN:
        return None

    match = (sub == rsub[None, :]) & (sub != GAP) & (rsub != GAP)
    valid = (sub != GAP) & (rsub != GAP)
    # windowed concordance per copy
    ker = np.ones(WIN)
    left_pos, right_pos = [], []
    for i in range(sub.shape[0]):
        m = np.convolve(match[i].astype(float), ker, "same")
        v = np.convolve(valid[i].astype(float), ker, "same")
        with np.errstate(invalid="ignore", divide="ignore"):
            rate = np.where(v >= WIN * 0.5, m / np.maximum(v, 1), np.nan)
        good = np.nan_to_num(rate, nan=0.0) >= HIT
        if not good[el_lo:el_hi + 1].any():
            continue
        # walk outward from the middle of the element to the last good window
        mid = (el_lo + el_hi) // 2
        i5 = mid
        while i5 - 1 >= 0 and good[i5 - 1]:
            i5 -= 1
        i3 = mid
        while i3 + 1 < W and good[i3 + 1]:
            i3 += 1
        left_pos.append(i5 - el_lo)      # negative: starts before the consensus
        right_pos.append(i3 - el_hi)     # positive: runs past the consensus
    if len(left_pos) < MIN_COPIES:
        return None

    def iqr(x):
        return float(np.percentile(x, 75) - np.percentile(x, 25))

    return {"spread5": round(iqr(left_pos), 1), "spread3": round(iqr(right_pos), 1),
            "offset5": round(float(np.median(left_pos)), 1),
            "offset3": round(float(np.median(right_pos)), 1),
            "n": len(left_pos)}


def _one(f):
    try:
        return os.path.basename(f).replace(".aln.fa", ""), measure(f)
    except Exception as exc:
        return os.path.basename(f).replace(".aln.fa", ""), {"error": str(exc)}


def main():
    src, out = sys.argv[1], sys.argv[2]
    files = sorted(glob.glob(os.path.join(src, "*.aln.fa")))
    import multiprocessing as mp
    res = {}
    with mp.Pool(max(1, mp.cpu_count() - 4)) as pool:
        for name, r in pool.imap_unordered(_one, files, chunksize=4):
            if r and "error" not in r:
                res[name] = r
    json.dump(res, open(out, "w"), separators=(",", ":"))
    print("measured %d of %d" % (len(res), len(files)))
    print("\n%-28s %8s %8s %8s %8s" % ("set", "spread5", "spread3", "off5", "off3"))
    for k in sorted(res):
        if "top100" not in k:
            continue
        v = res[k]
        print("%-28s %8.1f %8.1f %8.1f %8.1f"
              % (k.replace("__top100", "")[:28], v["spread5"], v["spread3"],
                 v["offset5"], v["offset3"]))


if __name__ == "__main__":
    main()
