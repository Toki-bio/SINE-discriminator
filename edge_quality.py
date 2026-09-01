#!/usr/bin/env python3
"""Are the element's ends firm?

Sergei uses this constantly and the tool has never measured it:
  hyd_SINE_5  "weak right end, mosaic left end"
  hyd_SINE_6  "not firm edges"
  aca_SINE_1  "almost no boundary"

A real SINE ends. Copies agree with each other right up to the last element
column and are at background one column later. A weak end declines gradually:
the copies stop agreeing well before the consensus says the element stops, so
the boundary the consensus claims is not the boundary the copies have.

Measured per column as pairwise identity among the copies present, then:

  edge5 / edge3   mean identity over the first / last 25 element columns
  mid             mean over the middle of the element
  step5 / step3   how far identity falls in the 25 flank columns just outside,
                  i.e. how sharp the boundary is
  ratio5 / ratio3 edge / mid - below 1 means the end is weaker than the body

The 3' end needs the composition guard already established: identity falling
while A+T rises is a poly-A tail, which is part of a SINE, not a fault. That is
reported separately (at3) rather than silently excused.
"""
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np

GAP = "-"
W = 25
MIN_COPIES = 8


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


def col_identity(block):
    """Pairwise identity per column, NaN where too few copies are present."""
    out = []
    for j in range(block.shape[1]):
        col = block[:, j]
        col = col[np.isin(col, list("ACGT"))]
        n = len(col)
        if n < MIN_COPIES:
            out.append(np.nan)
            continue
        _, cnt = np.unique(col, return_counts=True)
        out.append(float((cnt * (cnt - 1)).sum()) / (n * (n - 1)))
    return np.asarray(out)


def at_fraction(block):
    out = []
    for j in range(block.shape[1]):
        col = block[:, j]
        col = col[np.isin(col, list("ACGT"))]
        if len(col) < MIN_COPIES:
            out.append(np.nan)
            continue
        out.append(float(np.isin(col, list("AT")).sum()) / len(col))
    return np.asarray(out)


def measure(path):
    names, seqs = read_fa(path)
    ci = [i for i, n in enumerate(names) if "CONSENSUS_" in n.upper()]
    if not ci or len(seqs) < 12:
        return None
    k = ci[0]
    cons = seqs[k]
    nz = [i for i, c in enumerate(cons) if c != GAP]
    if len(nz) < 60:
        return None
    lo, hi = nz[0], nz[-1]
    rows = [s.upper() for i, s in enumerate(seqs) if i != k]
    A = np.array([list(r) for r in rows])

    el = A[:, lo:hi + 1]
    idc = col_identity(el)
    atc = at_fraction(el)
    if np.isnan(idc).all():
        return None
    L = len(idc)
    if L < 3 * W:
        return None

    def m(x):
        x = x[np.isfinite(x)]
        return float(np.mean(x)) if len(x) else np.nan

    edge5, edge3 = m(idc[:W]), m(idc[-W:])
    mid = m(idc[W:-W])
    at3, atmid = m(atc[-W:]), m(atc[W:-W])

    # the flank just outside, to see how far identity actually drops
    outL = A[:, max(0, lo - W):lo]
    outR = A[:, hi + 1:hi + 1 + W]
    out5 = m(col_identity(outL)) if outL.shape[1] >= 5 else np.nan
    out3 = m(col_identity(outR)) if outR.shape[1] >= 5 else np.nan

    def r(x):
        return None if not np.isfinite(x) else round(float(x), 3)

    return {"edge5": r(edge5), "edge3": r(edge3), "mid": r(mid),
            "ratio5": r(edge5 / mid) if np.isfinite(mid) and mid > 0 else None,
            "ratio3": r(edge3 / mid) if np.isfinite(mid) and mid > 0 else None,
            "step5": r(edge5 - out5), "step3": r(edge3 - out3),
            "at3": r(at3), "at_mid": r(atmid),
            "tail3": bool(np.isfinite(at3) and np.isfinite(atmid) and at3 - atmid > 0.15),
            "elem_cols": int(L)}


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

    by = defaultdict(list)
    for k, v in res.items():
        by[k.split("__")[0]].append(v)
    print("\n%-16s %5s %8s %8s %8s %8s %8s"
          % ("group", "n", "ratio5", "ratio3", "step5", "step3", "tail3"))
    print("-" * 66)
    g = lambda v, key: np.median([x[key] for x in v if x.get(key) is not None]) \
        if any(x.get(key) is not None for x in v) else float("nan")
    for c in sorted(by, key=lambda c: g(by[c], "ratio3")):
        v = by[c]
        print("%-16s %5d %8.3f %8.3f %8.3f %8.3f %8s"
              % (c[:16], len(v), g(v, "ratio5"), g(v, "ratio3"),
                 g(v, "step5"), g(v, "step3"),
                 sum(1 for x in v if x.get("tail3"))))


if __name__ == "__main__":
    main()
