#!/usr/bin/env python3
"""Do the same copies group together in every part of the element?

This is the distinction Sergei keeps drawing and the tool has never measured.
Two sets can both split their copies into two groups:

  subfamilies   the SAME copies are in the same group everywhere along the
                element. Two real families, and he is explicit that this does
                not matter: "subfamily ambiguity doesn't matter, we answer
                sine/not-sine".
  a mosaic      group membership SHUFFLES between regions - a copy sides with
                one set of copies at the 5' end and a different set in the
                middle. That is not two families, it is patchwork.

His second reading of hyd_SINE_5, with a screenshot of the left end: "looks like
mosaic, i still cant decide but its not a good sine if at all." The screenshot
shows a short conserved block at the 5' end that only some copies carry.

Method: cut the element into thirds, and in each third split the copies at the
median of their identity to the consensus. Then compare the partitions pairwise
with the adjusted Rand index. Near 1 means the same split everywhere -
subfamilies. Near 0 means membership is independent between regions - mosaic.

Four earlier attempts at this property all failed (edge sharpness, regional
patch2d, whole-element patch2d, per-copy boundary walking). They all asked how
MUCH copies differ. This asks whether they differ CONSISTENTLY, which is the
thing he is actually looking at.
"""
import glob
import itertools
import json
import os
import sys
from collections import defaultdict

import numpy as np

GAP = "-"
MIN_COPIES = 20
MIN_REGION = 40


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


def adjusted_rand(a, b):
    """Agreement between two 2-way partitions of the same copies."""
    a, b = np.asarray(a), np.asarray(b)
    n = len(a)
    tab = np.zeros((2, 2))
    for i in (0, 1):
        for j in (0, 1):
            tab[i, j] = np.sum((a == i) & (b == j))
    comb = lambda x: x * (x - 1) / 2.0
    sum_ij = comb(tab).sum()
    sum_i = comb(tab.sum(axis=1)).sum()
    sum_j = comb(tab.sum(axis=0)).sum()
    total = comb(float(n))
    exp = sum_i * sum_j / total if total else 0.0
    mx = (sum_i + sum_j) / 2.0
    return (sum_ij - exp) / (mx - exp) if mx != exp else 1.0


def measure(path, n_regions=3):
    names, seqs = read_fa(path)
    ci = [i for i, n in enumerate(names) if "CONSENSUS_" in n.upper()]
    if not ci or len(seqs) < MIN_COPIES + 1:
        return None
    k = ci[0]
    cons = seqs[k].upper()
    nz = [i for i, c in enumerate(cons) if c != GAP]
    if len(nz) < n_regions * MIN_REGION:
        return None
    lo, hi = nz[0], nz[-1]
    rows = [s.upper() for i, s in enumerate(seqs) if i != k]
    A = np.array([list(r[lo:hi + 1]) for r in rows])
    C = np.array(list(cons[lo:hi + 1]))
    L = A.shape[1]
    edges = [round(L * i / float(n_regions)) for i in range(n_regions + 1)]

    parts, ids = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        if b - a < MIN_REGION:
            return None
        sub, csub = A[:, a:b], C[a:b]
        ok = (sub != GAP) & (csub[None, :] != GAP)
        hit = (sub == csub[None, :]) & ok
        with np.errstate(invalid="ignore", divide="ignore"):
            frac = np.where(ok.sum(1) >= MIN_REGION * 0.4,
                            hit.sum(1) / np.maximum(ok.sum(1), 1), np.nan)
        if np.isfinite(frac).sum() < MIN_COPIES:
            return None
        ids.append(frac)
        med = np.nanmedian(frac)
        parts.append(np.where(np.nan_to_num(frac, nan=med) >= med, 1, 0))

    # only compare copies measurable in every region
    good = np.all([np.isfinite(f) for f in ids], axis=0)
    if good.sum() < MIN_COPIES:
        return None
    parts = [p[good] for p in parts]
    ars = [adjusted_rand(x, y) for x, y in itertools.combinations(parts, 2)]

    return {"consistency": round(float(np.mean(ars)), 3),
            "min_pair": round(float(min(ars)), 3),
            "region_id": [round(float(np.nanmedian(f[good])), 3) for f in ids],
            "n": int(good.sum())}


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
        by[k.split("__")[0]].append(v["consistency"])
    print("\n%-16s %5s %9s %9s %9s" % ("class", "n", "median", "min", "max"))
    print("-" * 54)
    for c in sorted(by, key=lambda c: np.median(by[c])):
        v = by[c]
        print("%-16s %5d %9.3f %9.3f %9.3f"
              % (c[:16], len(v), np.median(v), min(v), max(v)))


if __name__ == "__main__":
    main()
