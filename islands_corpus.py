#!/usr/bin/env python3
"""Do flank similarity islands separate anything, or are they an aca_SINE_0 quirk?

One candidate is not a finding. This runs the same coverage-controlled island
scan over the whole labelled corpus, where the classes are known, and asks
whether island columns tell the classes apart at all - and specifically whether
they catch NEGSAT / NEGSEGDUP / NEGCHIM, the classes whose whole problem is that
the copies are not in independent places.

The measurement must be on the raw long flanks (aln_c), never the justified or
trimmed ones: trimming drops aca_SINE_0 from 426 island columns to 6.
"""
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np

GAP = "-"
MIN_RUN = 6
Z_CUT = 8.0
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


def flank_blocks(path):
    names, seqs = read_fa(path)
    ci = [i for i, n in enumerate(names) if "CONSENSUS_" in n.upper()]
    if not ci or len(seqs) < 12:
        return None
    k = ci[0]
    cons = seqs[k]
    nz = [i for i, c in enumerate(cons) if c != GAP]
    if len(nz) < 40:
        return None
    lo, hi = nz[0], nz[-1]
    rows = [s.upper() for i, s in enumerate(seqs) if i != k]
    out = []
    if lo > 0:
        out.append(np.array([list(r[:lo]) for r in rows]))
    if hi + 1 < len(rows[0]):
        out.append(np.array([list(r[hi + 1:]) for r in rows]))
    return out


def scan(path):
    blocks = flank_blocks(path)
    if not blocks:
        return None
    blocks = [b for b in blocks if b.size]
    if not blocks:
        return None

    counts = {}
    for b in blocks:
        for base in "ACGT":
            counts[base] = counts.get(base, 0) + int((b == base).sum())
    tot = sum(counts.values())
    if tot < 500:
        return None
    p0 = sum((c / float(tot)) ** 2 for c in counts.values())

    zs, width = [], 0
    for b in blocks:
        width += b.shape[1]
        for j in range(b.shape[1]):
            col = b[:, j]
            col = col[np.isin(col, list("ACGT"))]
            n = len(col)
            if n < MIN_COPIES:
                zs.append(0.0)
                continue
            pairs = n * (n - 1) / 2.0
            _, cnt = np.unique(col, return_counts=True)
            match = float((cnt * (cnt - 1) / 2.0).sum())
            sd = np.sqrt(p0 * (1 - p0) / pairs)
            zs.append((match / pairs - p0) / sd if sd > 0 else 0.0)

    zs = np.asarray(zs)
    n_isl, n_col, run = 0, 0, 0
    for h in zs > Z_CUT:
        if h:
            run += 1
        else:
            if run >= MIN_RUN:
                n_isl += 1
                n_col += run
            run = 0
    if run >= MIN_RUN:
        n_isl += 1
        n_col += run
    # the fraction matters as much as the count: a set with a 2000 bp flank has
    # more room for islands than one with 200
    return {"islands": n_isl, "max_z": round(float(zs.max()) if zs.size else 0.0, 1),
            "cols": n_col, "flank_width": width,
            "frac": round(n_col / float(width), 4) if width else 0.0}


def _one(f):
    try:
        r = scan(f)
    except Exception as exc:
        return os.path.basename(f), {"error": str(exc)}
    return os.path.basename(f).replace(".aln.fa", ""), r


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "aln_c"
    out = sys.argv[2] if len(sys.argv) > 2 else "islands_corpus.json"
    files = sorted(glob.glob(os.path.join(src, "*.aln.fa")))
    import multiprocessing as mp
    res = {}
    with mp.Pool(max(1, mp.cpu_count() - 4)) as pool:
        for name, r in pool.imap_unordered(_one, files, chunksize=4):
            if r and "error" not in r:
                res[name] = r
    json.dump(res, open(out, "w"), separators=(",", ":"))
    print("scanned %d of %d" % (len(res), len(files)))

    by = defaultdict(list)
    for k, v in res.items():
        by[k.split("__")[0]].append(v)
    print("\n%-12s %5s %8s %8s %8s %8s %8s"
          % ("class", "n", "med cols", "max cols", "med frac", "max frac", "med z"))
    print("-" * 62)
    order = sorted(by, key=lambda c: -np.median([v["cols"] for v in by[c]]))
    for c in order:
        v = by[c]
        print("%-12s %5d %8.0f %8d %8.3f %8.3f %8.1f"
              % (c, len(v),
                 np.median([x["cols"] for x in v]), max(x["cols"] for x in v),
                 np.median([x["frac"] for x in v]), max(x["frac"] for x in v),
                 np.median([x["max_z"] for x in v])))


if __name__ == "__main__":
    main()
