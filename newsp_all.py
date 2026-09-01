#!/usr/bin/env python3
"""Rebuild every new-species result from the raw alignments, in one pass.

The earlier run scored 17 of hydra's 22 candidates because the last five
alignments finished seven minutes after the scoring did. This redoes the whole
thing so nothing is stale:

  raw alignment -> flanks justified -> flank width trimmed -> scored
                                                           -> flank islands

Writes clean alignments, newsp_verdicts_all.json, flank_islands_all.json.
"""
import glob
import json
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, "/staging/tmp/sinedisc")
import justify_all as J
import trim_flanks as T

SRC = "/staging/tmp/newsp"
CLEAN = "/staging/tmp/newsp/clean"
GAP = "-"
MIN_RUN = 6          # an island is at least this many consecutive columns
Z_CUT = 8.0          # ... each above this many SDs of the coverage-aware null
MIN_COPIES = 8       # a column with fewer copies present tells us nothing


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


def flank_columns(path):
    """Left and right flank blocks as arrays, copies only, consensus dropped."""
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
    # raw alignments carry soft-masked lowercase straight from getfasta
    rows = [s.upper() for i, s in enumerate(seqs) if i != k]
    L = np.array([list(r[:lo]) for r in rows]) if lo > 0 else None
    R = np.array([list(r[hi + 1:]) for r in rows]) if hi + 1 < len(rows[0]) else None
    return L, R


def island_scan(path):
    """Columns where pairwise flank identity beats chance, controlling for how
    many copies are actually present in that column.

    The point of the control: a column where 9 copies of 100 happen to be
    present can show a high identity by chance alone. Scoring it against the
    same null as a column with 98 copies present is what made me call the
    far-upstream islands artifacts the first time. They are not."""
    fl = flank_columns(path)
    if fl is None:
        return None
    zs, tot_bases, tot_match = [], 0, 0
    blocks = [b for b in fl if b is not None and b.size]

    # background composition, pooled over both flanks of this set
    counts = {}
    for b in blocks:
        for base in "ACGT":
            counts[base] = counts.get(base, 0) + int((b == base).sum())
    tot = sum(counts.values())
    if tot < 500:
        return None
    p0 = sum((c / float(tot)) ** 2 for c in counts.values())

    for b in blocks:
        for j in range(b.shape[1]):
            col = b[:, j]
            col = col[(col != GAP) & np.isin(col, list("ACGT"))]
            n = len(col)
            if n < MIN_COPIES:
                zs.append(0.0)
                continue
            pairs = n * (n - 1) / 2.0
            _, cnt = np.unique(col, return_counts=True)
            match = float((cnt * (cnt - 1) / 2.0).sum())
            p = match / pairs
            sd = np.sqrt(p0 * (1 - p0) / pairs)
            zs.append((p - p0) / sd if sd > 0 else 0.0)
            tot_bases += n
            tot_match += match
    zs = np.asarray(zs)
    hot = zs > Z_CUT
    n_isl, n_col, run = 0, 0, 0
    for h in hot:
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
    return [n_isl, round(float(zs.max()) if zs.size else 0.0, 1), n_col]


def main():
    os.makedirs(CLEAN, exist_ok=True)
    made = []
    for sp in ("pom", "aca", "hyd"):
        for f in sorted(glob.glob(os.path.join(SRC, sp, "aln", "*.aln.fa"))):
            base = os.path.basename(f).replace(".aln.fa", "")
            if ".clean" in base or ".degap" in base:
                continue
            tmp = os.path.join(CLEAN, "_j_" + base + ".aln.fa")
            out = os.path.join(CLEAN, "NEW__" + base + ".clean.aln.fa")
            if not J.justify(f, tmp):
                print("  no consensus:", base)
                continue
            if not T.trim(tmp, out):
                os.replace(tmp, out)
            else:
                os.remove(tmp)
            made.append(out)
    print("clean alignments: %d" % len(made))

    r = subprocess.run([sys.executable, "/staging/tmp/sinedisc/verdict.py",
                        CLEAN, "/staging/tmp/sinedisc/newsp_verdicts_all.json"],
                       capture_output=True, text=True, cwd="/staging/tmp/sinedisc")
    print(r.stdout[-1500:], r.stderr[-800:])

    isl = {}
    for out in made:
        base = os.path.basename(out).replace("NEW__", "").replace(".clean.aln.fa", "")
        v = island_scan(out)
        if v:
            isl[base] = v
    json.dump(isl, open("/staging/tmp/sinedisc/flank_islands_all.json", "w"),
              indent=1, sort_keys=True)
    print("islands scanned: %d" % len(isl))


if __name__ == "__main__":
    main()
