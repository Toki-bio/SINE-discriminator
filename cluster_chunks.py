#!/usr/bin/env python3
"""Group SubFam chunk consensuses into candidate families, and build a consensus
for each.

This is MANUAL section 6.1 done as code rather than by eye. SubFam collapses
thousands of loci into one consensus per ~50 of them; those chunk consensuses
then fall into groups, and each group is a candidate subfamily. The manual is
explicit that SubFam's own final pass is not a converged alignment, so the input
here must already have been degapped and realigned.

The clustering is deliberately plain: pairwise identity over aligned columns
both rows occupy, then single-link agglomeration at a fixed identity. Nothing
here is tuned - the automated assist that ships with SINEderella is documented
as recovering only the coarsest structure, so this is a first pass whose output
is meant to be looked at, not trusted.
"""
import io
import json
import os
import sys
from collections import Counter

import numpy as np

GAP = "-"


def read_fa(p):
    names, seqs, cur, buf = [], [], None, []
    for line in io.open(p, encoding="utf-8", errors="replace"):
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


def identity_matrix(A):
    n = A.shape[0]
    M = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            ok = (A[i] != GAP) & (A[j] != GAP)
            k = ok.sum()
            M[i, j] = M[j, i] = (A[i][ok] == A[j][ok]).mean() if k >= 40 else 0.0
    return M


def single_link(M, thr):
    n = M.shape[0]
    lab = list(range(n))

    def find(x):
        while lab[x] != x:
            lab[x] = lab[lab[x]]
            x = lab[x]
        return x

    for i in range(n):
        for j in range(i + 1, n):
            if M[i, j] >= thr:
                a, b = find(i), find(j)
                if a != b:
                    lab[a] = b
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    return sorted(groups.values(), key=len, reverse=True)


def consensus(A, rows, min_cov=0.5):
    sub = A[rows]
    out = []
    for j in range(sub.shape[1]):
        col = sub[:, j]
        col = col[col != GAP]
        if len(col) < max(2, min_cov * len(rows)):
            continue
        out.append(Counter(col).most_common(1)[0][0])
    return "".join(out)


def main():
    aln, out_dir, tag = sys.argv[1], sys.argv[2], sys.argv[3]
    thr = float(sys.argv[4]) if len(sys.argv) > 4 else 0.70
    min_members = int(sys.argv[5]) if len(sys.argv) > 5 else 3
    os.makedirs(out_dir, exist_ok=True)

    names, seqs = read_fa(aln)
    A = np.array([list(s.upper()) for s in seqs])
    print("%d chunk consensuses, %d columns" % A.shape)

    M = identity_matrix(A)
    off = M[np.triu_indices(len(names), 1)]
    print("pairwise identity: median %.3f, 90th %.3f, max %.3f"
          % (np.median(off), np.percentile(off, 90), off.max()))

    groups = single_link(M, thr)
    kept = [g for g in groups if len(g) >= min_members]
    print("at identity %.2f: %d groups, %d with %d+ members"
          % (thr, len(groups), len(kept), min_members))

    fa = os.path.join(out_dir, "%s_families.fa" % tag)
    meta = {}
    with io.open(fa, "w", encoding="utf-8") as fh:
        for i, g in enumerate(kept):
            seq = consensus(A, g)
            if len(seq) < 60:
                continue
            name = "%s_SINE_%d" % (tag, i)
            fh.write(u">%s\n%s\n" % (name, seq))
            within = [M[a, b] for a in g for b in g if a < b]
            meta[name] = {"members": len(g), "length": len(seq),
                          "within_identity": round(float(np.mean(within)), 3) if within else None}
            print("  %-14s %3d chunks, %4d bp, within-group identity %s"
                  % (name, len(g), len(seq),
                     "-" if not within else "%.3f" % np.mean(within)))
    json.dump(meta, open(os.path.join(out_dir, "%s_families.json" % tag), "w"), indent=1)
    print("wrote %s" % fa)


if __name__ == "__main__":
    main()
