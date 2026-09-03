#!/usr/bin/env python3
"""MAFFT's 6-mer distance, ported exactly from mafft/core.

From disttbfast.c and mltaln9.c (GSLBiotech/mafft master):

    tuplesize = 6                                    disttbfast.c:214
    shared    = commonsextet_p(A, B)                 mltaln9.c:14871
    bunbo     = MIN(selfscore_i, selfscore_j)        disttbfast.c:4053
    lenfac    = 1/( shorter/longer*0.1 + 2500/(longer+2500) + 0.01 )
    dist      = (1 - shared/bunbo) * lenfac * 2.0    disttbfast.c:4057

with D6LENFACA 0.01, D6LENFACB 2500, D6LENFACC 2500, D6LENFACD 0.1.

`commonsextet_p` is a MULTISET intersection: walking B's k-mer occurrences, an
occurrence counts only while B has not yet used up A's copies of that k-mer.

Why this matters for subfamily calling, and it is the reason he pointed at
MAFFT's ordering:

1. It normalises by **MIN(self)**, not by a union. A short sequence wholly
   contained in a longer one scores 1.0. That is containment, and it is tolerant
   of truncation - which is what SINE copies of varying completeness need.
2. It **aggregates weak evidence** across every k-mer. No single column has to be
   diagnostic, so a difference that is faint at every position but consistent
   still accumulates. That is "faint and real at the same time".
3. **An indel is weighted by its k-mer footprint, not its column count.** A
   deletion of length L destroys L+k-1 overlapping k-mers; a substitution
   destroys at most k. So a 5 bp indel counts about 10 against a SNP's 6 -
   indels outweigh substitutions automatically, which is the weighting MANUAL
   6.1.6 asks for ("small indels/SNPs", indels first). My peel counts a 5 bp
   indel as 5 features and a SNP as 1, under-weighting indels the other way.

His edge rule falls out of the same arithmetic: k-mers spanning a ragged edge are
mostly unique, so they inflate `bunbo` without ever contributing to `shared`.
Ragged edges therefore inflate distance systematically, which is why they must be
trimmed before the distance is computed rather than after.
"""
import sys
from collections import Counter

K = 6
LENFACA, LENFACB, LENFACC, LENFACD = 0.01, 2500.0, 2500.0, 0.1
CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def kmers(seq):
    """k-mer multiset. Non-ACGT is dropped, as seq_grp_nuc drops tmp >= 4."""
    s = [CODE[c] for c in seq.upper() if c in CODE]
    if len(s) < K:
        return Counter(), 0
    c = Counter()
    v = 0
    for i in range(K - 1):
        v = v * 4 + s[i]
    m = 4 ** (K - 1)          # place value of the OLDEST digit in a K-digit number
    for i in range(K - 1, len(s)):
        v = v * 4 + s[i]
        c[v] += 1
        v -= s[i - K + 1] * m  # drop the oldest digit; *m, not *m*4
    return c, len(s)


def selfscore(c):
    return sum(c.values())


def shared(ca, cb):
    """commonsextet_p: multiset intersection size."""
    return sum(min(v, ca.get(k, 0)) for k, v in cb.items())


def lenfac(li, lj):
    longer, shorter = (li, lj) if li >= lj else (lj, li)
    if longer == 0:
        return 1.0
    return 1.0 / (shorter / float(longer) * LENFACD
                  + LENFACB / (longer + LENFACC) + LENFACA)


def distance(sa, sb, ca=None, cb=None, la=None, lb=None):
    if ca is None:
        ca, la = kmers(sa)
    if cb is None:
        cb, lb = kmers(sb)
    bunbo = min(selfscore(ca), selfscore(cb))
    if bunbo == 0:
        return 2.0
    return (1.0 - shared(ca, cb) / float(bunbo)) * lenfac(la, lb) * 2.0


def main():
    import io
    path = sys.argv[1]
    names, seqs, cur, buf = [], [], None, []
    for line in io.open(path, encoding="utf-8", errors="replace"):
        line = line.rstrip()
        if line.startswith(">"):
            if cur is not None:
                seqs.append("".join(buf))
            cur = line[1:]
            names.append(cur)
            buf = []
        else:
            buf.append(line.strip().replace("-", ""))
    seqs.append("".join(buf))

    groups = {}
    for i, h in enumerate(names):
        g = h.split("__")[0]
        groups.setdefault(g, []).append(i)
    pre = [kmers(s) for s in seqs]

    keys = [g for g in sorted(groups) if not g.startswith("CONSENSUS")]
    print("MAFFT 6-mer distance, group means\n")
    print("%-10s %5s   %s" % ("group", "n", "  ".join("%9s" % k for k in keys)))
    for a in keys:
        row = []
        for b in keys:
            ds = []
            for i in groups[a]:
                for j in groups[b]:
                    if i < j or a != b:
                        ds.append(distance(None, None, pre[i][0], pre[j][0],
                                           pre[i][1], pre[j][1]))
            row.append(sum(ds) / len(ds) if ds else 0.0)
        print("%-10s %5d   %s" % (a, len(groups[a]),
                                  "  ".join("%9.4f" % v for v in row)))
    print("\nseparation = between-group mean minus the larger within-group mean")
    for x in range(len(keys)):
        for y in range(x + 1, len(keys)):
            a, b = keys[x], keys[y]
            def mean(g1, g2):
                ds = [distance(None, None, pre[i][0], pre[j][0], pre[i][1], pre[j][1])
                      for i in groups[g1] for j in groups[g2] if i < j or g1 != g2]
                return sum(ds) / len(ds) if ds else 0.0
            wa, wb, bt = mean(a, a), mean(b, b), mean(a, b)
            print("   %-8s vs %-8s   within %.4f / %.4f   between %.4f   sep %+.4f"
                  % (a, b, wa, wb, bt, bt - max(wa, wb)))


if __name__ == "__main__":
    main()
