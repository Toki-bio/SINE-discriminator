#!/usr/bin/env python3
"""Test 3 - do the new measurements earn their place?

Gate from PLAN.md: a new measurement earns its place only if it captures a
property nothing captured before, without spoiling one that already worked.

So: join Test 1 and Test 2 measurements, and for every pair of things Sergei
described, ask which single measurement separates that pair best and by how
much. Separation is reported as the gap between the two groups' means divided
by their pooled spread - so a big number means the two do not overlap.
"""
import io
import sys
from collections import defaultdict
from itertools import combinations

import numpy as np


def load(path):
    lines = io.open(path, encoding="utf-8").read().rstrip("\n").split("\n")
    hdr = lines[0].split("\t")
    out = {}
    for l in lines[1:]:
        f = l.split("\t")
        d = {}
        for k, v in zip(hdr[2:], f[2:]):
            try:
                d[k] = float(v)
            except Exception:
                pass
        out[f[0]] = (f[1], d)
    return out


MIN_N = 3   # a group of one has zero spread, which makes any gap look infinite


def sep(a, b):
    """Gap between two groups relative to their spread.

    Requires MIN_N on both sides. A first version did not, and every pair
    involving a single-example group returned ~1e11 - the group's variance is
    zero, so the pooled spread collapses to the epsilon. It reported "0 of 120
    pairs not told apart", which was meaningless.
    """
    a = np.asarray([x for x in a if x == x], float)
    b = np.asarray([x for x in b if x == x], float)
    if len(a) < MIN_N or len(b) < MIN_N:
        return np.nan
    pooled = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    if pooled < 1e-6:
        return np.nan
    return abs(a.mean() - b.mean()) / pooled


def main():
    t1 = load("test1_vars.tsv")
    t2 = load("test2_vars.tsv")
    merged = {}
    for st in set(t1) | set(t2):
        call = (t1.get(st) or t2.get(st))[0]
        d = {}
        d.update((t1.get(st) or (None, {}))[1])
        d.update((t2.get(st) or (None, {}))[1])
        merged[st] = (call, d)

    NEW = {"rowsplit", "rowsplit_frac", "rowsplit_dlen", "rowsplit_dident",
           "rowsplit_small", "redge_shift", "redge_naive_shift",
           "redge_tail_len", "redge_pos", "redge_usable",
           "dip_depth", "dip_width", "dip_pos"}

    g = defaultdict(lambda: defaultdict(list))
    for st, (c, d) in merged.items():
        for k, v in d.items():
            g[c][k].append(v)

    cats = [c for c in g if len(set(g[c].keys())) and
            sum(len(v) for v in g[c].values())]
    keys = sorted({k for c in g for k in g[c]})

    print("alignments %d   measurements %d   things you described %d"
          % (len(merged), len(keys), len(cats)))
    print()
    print("For each pair, the single measurement that separates it best.")
    print("'gap' is how far apart the two groups are, in units of their own spread.")
    print("A gap under 1.0 means they overlap and are not really told apart.")
    print()
    big = [c for c in cats if max((len(v) for v in g[c].values()), default=0) >= MIN_N]
    print("Only groups with at least %d judged alignments can be compared." % MIN_N)
    print("That is %d of your %d: %s" % (len(big), len(cats), ", ".join(sorted(big))))
    print("The other %d rest on one or two alignments each - Test 4 exists to fix that."
          % (len(cats) - len(big)))
    print()
    print("%-24s %-24s %-20s %6s %5s" % ("", "", "best measurement", "gap", "new?"))
    rows = []
    for a, b in combinations(sorted(big), 2):
        best = (None, -1)
        for k in keys:
            va, vb = g[a].get(k, []), g[b].get(k, [])
            if len(va) < 1 or len(vb) < 1:
                continue
            s = sep(va, vb)
            if s == s and s > best[1]:
                best = (k, s)
        if best[0]:
            rows.append((a, b, best[0], best[1], best[0] in NEW))
    rows.sort(key=lambda r: -r[3])
    for a, b, k, s, isnew in rows:
        print("%-24s %-24s %-20s %6.2f %5s" % (a[:24], b[:24], k, s, "NEW" if isnew else ""))

    weak = [r for r in rows if r[3] < 1.0]
    print()
    print("pairs still not told apart (gap < 1.0): %d of %d" % (len(weak), len(rows)))
    for a, b, k, s, _ in weak:
        print("   %-24s vs %-24s  best %s = %.2f" % (a[:24], b[:24], k, s))

    newbest = [r for r in rows if r[4]]
    print()
    print("pairs where a NEW measurement is the best one: %d of %d" % (len(newbest), len(rows)))
    from collections import Counter
    print("which new measurements are doing the work:",
          dict(Counter(r[2] for r in newbest).most_common()))


if __name__ == "__main__":
    main()
