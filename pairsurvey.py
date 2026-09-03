#!/usr/bin/env python3
"""Fraction-variable across every pair of curated subfamilies.

The question: does "what fraction of the element differs between these two
groups" separate the splits he keeps from the one he wants merged?

Calibration so far, from his own judgements:

    t3-1 / t3-2    9 % variable   he named the separation at a glance   ACCEPT
    t1-4 / t1-2   81 % variable   "more plausible to merge"             REJECT

His instruction on the bottom anchor: reject anything below the t1-4 case. So the
threshold sits between those two, and every other pair should be readable against
it.

All pairs are measured in ONE common alignment — the general subfam alignment —
because extracting a pair and realigning it alone destroys the comparison that
makes a separation visible. Rows are subset per pair; columns are not realigned.

Two numbers per pair:

  frac_var   fraction of the shared element span that is a variable site at 15 %,
             using the viewer's rule: minority count over sequences whose own
             span covers the column, gaps counted as a character
  best       the strongest single column, max over |gap-fraction difference| and
             |majority-base-fraction difference|, occupancy-guarded at 30 %

frac_var says how pervasive the difference is; best says whether any of it is
consistent. A real split is SPARSE and CONSISTENT: low frac_var, high best.
"""
import io
import json
import sys
from collections import Counter

import numpy as np

GAPS = "-."
PCT = 15
MIN_OCC = 0.30
MIN_CHUNKS = 8


def read_fa(p):
    n, s, cur, buf = [], [], None, []
    for line in io.open(p, encoding="utf-8", errors="replace"):
        line = line.rstrip()
        if line.startswith(">"):
            if cur is not None:
                s.append("".join(buf))
            cur = line[1:].split()[0]
            n.append(cur)
            buf = []
        else:
            buf.append(line.strip())
    if cur is not None:
        s.append("".join(buf))
    return n, s


def measure(rows_a, rows_b):
    rows = rows_a + rows_b
    N = len(rows)
    L = len(rows[0])
    spans = []
    for r in rows:
        f = l = -1
        for i, c in enumerate(r):
            if c not in GAPS:
                if f < 0:
                    f = i
                l = i
        spans.append((f, l))
    thr = int(np.ceil(PCT / 100.0 * N))

    span_cols = 0
    var_cols = 0
    best = 0.0
    na, nb = len(rows_a), len(rows_b)
    for pos in range(L):
        cnt = {}
        cov = 0
        for i, r in enumerate(rows):
            f, l = spans[i]
            if f < 0 or pos < f or pos > l:
                continue
            ch = r[pos].upper()
            ch = "-" if ch in GAPS else ch
            cnt[ch] = cnt.get(ch, 0) + 1
            cov += 1
        if cov < 0.5 * N:          # outside the shared element span
            continue
        span_cols += 1
        if cov - max(cnt.values()) >= thr:
            var_cols += 1
        # strongest single column
        ga = sum(1 for r in rows_a if r[pos] in GAPS) / float(na)
        gb = sum(1 for r in rows_b if r[pos] in GAPS) / float(nb)
        best = max(best, abs(ga - gb))
        oa, ob = 1 - ga, 1 - gb
        if oa >= MIN_OCC and ob >= MIN_OCC:
            ca = Counter(r[pos].upper() for r in rows_a if r[pos] not in GAPS)
            cb = Counter(r[pos].upper() for r in rows_b if r[pos] not in GAPS)
            t = ca.most_common(1)[0][0]
            fa = ca[t] / float(sum(ca.values()))
            fb = cb.get(t, 0) / float(sum(cb.values()))
            best = max(best, abs(fa - fb))
    frac = var_cols / float(span_cols) if span_cols else 0.0
    return frac, best, span_cols


def main():
    aln, truth_json = sys.argv[1], sys.argv[2]
    names, seqs = read_fa(aln)
    truth = json.load(open(truth_json))
    lab = {}
    for i, nm in enumerate(names):
        g = truth.get(nm) or truth.get(nm + ".bnk")
        if g:
            lab.setdefault(g, []).append(i)
    groups = {g: v for g, v in lab.items() if len(v) >= MIN_CHUNKS}
    print("groups with >= %d chunks: %s\n"
          % (MIN_CHUNKS, {g: len(v) for g, v in sorted(groups.items())}))
    ks = sorted(groups)
    rows = []
    for x in range(len(ks)):
        for y in range(x + 1, len(ks)):
            a, b = ks[x], ks[y]
            fr, be, sp = measure([seqs[i] for i in groups[a]],
                                 [seqs[i] for i in groups[b]])
            rows.append((fr, be, a, b, len(groups[a]), len(groups[b]), sp))
    rows.sort()
    print("%-10s %-10s %6s %6s %8s %8s   %s"
          % ("group A", "group B", "nA", "nB", "frac_var", "best", "reading"))
    for fr, be, a, b, na, nb, sp in rows:
        if fr <= 0.25 and be >= 0.80:
            call = "SPLIT   sparse and consistent"
        elif fr >= 0.60 and be < 0.70:
            call = "MERGE   pervasive and inconsistent"
        else:
            call = "grey"
        print("%-10s %-10s %6d %6d %8.2f %8.2f   %s"
              % (a, b, na, nb, fr, be, call))


if __name__ == "__main__":
    main()
