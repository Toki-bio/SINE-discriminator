#!/usr/bin/env python3
"""Derive the element window from the copies, not from the consensus.

The defect this replaces, in measure_c.py:91 and justify_all.py:24:

    nz = np.where(cons != GAP)[0]
    lo, hi = int(nz[0]), int(nz[-1])

The window is the consensus's first and last non-gap column, taken on faith.
Nothing checks it against the copies, so an over-extended consensus puts real
flank inside the element (where justify never degaps it) and an under-extended
one puts real element outside (where justify degaps it away).

The rule here is step7's, applied per column of an existing alignment instead of
by re-extracting windows from the genome:

  occupancy(col)    fraction of copies with a base there
  conservation(col) fraction of those bases that are the majority base

step7 calls a window background when its ELEVATED fraction - pairs above 45 %
identity - falls to what random genomic regions give. The per-column analogue of
"elevated" is conservation > 0.45; random DNA sits near 0.25-0.30, and measured
on soft-masked Timema the true genomic background is 0.272.

A column is part of the element when occupancy >= OCC and conservation >= CONS.
Starting from the consensus's own supported core, walk outward on each side and
stop after MISS consecutive failing columns, so a single ragged column does not
end the element. That both trims an over-extended consensus and extends an
under-extended one, with the same rule and no per-set tuning.
"""
import numpy as np

OCC = 0.50      # half the copies must reach the column
CONS = 0.45     # step7's elevated cutoff
MISS = 8        # consecutive failures that end the element
GAPS = set("-.")


def column_stats(rows, j):
    b = [r[j] for r in rows if j < len(r) and r[j] not in GAPS]
    if not b:
        return 0.0, 0.0
    occ = len(b) / float(len(rows))
    u = [x.upper() for x in b]
    top = max(u.count(c) for c in set(u))
    return occ, top / float(len(u))


def element_window(cons, rows):
    """Return (lo, hi, diagnostics) for the element, measured on the copies."""
    L = len(cons)
    nz = [i for i, c in enumerate(cons) if c not in GAPS]
    if not nz or not rows:
        return (0, L - 1, {})
    old_lo, old_hi = nz[0], nz[-1]

    ok = {}

    def good(j):
        if j not in ok:
            o, c = column_stats(rows, j)
            ok[j] = (o >= OCC and c >= CONS)
        return ok[j]

    # start from the consensus columns that are themselves supported; if none
    # are, fall back to the consensus's own span so we never return nothing
    core = [j for j in nz if good(j)]
    if not core:
        return (old_lo, old_hi, {"note": "no supported consensus column"})
    lo, hi = core[0], core[-1]

    miss = 0
    j = lo - 1
    while j >= 0 and miss < MISS:
        if good(j):
            lo = j
            miss = 0
        else:
            miss += 1
        j -= 1

    miss = 0
    j = hi + 1
    while j < L and miss < MISS:
        if good(j):
            hi = j
            miss = 0
        else:
            miss += 1
        j += 1

    d = {
        "old_lo": old_lo, "old_hi": old_hi, "new_lo": lo, "new_hi": hi,
        "trimmed_left": max(0, lo - old_lo), "trimmed_right": max(0, old_hi - hi),
        "extended_left": max(0, old_lo - lo), "extended_right": max(0, hi - old_hi),
        "old_span": old_hi - old_lo + 1, "new_span": hi - lo + 1,
        "cons_bp": len(nz),
        "hit_left_edge": lo == 0, "hit_right_edge": hi == L - 1,
    }
    return (lo, hi, d)
