#!/usr/bin/env python3
"""Is the true element edge detectable, independently of the supplied consensus?

The question this answers is NOT "is the consensus wrong" - the tool already
flags that. It is whether the copy-to-copy identity profile carries a clean
enough transition to RELOCATE the boundary, which is what an edge-adjustment
loop would need.

Method: take profile()'s pairwise-identity curve, which is defined identically
in all three zones (left flank / element / right flank) with flanks indexed
outward from each copy's own ungapped edge. Estimate the element plateau and
the genomic background, then find where the curve crosses the midpoint between
them, scanning outward from the element centre on each side independently.

That crossing is the data's own opinion of where the element stops. Compare it
to the consensus edge, which is at 0 on the left axis and L-1 on the right.
"""
import json
import sys

import numpy as np

sys.path.insert(0, ".")
import profiles as P


def smooth(y, w=9):
    y = np.asarray(y, float)
    ok = np.isfinite(y)
    if ok.sum() < 3:
        return y
    filled = np.interp(np.arange(len(y)), np.flatnonzero(ok), y[ok])
    k = np.ones(w) / w
    return np.convolve(filled, k, mode="same")


def edges(path, plateau_q=0.75):
    pr = P.profile(path)
    if pr is None:
        return None
    xs = np.asarray(pr["x"], float)
    pid = smooth(pr["pair_id"])

    inside = (xs >= 0)
    if inside.sum() < 60:
        return None
    L = int(xs[inside].max()) + 1

    # element plateau: upper quartile of identity over the middle half of the
    # element, so a bad end cannot drag it down
    mid = inside & (xs > L * 0.25) & (xs < L * 0.75)
    plateau = float(np.nanpercentile(pid[mid], 100 * plateau_q))

    # background: the far half of each flank, away from any edge effect
    farL = xs <= -35
    farR = xs >= L + 35 - 1
    bg_vals = np.concatenate([pid[farL], pid[farR]])
    bg = float(np.nanmedian(bg_vals)) if bg_vals.size else 0.25

    if plateau - bg < 0.12:
        return {"path": path, "L": L, "plateau": round(plateau, 3),
                "bg": round(bg, 3), "usable": False}

    half = (plateau + bg) / 2.0

    # scan outward from the element centre; the crossing is the last position
    # still above the half-way level
    centre = L / 2.0
    left_x, right_x = None, None
    order = np.argsort(xs)
    xo, po = xs[order], pid[order]

    for i in range(len(xo)):
        if xo[i] >= centre:
            break
        if po[i] >= half and left_x is None:
            left_x = xo[i]
    # take the LAST crossing before the centre, i.e. the innermost sustained rise
    cand = [xo[i] for i in range(len(xo)) if xo[i] < centre and po[i] >= half]
    left_x = min(cand) if cand else None
    cand = [xo[i] for i in range(len(xo)) if xo[i] >= centre and po[i] >= half]
    right_x = max(cand) if cand else None

    return {"path": path, "L": L,
            "plateau": round(plateau, 3), "bg": round(bg, 3),
            "half": round(half, 3),
            "left_edge": None if left_x is None else float(left_x),
            "right_edge": None if right_x is None else float(right_x),
            "left_shift": None if left_x is None else float(left_x - 0),
            "right_shift": None if right_x is None else float(right_x - (L - 1)),
            "usable": True}


def main():
    for p in sys.argv[1:]:
        r = edges(p)
        if r is None:
            print("%-46s  no profile" % p)
            continue
        if not r["usable"]:
            print("%-46s  L=%-4d plateau=%.3f bg=%.3f  NO CONTRAST"
                  % (p, r["L"], r["plateau"], r["bg"]))
            continue
        print("%-46s  L=%-4d plateau=%.3f bg=%.3f half=%.3f  left_shift=%+.0f  right_shift=%+.0f"
              % (p, r["L"], r["plateau"], r["bg"], r["half"],
                 r["left_shift"], r["right_shift"]))


if __name__ == "__main__":
    main()
