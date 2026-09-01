#!/usr/bin/env python3
"""Do the copies agree with EACH OTHER more than with the consensus they were given?

This separates two failures the tool currently conflates under NO_ELEMENT:

  (a) there is no element - the copies do not agree with each other either
  (b) there IS an element, but the supplied consensus does not describe it

(b) is recoverable and is the normal state of an upstream de novo candidate.
(a) is a rejection. Telling them apart needs two numbers that are already
computed but never compared:

  pair_id   plateau of copy-to-copy identity across the element (from profiles)
  cons_id   identity of copies to the supplied consensus

A real family with a bad consensus shows pair_id well above cons_id. A genuine
non-family shows both at background.
"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, ".")
import profiles as P
import verdict as V


def smooth(y, w=9):
    y = np.asarray(y, float)
    ok = np.isfinite(y)
    if ok.sum() < 3:
        return y
    filled = np.interp(np.arange(len(y)), np.flatnonzero(ok), y[ok])
    return np.convolve(filled, np.ones(w) / w, mode="same")


def check(path):
    pr = P.profile(path)
    if pr is None:
        return None
    xs = np.asarray(pr["x"], float)
    pid = smooth(pr["pair_id"])
    cid = smooth(pr["cons_id"])
    inside = xs >= 0
    if inside.sum() < 60:
        return None
    L = int(xs[inside].max()) + 1
    mid = inside & (xs > L * 0.25) & (xs < L * 0.75)
    if mid.sum() < 20:
        return None

    pair = float(np.nanpercentile(pid[mid], 75))
    cons = float(np.nanmedian(cid[mid]))
    far = (xs <= -35) | (xs >= L + 34)
    bg = float(np.nanmedian(pid[far])) if far.sum() else 0.25

    try:
        v = V.verdict(path)
        score = v["score"] if v else None
        flags = ",".join(f["code"] for f in v["flags"]) if v else ""
    except Exception:
        score, flags = None, ""

    return {"set": os.path.basename(path).replace(".aln.fa", ""),
            "L": L, "pair": pair, "cons": cons, "bg": bg,
            "pair_minus_bg": pair - bg,
            "pair_minus_cons": pair - cons,
            "score": score, "flags": flags}


def main():
    args = sys.argv[1:]
    files = []
    for a in args:
        files.extend(sorted(glob.glob(a)) if any(c in a for c in "*?") else [a])
    rows = []
    for f in files:
        try:
            r = check(f)
        except Exception:
            r = None
        if r:
            rows.append(r)
    rows.sort(key=lambda r: -r["pair_minus_cons"])
    print("%-30s %6s %6s %6s %8s %8s %6s  %s"
          % ("set", "pair", "cons", "bg", "pair-bg", "pair-cons", "score", "flags"))
    for r in rows:
        print("%-30s %6.3f %6.3f %6.3f %8.3f %8.3f %6s  %s"
              % (r["set"][:30], r["pair"], r["cons"], r["bg"],
                 r["pair_minus_bg"], r["pair_minus_cons"],
                 "-" if r["score"] is None else "%.1f" % r["score"], r["flags"][:44]))


if __name__ == "__main__":
    main()
