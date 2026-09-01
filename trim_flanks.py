#!/usr/bin/env python3
"""Trim flank display width so a few long copies stop swamping the view.

Sergei, on the review batch: "flanks are badly degapped and it hinders my
estimates".

The flanks ARE justified - every copy's flank sits hard against the element,
gap_to_element is 0 for all of them. The problem is the column WIDTH: it is set
by the single longest copy, while the median copy has far fewer bases. On
NEGCHIM__ccr__g3_71seqs the right flank is 125 columns wide, the median copy has
10 bases in it, and the panel is 85% gaps. One outlier pads seventy other copies
with whitespace.

Fix: cap each flank at a percentile of the copies' own flank lengths, so the
panel is sized for the copies that are actually there. Copies with more flank
than that are cut - their extra sequence is unalignable anyway and carries no
information for judging the boundary.

The element itself is never touched.
"""
import io
import os
import sys

import numpy as np

sys.path.insert(0, ".")
import measure_c as M

SYM = "ACGT-"
PCTL = 75        # keep this percentile of flank length
MINCOL = 25      # never trim below this, even if most copies have no flank
MAXGAP = 0.25    # the flank panel may be at most this fraction gaps


def read_named(path):
    names, seqs, cur, buf = [], [], None, []
    for line in io.open(path, encoding="utf-8", errors="replace"):
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


def trim(path, out, pctl=PCTL):
    names, seqs = read_named(path)
    if not seqs:
        return False
    ci = [i for i, n in enumerate(names) if "CONSENSUS_" in n.upper()]
    if not ci:
        return False
    k = ci[0]
    cons = seqs[k]
    nz = [i for i, c in enumerate(cons) if c != "-"]
    if len(nz) < 40:
        return False
    lo, hi = nz[0], nz[-1]

    others = [s for i, s in enumerate(seqs) if i != k]
    # how much flank does each copy actually carry?
    lb = [sum(1 for c in s[:lo] if c != "-") for s in others]
    rb = [sum(1 for c in s[hi + 1:] if c != "-") for s in others]
    # A percentile is not enough when flank lengths are very skewed. On
    # NEGCHIM__ccr__g1_180seqs the right flank has median 20 bases but 75th
    # percentile 58, so a p75 cap still left the panel 53% gaps and Sergei
    # rightly said the right flank was still not degapped.
    #
    # Choose instead the widest cap whose panel is at most MAXGAP gaps: the
    # display is then sized for the copies that are actually there.
    def width(lens, avail):
        if not lens:
            return min(MINCOL, avail)
        lens = np.asarray(lens)
        # The floor must follow the copies, not a constant. With MINCOL=25 and a
        # median right flank of 10 bases the panel is 60% gaps no matter what -
        # which is what Sergei was still seeing on NEGCHIM__ccr__g3_71seqs.
        floor = int(max(5, min(MINCOL, np.median(lens))))
        best = min(floor, avail)
        for w in range(min(floor, avail), avail + 1, 5):
            filled = np.minimum(lens, w).sum()
            if filled / float(w * len(lens)) >= 1.0 - MAXGAP:
                best = w
        return best

    keepL = width(lb, lo)
    keepR = width(rb, len(cons) - hi - 1)

    # flanks are already justified against the element, so the informative
    # columns are the ones nearest the element: keep the inner slice
    with io.open(out, "w", encoding="utf-8") as fh:
        for n, s in zip(names, seqs):
            left = s[lo - keepL:lo] if keepL else ""
            right = s[hi + 1:hi + 1 + keepR] if keepR else ""
            fh.write(">%s\n%s%s%s\n" % (n, left, s[lo:hi + 1], right))
    return True


def main():
    src = sys.argv[1]
    dst = sys.argv[2]
    names = sys.argv[3:]
    os.makedirs(dst, exist_ok=True)
    import glob
    files = ([os.path.join(src, n + ".aln.fa") for n in names] if names
             else sorted(glob.glob(os.path.join(src, "*.aln.fa"))))
    ok = 0
    for p in files:
        if not os.path.exists(p):
            continue
        o = os.path.join(dst, os.path.basename(p))
        try:
            if trim(p, o):
                ok += 1
        except Exception as exc:
            print("  fail %s %s" % (os.path.basename(p), exc))
    print("trimmed %d of %d" % (ok, len(files)))


if __name__ == "__main__":
    main()
