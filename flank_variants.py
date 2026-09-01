#!/usr/bin/env python3
"""Produce several flank presentations of the same alignment, to be compared.

He keeps saying the flanks are wrong and my guesses about why have been wrong
twice. So rather than guess a third time: build the plausible presentations of
the same data and let him say which one is right.

A FASTA alignment is rectangular, so a copy with a short flank must be padded
somewhere. The variants differ in HOW MUCH flank is shown and WHERE the padding
goes. The element alignment is byte-identical in all of them.

  A  flush     flanks de-gapped and butted against the element, padded on the
               OUTER side to the longest copy. This is the current aln_v2.
  B  trimmed   same, but the width is capped so the panel is at most 25% gaps -
               a few long copies no longer pad everyone else.
  C  short     every flank cut to at most 30 bp. Nearly no padding; you see only
               the sequence immediately outside the element, which is where the
               boundary shows.
  D  ragged    no padding at all - each copy's flank written at its own length,
               so rows have different lengths. NOT a valid alignment file; some
               viewers will refuse it. Included because it is the only true
               "unaligned" presentation.
"""
import io
import os
import sys

GAP = "-"


def read_fa(p):
    names, seqs, cur, buf = [], [], None, []
    for l in io.open(p, encoding="utf-8", errors="replace"):
        l = l.rstrip("\n\r")
        if l.startswith(">"):
            if cur is not None:
                seqs.append("".join(buf))
            cur = l[1:]
            names.append(cur)
            buf = []
        else:
            buf.append(l.strip())
    if cur is not None:
        seqs.append("".join(buf))
    return names, seqs


def split(names, seqs):
    ci = [i for i, n in enumerate(names) if "CONSENSUS_" in n.upper()]
    if not ci:
        return None
    k = ci[0]
    cons = seqs[k]
    nz = [i for i, c in enumerate(cons) if c != GAP]
    if len(nz) < 40:
        return None
    lo, hi = nz[0], nz[-1]
    L = [s[:lo].replace(GAP, "") for s in seqs]
    E = [s[lo:hi + 1] for s in seqs]          # element alignment kept as-is
    R = [s[hi + 1:].replace(GAP, "") for s in seqs]
    return L, E, R


def write(path, names, L, E, R, capL, capR, pad=True):
    with io.open(path, "w", encoding="utf-8") as fh:
        for n, l, e, r in zip(names, L, E, R):
            l = l[-capL:] if capL else ""
            r = r[:capR] if capR else ""
            if pad:
                l = l.rjust(capL, GAP)
                r = r.ljust(capR, GAP)
            fh.write(">%s\n%s%s%s\n" % (n, l, e, r))


def width_for_gap(lens, avail, maxgap=0.25, floor=5):
    if not lens or avail <= 0:
        return 0
    import numpy as np
    lens = np.asarray(lens)
    start = int(max(floor, min(25, np.median(lens))))
    best = min(start, avail)
    for w in range(min(start, avail), avail + 1, 5):
        if np.minimum(lens, w).sum() / float(w * len(lens)) >= 1.0 - maxgap:
            best = w
    return best


def main():
    src, dst = sys.argv[1], sys.argv[2]
    os.makedirs(dst, exist_ok=True)
    made = 0
    for name in sys.argv[3:]:
        p = os.path.join(src, name + ".aln.fa")
        if not os.path.exists(p):
            print("  missing", name)
            continue
        names, seqs = read_fa(p)
        sp = split(names, seqs)
        if sp is None:
            print("  no consensus", name)
            continue
        L, E, R = sp
        lb = [len(x) for i, x in enumerate(L) if "CONSENSUS_" not in names[i].upper()]
        rb = [len(x) for i, x in enumerate(R) if "CONSENSUS_" not in names[i].upper()]
        maxL, maxR = max(lb or [0]), max(rb or [0])

        write(os.path.join(dst, name + ".A_flush.aln.fa"), names, L, E, R, maxL, maxR)
        write(os.path.join(dst, name + ".B_trimmed.aln.fa"), names, L, E, R,
              width_for_gap(lb, maxL), width_for_gap(rb, maxR))
        write(os.path.join(dst, name + ".C_short.aln.fa"), names, L, E, R,
              min(30, maxL), min(30, maxR))
        write(os.path.join(dst, name + ".D_ragged.aln.fa"), names, L, E, R,
              maxL, maxR, pad=False)
        made += 1

        import numpy as np
        def gapfrac(cap, lens):
            if not cap or not lens:
                return 0.0
            return 1.0 - np.minimum(np.asarray(lens), cap).sum() / float(cap * len(lens))
        print("%-38s L/R shown: A %d/%d (gaps %.2f/%.2f)  B %d/%d (%.2f/%.2f)  C %d/%d (%.2f/%.2f)"
              % (name[:38], maxL, maxR, gapfrac(maxL, lb), gapfrac(maxR, rb),
                 width_for_gap(lb, maxL), width_for_gap(rb, maxR),
                 gapfrac(width_for_gap(lb, maxL), lb), gapfrac(width_for_gap(rb, maxR), rb),
                 min(30, maxL), min(30, maxR),
                 gapfrac(min(30, maxL), lb), gapfrac(min(30, maxR), rb)))
    print("built %d sets x 4 variants" % made)


if __name__ == "__main__":
    main()
