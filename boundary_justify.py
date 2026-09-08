#!/usr/bin/env python3
"""Copy-supported element window + display justify for publish alignments.

Imports boundary.element_window and fix_alignments.justify — does not
reimplement either. Replaces extract_alignments.sh postprocess_flanks on
MAFFT output (must run on the raw alignment, before any consensus-only
degapping).
"""
import glob
import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import boundary as B
from fix_alignments import consensus_index, justify, read_fa


def flanked_element_window(cons, rows):
    """Copy-supported window for 50L/70R publish alignments.

    boundary.element_window walks outward on copy support. Border loop + rebuilt
    consensus must fix the anchor row before publish; do not cap 3' extension
    back to the seed consensus span (that blocked SINE10 GAC recovery).
    """
    lo_b, hi_b, diag = B.element_window(cons, rows)
    lo = lo_b
    hi = hi_b
    diag["display_lo"] = lo
    diag["display_hi"] = hi
    diag["final_span"] = hi - lo + 1
    return lo, hi, diag


def apply(path):
    names, seqs = read_fa(path)
    if len(seqs) < 2:
        return None
    ci = consensus_index(names)
    cons = seqs[ci]
    others = [s for i, s in enumerate(seqs) if i != ci]
    lo, hi, diag = flanked_element_window(cons, others)

    out = []
    for i, s in enumerate(seqs):
        s = s.ljust(len(cons), "-")
        out.append(justify(s, lo, hi))

    order = [ci] + [i for i in range(len(names)) if i != ci]
    with io.open(path, "w", encoding="utf-8") as fh:
        for i in order:
            h, s = names[i], out[i]
            fh.write(u">%s\n" % h)
            for k in range(0, len(s), 80):
                fh.write(s[k:k + 80] + u"\n")
    return diag


def main():
    paths = sys.argv[1:] if len(sys.argv) > 1 else sorted(
        glob.glob(os.path.join("alignments", "*.aln.fa")))
    n = 0
    ext_l = ext_r = 0
    for p in paths:
        d = apply(p)
        if not d:
            continue
        n += 1
        ext_l += d.get("extended_left", 0)
        ext_r += d.get("extended_right", 0)
        el = d.get("extended_left", 0)
        er = d.get("extended_right", 0)
        cap = d.get("capped_right", 0)
        if el or er or d.get("trimmed_left") or d.get("trimmed_right") or cap:
            print("%-40s extend L=%d R=%d  trim L=%d R=%d  capR=%d  span %d->%d"
                  % (os.path.basename(p), el, er,
                     d.get("trimmed_left", 0), d.get("trimmed_right", 0),
                     cap, d.get("old_span", 0), d.get("final_span", 0)))
    print("boundary_justify: %d files" % n)


if __name__ == "__main__":
    main()
