#!/usr/bin/env python3
"""Replace row-0 consensus with per-column majority from copies in this MSA.

step8a pastes a pre-MAFFT conse sequence as row 0; after MAFFT + --adjustdirection
it often no longer matches copy columns.  Rebuild row 0 from the copies in the
same file, then boundary_justify adjusts display case/flanks.

Universal rule (every species / subfamily — no per-set tuning):

  WHERE — element span comes from the alignment geometry already in the file:
          row-0 uppercase columns (border loop BED + step8a extract + any prior
          justify).  Re-walking copy columns to choose span pulls in MAFFT
          column noise (oma_SINE18: an A-rich block left of the real element).

  WHAT  — letters inside that span are copy-column majority, included only when
          CONS_EDGE=0.50 (boundary.py).  Columns at 0.45-0.49 are elevated but
          ambiguous vs random DNA (~0.25-0.30) and must not become consensus
          bases (oma_SINE10 leading T; oma_SINE18 trailing ATT).

  TRIM  — outermost letters still failing CONS_EDGE are dropped (handles internal
          gaps at the 3 prime edge).

Run after step8a MAFFT, before boundary_justify (or re-run justify after).
"""
import glob
import io
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import boundary as B
from fix_alignments import consensus_index, read_fa

GAPS = set("-.")
CONS_EDGE = B.CONS_EDGE


def element_rows(rows):
    """Uppercase element only; lowercase flanks treated as gaps for voting."""
    out = []
    for s in rows:
        out.append("".join(c.upper() if c.isupper() else "-" for c in s))
    return out


def majority_base(rows, j):
    bases = [r[j] for r in rows if j < len(r) and r[j] not in GAPS]
    if not bases:
        return None, 0.0, 0.0
    occ = len(bases) / float(len(rows))
    top = Counter(bases).most_common(1)[0][0]
    c = bases.count(top) / float(len(bases))
    return top, occ, c


def consensus_span(seqs, ci, others):
    """Row-0 uppercase span; fall back to copy-supported window only if empty."""
    upper = [i for i, c in enumerate(seqs[ci]) if c.isupper()]
    if upper:
        return upper[0], upper[-1]
    seed = []
    for j in range(len(seqs[ci])):
        if B.column_supported(others, j, CONS_EDGE):
            b, _, _ = majority_base(others, j)
            seed.append(b or "-")
        else:
            seed.append("-")
    seed = "".join(seed)
    lo, hi, _ = B.element_window(seed, others, cons_min=CONS_EDGE)
    return lo, hi


def apply(path):
    names, seqs = read_fa(path)
    if len(seqs) < 2:
        return None
    ci = consensus_index(names)
    length = max(len(s) for s in seqs)
    seqs = [s.ljust(length, "-") for s in seqs]
    others = element_rows([seqs[i] for i in range(len(seqs)) if i != ci])
    lo, hi = consensus_span(seqs, ci, others)

    new = list("-" * length)
    agree = cols = 0
    old_elem = "".join(c for c in seqs[ci] if c.isupper())
    for j in range(lo, hi + 1):
        cols += 1
        if B.column_supported(others, j, CONS_EDGE):
            b, _, _ = majority_base(others, j)
            if b:
                new[j] = b
                if seqs[ci][j].upper() == b:
                    agree += 1
    lo, hi = B.trim_consensus_edges(new, others, lo, hi, CONS_EDGE)
    seqs[ci] = "".join(new)

    order = [ci] + [i for i in range(len(names)) if i != ci]
    with io.open(path, "w", encoding="utf-8") as fh:
        for i in order:
            h, s = names[i], seqs[i]
            fh.write(u">%s\n" % h)
            for k in range(0, len(s), 80):
                fh.write(s[k:k + 80] + u"\n")

    new_elem = "".join(c for c in seqs[ci] if c.isupper())
    return {
        "lo": lo, "hi": hi,
        "old_elem_len": len(old_elem),
        "elem_len": len(new_elem),
        "elem_cols": cols,
        "col_agree_before": agree,
    }


def main():
    paths = sys.argv[1:] if len(sys.argv) > 1 else sorted(
        glob.glob(os.path.join("alignments", "*.aln.fa")))
    n = 0
    for p in paths:
        d = apply(p)
        if not d:
            continue
        n += 1
        ch = "" if d["old_elem_len"] == d["elem_len"] else (
            "  TRIM" if d["elem_len"] < d["old_elem_len"] else "  EXT")
        print("%-40s elem %d->%d%s  agree %d/%d"
              % (os.path.basename(p), d["old_elem_len"], d["elem_len"], ch,
                 d["col_agree_before"], d["elem_cols"]))
    print("rebuild_consensus_row: %d files" % n)


if __name__ == "__main__":
    main()
