#!/usr/bin/env python3
"""Apply the chosen flank strategy (v2) to the whole corpus.

De-gap each copy's flanks and push them against the element; leave the element
alignment byte-identical. See HANDOFF.md sections 14 and 14a for why this beats
re-aligning the flanks at any gap penalty.
"""
import os, glob, sys
import numpy as np
import measure_c as M

SRC, DST = "aln_c", "aln_v2"


def justify(path, out):
    names, A = M.read_aln(path)
    ci = [i for i, n in enumerate(names) if "CONSENSUS_" in n]
    if not ci:
        return False
    k = ci[0]
    nz = np.where(A[k] != M.GAP)[0]
    if len(nz) < 40:
        return False
    lo, hi = int(nz[0]), int(nz[-1])
    sym = "ACGT-"
    rows = []
    for i in range(len(names)):
        if i == k:
            rows.append((names[i], "", "".join(sym[c] for c in A[i][lo:hi + 1]), ""))
            continue
        r = A[i]
        l = r[:lo]
        l = l[l != M.GAP]
        rr = r[hi + 1:]
        rr = rr[rr != M.GAP]
        rows.append((names[i], "".join(sym[c] for c in l),
                     "".join(sym[c] for c in r[lo:hi + 1]),
                     "".join(sym[c] for c in rr)))
    wL = max(len(x[1]) for x in rows)
    wR = max(len(x[3]) for x in rows)
    with open(out, "w") as fh:
        for nm, l, e, rr in rows:
            fh.write(">%s\n%s%s%s\n" % (nm, l.rjust(wL, "-"), e, rr.ljust(wR, "-")))
    return True


def main():
    os.makedirs(DST, exist_ok=True)
    files = sorted(glob.glob(os.path.join(SRC, "*.aln.fa")))
    n = 0
    for f in files:
        if justify(f, os.path.join(DST, os.path.basename(f))):
            n += 1
    print("justified %d / %d" % (n, len(files)))


if __name__ == "__main__":
    main()
