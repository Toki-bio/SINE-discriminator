#!/usr/bin/env python3
"""Turn a cluster of SubFam chunks into a candidate consensus - from the MEMBERS.

Sergei's point, and the reason this is not just clustering:

  "the algorithm should not create consensus from these 'consensi of 50s' used
   in subfam but look back into those original 50sequence chunks used for
   subfam consensuses"

A consensus of chunk-consensuses inherits each chunk's smoothing and loses the
real variation. The members are on disk (input_NNN.bnk, 50 sequences each), so
there is no reason to accept that loss.

For each cluster: gather the union of member sequences from its chunks, align
them with MAFFT, and take a plurality consensus over columns present in enough
copies. Then blastn the result against the curated tim consensuses to see
whether the automated grouping recovered what Sergei found by eye.
"""
import os
import subprocess
import sys
import tempfile
from collections import Counter

import numpy as np

sys.path.insert(0, ".")
import subfam_cluster as SC

MIN_COV = 0.50        # a consensus column needs this fraction of copies present
MAX_MEMBERS = 300     # cap per cluster; MAFFT on thousands is slow and needless


def read_fa(path):
    return SC.read_fa(path)


def members_of(chunks, workdir):
    names, seqs = [], []
    for c in chunks:
        p = os.path.join(workdir, c)
        if not os.path.exists(p):
            continue
        n, s = read_fa(p)
        names.extend(n)
        seqs.extend(s)
    return names, seqs


def mafft(seqs, names):
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        for n, s in zip(names, seqs):
            fh.write(">%s\n%s\n" % (n, s))
        inp = fh.name
    out = inp + ".aln"
    cmd = ["mafft", "--retree", "2", "--maxiterate", "0",
           "--adjustdirection", "--quiet", "--thread", "4", inp]
    try:
        with open(out, "w") as oh:
            subprocess.run(cmd, stdout=oh, stderr=subprocess.DEVNULL, check=True)
    except Exception:
        os.unlink(inp)
        return None, None
    n, s = read_fa(out)
    os.unlink(inp)
    os.unlink(out)
    return n, s


def consensus(seqs, min_cov=MIN_COV):
    if not seqs:
        return ""
    L = max(len(s) for s in seqs)
    seqs = [s.upper().ljust(L, "-") for s in seqs]
    n = len(seqs)
    cols = []
    for i in range(L):
        col = [s[i] for s in seqs]
        present = [c for c in col if c in "ACGT"]
        if len(present) < min_cov * n:
            continue
        cols.append(Counter(present).most_common(1)[0][0])
    return "".join(cols)


def main():
    workdir = sys.argv[1] if len(sys.argv) > 1 else "subfam"
    thr = float(sys.argv[2]) if len(sys.argv) > 2 else 0.40
    minsize = int(sys.argv[3]) if len(sys.argv) > 3 else 3

    names, seqs = read_fa(os.path.join(workdir, "allcons.fa"))
    S = SC.similarity(seqs)
    cl = SC.agglomerate(S, thr)
    cl = [c for c in cl if len(c) >= minsize]
    cl.sort(key=lambda m: (-SC.tightness(S, m) * len(m)))

    out = []
    print("clusters at thr %.2f with >=%d chunks: %d" % (thr, minsize, len(cl)))
    for gi, members in enumerate(cl, 1):
        chunks = [names[i] for i in members]
        mn, ms = members_of(chunks, workdir)
        if len(ms) > MAX_MEMBERS:
            step = len(ms) / float(MAX_MEMBERS)
            keep = [int(i * step) for i in range(MAX_MEMBERS)]
            mn = [mn[i] for i in keep]
            ms = [ms[i] for i in keep]
        an, asq = mafft(ms, mn)
        if asq is None:
            print("  g%02d  MAFFT failed" % gi)
            continue
        cons = consensus(asq)
        print("  g%02d  chunks=%-3d members=%-4d tight=%.3f  consensus=%d bp"
              % (gi, len(chunks), len(ms), SC.tightness(S, members), len(cons)))
        out.append(("g%02d" % gi, cons, chunks, len(ms)))

    with open(os.path.join(workdir, "auto_consensus.fa"), "w") as fh:
        for gid, cons, chunks, nm in out:
            if len(cons) >= 60:
                fh.write(">%s n_chunks=%d n_members=%d\n%s\n"
                         % (gid, len(chunks), nm, cons))
    with open(os.path.join(workdir, "auto_groups.tsv"), "w") as fh:
        fh.write("group\tn_chunks\tn_members\tcons_bp\tchunks\n")
        for gid, cons, chunks, nm in out:
            fh.write("%s\t%d\t%d\t%d\t%s\n"
                     % (gid, len(chunks), nm, len(cons), ",".join(chunks)))
    print("\nwrote %s/auto_consensus.fa and auto_groups.tsv"
          % workdir)


if __name__ == "__main__":
    main()
