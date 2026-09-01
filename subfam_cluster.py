#!/usr/bin/env python3
"""Cluster SubFam chunk-consensuses - the step Sergei currently does by eye.

The subfam alignment is ~161 chunk-consensuses, each summarising 50 original
sequences. He reads it, picks out groups of similar sequences, and turns each
group into a consensus for the next rescan.

This clusters them instead. It is validated against a known answer: he grouped
these same 161 by hand into t1-t8, so the automated partition can be compared to
his rather than merely admired.

Similarity is ungapped k-mer cosine on both strands - chunk consensuses are
short and of uneven length, so a full pairwise alignment is neither necessary
nor robust here. Clustering is single-pass agglomerative on that similarity.
"""
import sys
from collections import Counter

import numpy as np

K = 6


def read_fa(path):
    names, seqs, cur, buf = [], [], None, []
    for line in open(path):
        line = line.rstrip("\n\r")
        if line.startswith(">"):
            if cur is not None:
                seqs.append("".join(buf))
            cur = line[1:].split()[0]
            names.append(cur)
            buf = []
        else:
            buf.append(line.strip())
    if cur is not None:
        seqs.append("".join(buf))
    return names, seqs


COMP = str.maketrans("ACGTacgtNn-", "TGCAtgcaNn-")


def rc(s):
    return s.translate(COMP)[::-1]


def kmers(s, k=K):
    s = s.upper().replace("-", "")
    c = Counter()
    for i in range(len(s) - k + 1):
        w = s[i:i + k]
        if "N" in w:
            continue
        c[w] += 1
    return c


def vec(c, idx):
    v = np.zeros(len(idx), dtype=np.float32)
    for w, n in c.items():
        j = idx.get(w)
        if j is not None:
            v[j] = n
    nrm = np.linalg.norm(v)
    return v / nrm if nrm else v


def similarity(seqs):
    """Cosine on k-mer profiles, max over the two strands."""
    fwd = [kmers(s) for s in seqs]
    rev = [kmers(rc(s)) for s in seqs]
    vocab = set()
    for c in fwd + rev:
        vocab.update(c)
    idx = {w: i for i, w in enumerate(sorted(vocab))}
    F = np.vstack([vec(c, idx) for c in fwd])
    R = np.vstack([vec(c, idx) for c in rev])
    S = np.maximum(F @ F.T, F @ R.T)
    np.fill_diagonal(S, 1.0)
    return np.maximum(S, S.T)


def agglomerate(S, thr):
    """Average-linkage agglomerative clustering, stopping at `thr`."""
    n = S.shape[0]
    clusters = [[i] for i in range(n)]
    active = list(range(len(clusters)))
    D = S.copy().astype(np.float64)

    while len(active) > 1:
        best, bi, bj = -1.0, None, None
        for a in range(len(active)):
            for b in range(a + 1, len(active)):
                i, j = active[a], active[b]
                m = np.mean(D[np.ix_(clusters[i], clusters[j])])
                if m > best:
                    best, bi, bj = m, i, j
        if best < thr:
            break
        clusters[bi] = clusters[bi] + clusters[bj]
        active.remove(bj)
    return [clusters[i] for i in active]


def tightness(S, members):
    if len(members) < 2:
        return 1.0
    sub = S[np.ix_(members, members)]
    iu = np.triu_indices(len(members), 1)
    return float(np.mean(sub[iu]))


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "subfam/allcons.fa"
    thrs = [float(x) for x in sys.argv[2:]] or [0.55, 0.60, 0.65, 0.70, 0.75]
    names, seqs = read_fa(path)
    print("chunk consensuses: %d   length min=%d med=%d max=%d"
          % (len(seqs), min(map(len, seqs)),
             int(np.median([len(s) for s in seqs])), max(map(len, seqs))))
    S = similarity(seqs)
    iu = np.triu_indices(len(seqs), 1)
    v = S[iu]
    print("pairwise similarity: min=%.3f q1=%.3f median=%.3f q3=%.3f max=%.3f"
          % (v.min(), np.percentile(v, 25), np.median(v),
             np.percentile(v, 75), v.max()))
    print()
    for thr in thrs:
        cl = agglomerate(S, thr)
        cl.sort(key=lambda m: -len(m))
        big = [c for c in cl if len(c) >= 3]
        singles = sum(1 for c in cl if len(c) == 1)
        print("thr %.2f -> %3d clusters   (>=3 members: %2d, singletons: %3d)"
              % (thr, len(cl), len(big), singles))
        for c in big[:10]:
            print("      n=%-3d tightness=%.3f  e.g. %s"
                  % (len(c), tightness(S, c), ", ".join(names[i] for i in c[:3])))
    print()
    print("Sergei's manual pass on this same input produced 8 groups (t1-t8),")
    print("later refined to 6 and then to 14. That is the comparison target.")


if __name__ == "__main__":
    main()
