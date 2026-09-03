#!/usr/bin/env python3
"""Cut the MAFFT 6-mer guide tree, the way MSA-viewer's clusterByGuideTree does.

The mean 6-mer distance failed on t3, which he separated at a glance. That was my
error of measurement, not his of judgement: he said *sorting*, and MAFFT's order
comes from a tree that joins nearest neighbours. A group mean averages away
exactly the local structure the tree is built from — a heterogeneous group with a
tight sub-lineage looks bad by mean and fine by tree.

So: build the tree on MAFFT's own distance, cut it into k groups, and report

  cut height        the merge height at which the cut sits
  next merge        the height the next merge would have had

which is the viewer's own readout: "The next merge would have joined groups at
distance X, so a larger number of groups splits them further." That gap is a
ready-made isolation measure, and better than the one I invented for the peel.

Linkage is average (UPGMA), which is what MAFFT's guide tree uses and what the
peel already learned to use — single linkage chained t3/t4/t5 into one 143-chunk
cluster whose within-cluster identity was the highest of any while being 40 %
pure.

Usage:  tree_cut.py <fasta with GROUP__ prefixed headers> [max_k]
"""
import io
import sys

import numpy as np

sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import mafft_dist as M


def read_fa(path):
    names, seqs, cur, buf = [], [], None, []
    for line in io.open(path, encoding="utf-8", errors="replace"):
        line = line.rstrip()
        if line.startswith(">"):
            if cur is not None:
                seqs.append("".join(buf))
            cur = line[1:]
            names.append(cur)
            buf = []
        else:
            buf.append(line.strip().replace("-", ""))
    if cur is not None:
        seqs.append("".join(buf))
    return names, seqs


def matrix(seqs):
    pre = [M.kmers(s) for s in seqs]
    n = len(seqs)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = M.distance(None, None, pre[i][0], pre[j][0], pre[i][1], pre[j][1])
            D[i, j] = D[j, i] = d
    return D


def upgma(D):
    """Average-linkage agglomeration. Returns merges as (a, b, height, size)."""
    n = D.shape[0]
    act = {i: [i] for i in range(n)}
    S = D.astype(float).copy()
    np.fill_diagonal(S, np.inf)
    merges = []
    while len(act) > 1:
        ks = sorted(act)
        best, pair = np.inf, None
        for x in range(len(ks)):
            for y in range(x + 1, len(ks)):
                a, b = ks[x], ks[y]
                if S[a, b] < best:
                    best, pair = S[a, b], (a, b)
        a, b = pair
        na, nb = len(act[a]), len(act[b])
        for c in list(act):
            if c in (a, b):
                continue
            S[a, c] = S[c, a] = (S[a, c] * na + S[b, c] * nb) / float(na + nb)
        act[a] = act[a] + act[b]
        del act[b]
        S[b, :] = np.inf
        S[:, b] = np.inf
        merges.append((a, b, best, len(act[a])))
    return merges


def cut(merges, n, k):
    """Groups after undoing the last k-1 merges; plus cut and next-merge heights."""
    parent = {i: i for i in range(n)}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    use = merges[:len(merges) - (k - 1)] if k > 1 else merges
    for a, b, h, _ in use:
        parent[find(b)] = find(a)
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    cut_h = use[-1][2] if use else 0.0
    nxt = merges[len(use)][2] if len(use) < len(merges) else None
    return list(groups.values()), cut_h, nxt


def main():
    path = sys.argv[1]
    max_k = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    names, seqs = read_fa(path)
    keep = [i for i, h in enumerate(names) if not h.startswith("CONSENSUS")]
    names = [names[i] for i in keep]
    seqs = [seqs[i] for i in keep]
    truth = [h.split("__")[0] for h in names]
    n = len(seqs)
    print("%d sequences, %d groups in truth: %s"
          % (n, len(set(truth)), sorted(set(truth))))
    D = matrix(seqs)
    merges = upgma(D)
    print("\n%3s %10s %10s %8s  %s" % ("k", "cut", "next", "gap", "composition"))
    for k in range(2, max_k + 1):
        groups, ch, nx = cut(merges, n, k)
        groups.sort(key=len, reverse=True)
        comp = []
        correct = 0
        for g in groups:
            c = {}
            for i in g:
                c[truth[i]] = c.get(truth[i], 0) + 1
            top = max(c, key=c.get)
            correct += c[top]
            comp.append("%s:%d/%d" % (top, c[top], len(g)))
        gap = (nx - ch) if nx is not None else float("nan")
        print("%3d %10.4f %10.4f %8.4f  purity %.3f  %s"
              % (k, ch, nx if nx is not None else float("nan"), gap,
                 correct / float(n), " ".join(comp[:6])))


if __name__ == "__main__":
    main()
