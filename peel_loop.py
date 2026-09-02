#!/usr/bin/env python3
"""The peel loop, his way.

  gather the low hanging fruit - clusters cleanly isolated from the alignment
  remove them
  re-align the remainder when a lot has gone
  repeat until nothing is cleanly separable, leaving only ungroupable noise

Two gates on peeling, not one. His criterion is ISOLATION - a clear gap to
everything outside. But step 1 showed a 143-chunk cluster that was isolated
(+0.213) and yet only 40 % one of his groups: isolated from the rest, and a blur
inside. Peeling that would fuse three of his subfamilies into one consensus. So
a cluster is peeled only if it is BOTH

  isolated   within-group identity minus the best identity to anything outside
  coherent   tight enough inside that it is not itself several groups

An isolated but incoherent cluster is left in, deliberately: re-aligning without
the dominant groups is exactly what should make its internal structure separable
on a later round. That is his stated reason for re-aligning at all.

Scored every round against his own partition, so the loop can be judged rather
than admired.
"""
import io
import json
import os
import subprocess
import sys
from collections import Counter, defaultdict

import numpy as np

GAP = "-"
MIN_CLUSTER = 5
ISO_MIN = 0.10          # a clear gap to everything outside; t7 sits at 0.145
COH_MIN = 0.62          # median identity inside the cluster
REALIGN_FIRST = 100     # his trigger for round 1
REALIGN_LATER = 40      # "less are required for re-alignment in future steps"
MAX_ROUNDS = 8


def read_fa(p):
    names, seqs, cur, buf = [], [], None, []
    for line in io.open(p, encoding="utf-8", errors="replace"):
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


def write_fa(path, names, seqs, degap=True):
    with io.open(path, "w", encoding="utf-8") as fh:
        for n, s in zip(names, seqs):
            fh.write(u">%s\n%s\n" % (n, s.replace(GAP, "") if degap else s))


def identity_matrix(A):
    n = A.shape[0]
    M = np.eye(n)
    for i in range(n):
        ai = A[i]
        for j in range(i + 1, n):
            ok = (ai != GAP) & (A[j] != GAP)
            k = int(ok.sum())
            M[i, j] = M[j, i] = float((ai[ok] == A[j][ok]).mean()) if k >= 40 else 0.0
    return M


def average_link(M, thr):
    """Agglomerate on AVERAGE identity between clusters, not best pair.

    Single linkage chains: t3, t4 and t5 were strung into one 143-chunk cluster
    through intermediates, and because every adjacent pair is highly similar its
    median within-cluster identity came out at 0.983 - higher than any real
    group - while it was only 40 % pure. No within-cluster statistic can see a
    join made one link at a time. Average linkage cannot make that merge,
    because the average across the whole chain is low even when the links are
    high.
    """
    n = M.shape[0]
    clusters = {i: [i] for i in range(n)}
    S = M.astype(float).copy()
    np.fill_diagonal(S, -1.0)
    while len(clusters) > 1:
        keys = sorted(clusters)
        best, pair = -1.0, None
        for x in range(len(keys)):
            for y in range(x + 1, len(keys)):
                a, b = keys[x], keys[y]
                v = S[a, b]
                if v > best:
                    best, pair = v, (a, b)
        if pair is None or best < thr:
            break
        a, b = pair
        na, nb = len(clusters[a]), len(clusters[b])
        for c in list(clusters):
            if c in (a, b):
                continue
            S[a, c] = S[c, a] = (S[a, c] * na + S[b, c] * nb) / float(na + nb)
        clusters[a] = clusters[a] + clusters[b]
        del clusters[b]
        S[b, :] = -1.0
        S[:, b] = -1.0
    return sorted(clusters.values(), key=len, reverse=True)


def splits_further(M, members, thr, step=0.08):
    """Is this one group, or several that happen to sit together?

    His rule is to take the low hanging fruit first and leave the rest for a
    later round, after re-alignment. The test cannot use purity - that is the
    answer we are trying to recover. It has to be intrinsic: re-cluster the
    cluster on its own at a stricter threshold, and if it falls into two or more
    real pieces it is not yet a single group, so defer it.

    This is what keeps t3 and t4 in the alignment for round 2 instead of being
    peeled as one 143-chunk lump in round 1.
    """
    if len(members) < 2 * MIN_CLUSTER:
        return False
    sub = M[np.ix_(members, members)]
    parts = average_link(sub, min(0.97, thr + step))
    big = [g for g in parts if len(g) >= MIN_CLUSTER]
    return len(big) >= 2


def cluster_stats(M, members, others):
    within = [M[a, b] for i, a in enumerate(members) for b in members[i + 1:]]
    coh = float(np.median(within)) if within else 0.0
    best_out = float(np.max([M[a, b] for a in members for b in others])) if others else 0.0
    return coh, float(np.mean(within) - best_out) if within else 0.0


def main():
    step1 = sys.argv[1]          # working dir from peel_step1
    work = sys.argv[2]
    thr = float(sys.argv[3]) if len(sys.argv) > 3 else 0.80
    os.makedirs(work, exist_ok=True)

    truth = json.load(open(os.path.join(step1, "step1.json")))["truth"]
    cur = os.path.join(work, "round0.fa")
    subprocess.run("cp %s %s" % (os.path.join(step1, "chunks.oriented.fa"), cur), shell=True)

    peeled, log = [], []
    for rnd in range(1, MAX_ROUNDS + 1):
        aln = os.path.join(work, "round%d.aln.fa" % rnd)
        subprocess.run("mafft --retree 2 --maxiterate 0 --quiet --thread 32 %s > %s"
                       % (cur, aln), shell=True)
        names, seqs = read_fa(aln)
        if len(names) < MIN_CLUSTER * 2:
            print("round %d: %d chunks left, stopping" % (rnd, len(names)))
            break
        A = np.array([list(s.upper()) for s in seqs])
        M = identity_matrix(A)
        groups = average_link(M, thr)
        allidx = set(range(len(names)))

        take, report = [], []
        for g in groups:
            if len(g) < MIN_CLUSTER:
                continue
            coh, iso = cluster_stats(M, g, sorted(allidx - set(g)))
            labs = Counter(truth.get(names[i], "?") for i in g)
            lab, cnt = labs.most_common(1)[0]
            defer = splits_further(M, g, thr)
            ok = iso >= ISO_MIN and coh >= COH_MIN and not defer
            report.append((iso, coh, len(g), lab, cnt / float(len(g)), ok,
                           "splits" if defer else ""))
            if ok:
                take.append(g)

        print("round %d: %d chunks, %d clusters, %d peelable"
              % (rnd, len(names), len([g for g in groups if len(g) >= MIN_CLUSTER]), len(take)))
        for iso, coh, n, lab, pur, ok, why in sorted(report, reverse=True)[:8]:
            print("    %s iso %+.3f coh %.3f  %3d chunks  mostly %-4s %.0f %% pure %s"
                  % ("PEEL" if ok else "defer", iso, coh, n, lab, 100 * pur, why))

        if not take:
            print("round %d: nothing cleanly separable, stopping" % rnd)
            break

        removed = set()
        for g in take:
            labs = Counter(truth.get(names[i], "?") for i in g)
            lab, cnt = labs.most_common(1)[0]
            peeled.append({"round": rnd, "size": len(g), "majority": lab,
                           "purity": round(cnt / float(len(g)), 3),
                           "members": [names[i] for i in g]})
            removed |= set(g)
        log.append({"round": rnd, "chunks_in": len(names), "peeled": len(removed),
                    "clusters": len(take)})

        keep = [i for i in range(len(names)) if i not in removed]
        cur = os.path.join(work, "round%d.fa" % rnd)
        write_fa(cur, [names[i] for i in keep], [seqs[i] for i in keep])
        trig = REALIGN_FIRST if rnd == 1 else REALIGN_LATER
        print("    peeled %d chunks, %d remain%s"
              % (len(removed), len(keep),
                 "; re-aligning" if len(removed) >= trig else "; below re-align trigger"))
        if len(keep) < MIN_CLUSTER * 2:
            break

    names, seqs = read_fa(cur)
    resid = Counter(truth.get(n, "?") for n in names)
    print()
    print("peeled %d groups over %d rounds" % (len(peeled), len(log)))
    for p in peeled:
        print("   round %d  %3d chunks  mostly %-4s  %.0f %% pure"
              % (p["round"], p["size"], p["majority"], 100 * p["purity"]))
    print("residue: %d chunks %s" % (len(names), dict(resid.most_common())))
    json.dump({"peeled": peeled, "log": log,
               "residue": [n for n in names]},
              open(os.path.join(work, "peel.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
