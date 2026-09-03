#!/usr/bin/env python3
"""Peel subfamilies on shared diagnostic columns, not on identity.

MANUAL 6.1.6: a subfamily is "a specific, shared, diagnostic pattern of small
indels/SNPs common to that lineage's copies (a synapomorphy) - not the generic
accumulation of private per-copy mutations from ordinary post-insertion decay,
which is noise on top of the subfamily signal."

So identity is the wrong axis. This clusters on features.

Difference from his SINEClusterer, and the reason for it: there, a candidate
group is the EXACT set of sequences sharing one (pos, char), merged afterwards
at Jaccard >= 0.90. On 600 noisy chunk consensuses an exact co-occurrence set
almost never repeats - one sequencing-noise base splits it - so candidates
fragment and only the coarsest structure survives. That is consistent with the
documented result: 2 blobs + 1 small cluster on saq's 9 subfamilies.

Here a group is instead defined by a BLOCK of features that co-occur:

  1. every (pos, char) carried by between MIN_SET and 50% of the pool is a
     candidate feature, gaps excluded (never diagnostic), and columns where one
     base holds >80% of the WHOLE alignment are skipped as non-diagnostic
     (his rule, and his reason: measure it globally, not against the shrinking
     pool, or a clean remainder looks non-diagnostic once earlier groups leave)
  2. features are clustered against each other by Jaccard on their sequence
     sets - features marking the same lineage co-occur
  3. a feature block scores by how exclusive its features are: a member should
     carry most of the block, a non-member almost none
  4. the group is the set of sequences carrying >= CARRY of the block

Because membership is "most of the block", one noisy column no longer expels a
sequence. Then peel the winner, re-align the remainder when enough has gone
(his rule: >=100 in round 1, >=40 later), and repeat.
"""
import glob
import io
import json
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict

import numpy as np

GAP = "-"
MIN_SET = 5          # a feature must mark at least this many sequences
MAX_FRAC = 0.50      # ...and at most this fraction of the pool (his 50% cap)
GLOBAL_CONS = 0.80   # skip columns where one base holds >80% of the whole aln
FEAT_JACCARD = 0.45  # two features belong to the same block above this
MIN_BLOCK = 3        # a block needs this many features
CARRY = 0.60         # a sequence joins if it carries this fraction of the block
EXCL_MIN = 0.35      # required gap between member carriage and non-member carriage
MIN_GROUP = 5
REALIGN_FIRST = 100
REALIGN_LATER = 40
MAX_ROUNDS = 12


def read_fa(p):
    names, seqs, cur, buf = [], [], None, []
    for line in io.open(p, encoding="utf-8", errors="replace"):
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


def write_fa(path, names, seqs, degap=True):
    with io.open(path, "w", encoding="utf-8") as fh:
        for n, s in zip(names, seqs):
            fh.write(u">%s\n%s\n" % (n, s.replace(GAP, "") if degap else s))


def realign(src, dst, threads=32):
    subprocess.run("mafft --retree 2 --maxiterate 0 --quiet --thread %d %s > %s"
                   % (threads, src, dst), shell=True)


def features(A):
    """(pos, char) -> boolean membership vector. GAP is a state, not a skip.

    SINEClusterer skips gap columns outright - "gaps are not diagnostic" - and I
    copied that. It is wrong here, and it is exactly why t3-1 and t3-2 never
    separated.

    MANUAL 6.1.6 defines a subfamily by "a specific, shared, diagnostic pattern
    of small indels/SNPs" - indels named first - and t3 is that case: a ~15 bp
    insertion present in t3-2 and absent in t3-1. Measured on the v4 alignment,
    14 columns sit at 90-97 % gap in t3-1 against 2 % in t3-2, a cleaner
    discriminator than any base-level column. Skipping gaps makes the t3-1 side
    unrepresentable, because its defining feature IS the gap.

    The global-dominance filter still counts BASES only, as his does: a column
    that is mostly gap across the whole alignment but carries a distinctive base
    subset has to survive, or a minority insertion is filtered away before it can
    be used at all.
    """
    n, L = A.shape
    out = []
    for j in range(L):
        col = A[:, j]
        nongap = col != GAP
        n_ng = int(nongap.sum())
        vals, counts = np.unique(col[nongap], return_counts=True)
        if n_ng and counts.max() / float(n_ng) > GLOBAL_CONS:
            continue                                   # non-diagnostic column
        for v in list(vals) + [GAP]:
            m = (col == v)
            c = int(m.sum())
            if c < MIN_SET or c > MAX_FRAC * n:
                continue
            out.append(((j, v), m))
    return out


def blocks(feats):
    """Group features whose sequence sets agree - a block is one synapomorphy set.

    Seeded growth, NOT single linkage. Linking features transitively chains them
    exactly the way single linkage chained Timema t3/t4/t5 into one 143-chunk
    cluster: every adjacent pair is similar, so the whole feature set collapses
    into one blob and no group survives the carriage test. Measured on Timema's
    596 chunks, union-find gave 5 blocks and peeled 111 of 596 sequences.

    Instead: take the widest unused feature as a seed, attach only features that
    agree with THAT SEED, emit the block, remove them, repeat. Every member of a
    block is then similar to a common reference rather than to a neighbour.
    """
    if not feats:
        return []
    M = np.array([m for _, m in feats])
    sizes = M.sum(axis=1).astype(float)
    inter = M.astype(np.int16) @ M.T.astype(np.int16)
    union = sizes[:, None] + sizes[None, :] - inter
    J = np.where(union > 0, inter / np.maximum(union, 1.0), 0.0)

    order = np.argsort(-sizes)
    used = np.zeros(len(feats), dtype=bool)
    out = []
    for s in order:
        if used[s]:
            continue
        mem = np.where((J[s] >= FEAT_JACCARD) & (~used))[0]
        if s not in mem:
            mem = np.append(mem, s)
        used[mem] = True
        if len(mem) >= MIN_BLOCK:
            out.append(mem.tolist())
    return out


def group_from_block(feats, idxs, n, max_frac=MAX_FRAC):
    """Sequences carrying most of the block, plus how exclusive the block is."""
    M = np.array([feats[i][1] for i in idxs])
    carriage = M.mean(axis=0)
    members = np.where(carriage >= CARRY)[0]
    if len(members) < MIN_GROUP or len(members) > max_frac * n:
        return None
    others = np.setdiff1d(np.arange(n), members)
    inside = float(carriage[members].mean())
    outside = float(carriage[others].mean()) if len(others) else 0.0
    excl = inside - outside
    return {"members": members.tolist(), "inside": inside,
            "outside": outside, "excl": excl, "nfeat": len(idxs)}


def main():
    src = sys.argv[1]           # fasta of chunk consensuses (unaligned ok)
    work = sys.argv[2]
    truth_json = sys.argv[3] if len(sys.argv) > 3 else None
    os.makedirs(work, exist_ok=True)
    truth = json.load(open(truth_json)) if truth_json and os.path.exists(truth_json) else {}

    cur = os.path.join(work, "round0.fa")
    subprocess.run("cp %s %s" % (src, cur), shell=True)

    peeled, log = [], []
    for rnd in range(1, MAX_ROUNDS + 1):
        aln = os.path.join(work, "round%d.aln.fa" % rnd)
        realign(cur, aln)
        names, seqs = read_fa(aln)
        if len(names) < MIN_GROUP * 2:
            print("round %d: %d left, stopping" % (rnd, len(names)))
            break
        A = np.array([list(s.upper()) for s in seqs])
        n = len(names)

        fe = features(A)
        bl = blocks(fe)

        def collect(max_frac):
            out = []
            for idxs in bl:
                g = group_from_block(fe, idxs, n, max_frac)
                if g and g["excl"] >= EXCL_MIN:
                    out.append(g)
            out.sort(key=lambda g: (-g["excl"], -g["nfeat"]))
            return out

        cands = collect(MAX_FRAC)
        relaxed = ""
        if not cands:
            # his rule: retry with the upper bound relaxed rather than give up -
            # a genuinely dominant subfamily can be most of the remaining pool
            cands = collect(0.90)
            if cands:
                relaxed = "  [relaxed upper bound]"

        print("round %d: %d seqs, %d columns, %d features, %d blocks, %d candidate groups%s"
              % (rnd, n, A.shape[1], len(fe), len(bl), len(cands), relaxed))
        if not cands:
            print("  nothing separable, stopping")
            break

        # peel non-overlapping candidates, best first
        taken, used = [], set()
        for g in cands:
            ms = set(g["members"])
            if ms & used:
                continue
            taken.append(g)
            used |= ms

        for g in taken:
            labs = Counter(truth.get(names[i], "?") for i in g["members"])
            lab, cnt = labs.most_common(1)[0]
            pur = cnt / float(len(g["members"]))
            print("  PEEL %3d seqs  %2d feats  excl %.3f (in %.2f / out %.2f)  mostly %-6s %.0f%% pure"
                  % (len(g["members"]), g["nfeat"], g["excl"], g["inside"], g["outside"], lab, 100 * pur))
            peeled.append({"round": rnd, "size": len(g["members"]), "majority": lab,
                           "purity": round(pur, 3), "nfeat": g["nfeat"],
                           "excl": round(g["excl"], 3),
                           "members": [names[i] for i in g["members"]]})

        log.append({"round": rnd, "in": n, "peeled": len(used), "groups": len(taken)})
        keep = [i for i in range(n) if i not in used]
        cur = os.path.join(work, "round%d.fa" % rnd)
        write_fa(cur, [names[i] for i in keep], [seqs[i] for i in keep])
        trig = REALIGN_FIRST if rnd == 1 else REALIGN_LATER
        print("  peeled %d, %d remain%s" % (len(used), len(keep),
              "; re-aligning" if len(used) >= trig else "; below re-align trigger"))
        if len(keep) < MIN_GROUP * 2:
            break

    names, _ = read_fa(cur)
    resid = Counter(truth.get(nm, "?") for nm in names)
    print("\npeeled %d groups over %d rounds" % (len(peeled), len(log)))
    if truth:
        tot = sum(p["size"] for p in peeled)
        wsum = sum(p["size"] * p["purity"] for p in peeled)
        print("weighted mean purity: %.3f over %d seqs" % (wsum / tot if tot else 0, tot))
        best = {}
        for p in peeled:
            if p["majority"] not in best or p["size"] > best[p["majority"]]["size"]:
                best[p["majority"]] = p
        print("\n%-8s %6s %6s" % ("group", "size", "purity"))
        for k in sorted(best):
            print("%-8s %6d %6.2f" % (k, best[k]["size"], best[k]["purity"]))
    print("residue: %d %s" % (len(names), dict(resid.most_common())))
    json.dump({"peeled": peeled, "log": log, "residue": names},
              open(os.path.join(work, "peel_features.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
