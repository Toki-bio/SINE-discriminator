#!/usr/bin/env python3
"""Step 1 of the subfam peel loop: can clustering recover his groups at all?

This is the go/no-go. If automated clustering cannot approximately recover the
partition he made by hand from the same chunk consensuses, the rest of the loop
is not worth building.

Input   597 SubFam chunk consensuses from a SINEderella run
Truth   his 8 groups, recovered per chunk by asking which subfamily that
        chunk's 50 members were assigned to (majority vote over the run's own
        assignment_full.tsv)
Output  a clustering, scored against that truth chunk by chunk

Three things this does that the earlier design did not:

- ORIENT FIRST. About half of all families come out reverse-complemented,
  and a family whose chunks are in mixed orientation will not cluster with
  itself - it splits into two groups that look unrelated. His eye flips
  sequences without noticing; clustering cannot.
- REALIGN FIRST. MANUAL section 6.1.1 is explicit that input.clw is not a
  converged alignment, so it is degapped and realigned before anything is
  measured off it.
- Rank candidate groups by ISOLATION, not tightness - his words: "easily
  defineable clusters, cleanly isolated from alignment". A tight cluster can
  sit inside a blur; what he peels first has a clear gap to everything else.
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


def truth_per_chunk(run, subfam_dir):
    """Which of his groups does each chunk belong to? Majority over its members."""
    loc2grp = {}
    p = os.path.join(run, "results", "assignment_full.tsv")
    for line in io.open(p, encoding="utf-8", errors="replace"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 5 or f[0] == "Sequence" or f[4] != "assigned":
            continue
        loc2grp[f[0].replace("@U@", "_")] = f[1]

    out = {}
    for bnk in sorted(glob.glob(os.path.join(subfam_dir, "*.bnk"))):
        tag = os.path.basename(bnk)
        votes = Counter()
        for line in io.open(bnk, encoding="utf-8", errors="replace"):
            if not line.startswith(">"):
                continue
            h = line[1:].strip()
            # headers look like CM115237.1:39539441-39539676(+,-)
            m = re.match(r"([^:]+:\d+-\d+)", h)
            if not m:
                continue
            base = m.group(1)
            for suf in ("(+)", "(-)"):
                if base + suf in loc2grp:
                    votes[loc2grp[base + suf]] += 1
                    break
        if votes:
            top, n = votes.most_common(1)[0]
            out[tag] = (top, n / float(sum(votes.values())))
    return out


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


def single_link(M, thr):
    n = M.shape[0]
    lab = list(range(n))

    def find(x):
        while lab[x] != x:
            lab[x] = lab[lab[x]]
            x = lab[x]
        return x

    for i in range(n):
        for j in range(i + 1, n):
            if M[i, j] >= thr:
                a, b = find(i), find(j)
                if a != b:
                    lab[a] = b
    g = defaultdict(list)
    for i in range(n):
        g[find(i)].append(i)
    return sorted(g.values(), key=len, reverse=True)


def isolation(M, members, others):
    """His criterion: a clean gap between the group and everything else."""
    if len(members) < 2 or not others:
        return 0.0
    within = np.mean([M[a, b] for i, a in enumerate(members) for b in members[i + 1:]])
    between = np.max([M[a, b] for a in members for b in others])
    return float(within - between)


def rand_index(a, b):
    """Agreement between two labellings of the same items."""
    pairs_same_a = pairs_same_b = both = total = 0
    n = len(a)
    for i in range(n):
        for j in range(i + 1, n):
            sa, sb = a[i] == a[j], b[i] == b[j]
            pairs_same_a += sa
            pairs_same_b += sb
            both += sa and sb
            total += 1
    exp = pairs_same_a * pairs_same_b / float(total) if total else 0
    mx = (pairs_same_a + pairs_same_b) / 2.0
    return (both - exp) / (mx - exp) if mx != exp else 1.0


def main():
    run = sys.argv[1]
    work = sys.argv[2]
    os.makedirs(work, exist_ok=True)
    sub = os.path.join(run, "genome.clean_step1", "subfam_input")

    print("1. his groups, per chunk")
    truth = truth_per_chunk(run, sub)
    tally = Counter(g for g, _ in truth.values())
    print("   %d chunks labelled: %s" % (len(truth), dict(tally.most_common())))
    pure = np.mean([p for _, p in truth.values()])
    print("   mean purity of a chunk (fraction of members in its majority group): %.3f" % pure)

    print("2. orient, degap, realign")
    cons = os.path.join(work, "chunks.fa")
    with io.open(cons, "w", encoding="utf-8") as fh:
        for f in sorted(glob.glob(os.path.join(sub, "*.bnk.cons"))):
            n, s = read_fa(f)
            if s and len(s[0].replace(GAP, "")) >= 60:
                fh.write(u">%s\n%s\n" % (os.path.basename(f).replace(".bnk.cons", ""),
                                         s[0].replace(GAP, "").upper()))
    subprocess.run("python3 orient_consensus.py %s %s > %s/orient.log 2>&1"
                   % (cons, os.path.join(work, "chunks.oriented.fa"), work),
                   shell=True, cwd="/staging/tmp/sinedisc")
    ori = os.path.join(work, "chunks.oriented.fa")
    tail = subprocess.run("tail -1 %s/orient.log" % work, shell=True,
                          capture_output=True, text=True).stdout.strip()
    print("   " + tail)

    aln = os.path.join(work, "chunks.aln.fa")
    subprocess.run("mafft --retree 2 --maxiterate 0 --quiet --thread 32 %s > %s"
                   % (ori, aln), shell=True)
    names, seqs = read_fa(aln)
    A = np.array([list(s.upper()) for s in seqs])
    print("   %d chunks aligned, %d columns" % A.shape)

    print("3. cluster and score against his partition")
    M = identity_matrix(A)
    off = M[np.triu_indices(len(names), 1)]
    print("   pairwise identity: median %.3f, 90th %.3f, max %.3f"
          % (np.median(off), np.percentile(off, 90), off.max()))

    key = {n: truth.get(n + ".bnk", ("?", 0))[0] for n in names}
    tvec = [key[n] for n in names]
    best = None
    for thr in (0.60, 0.65, 0.70, 0.75, 0.80, 0.85):
        groups = single_link(M, thr)
        cvec = [0] * len(names)
        for gi, g in enumerate(groups):
            for i in g:
                cvec[i] = gi
        ari = rand_index(tvec, cvec)
        big = [len(g) for g in groups if len(g) >= 5]
        print("   thr %.2f -> %3d clusters (%d with 5+), agreement with his groups %.3f"
              % (thr, len(groups), len(big), ari))
        if best is None or ari > best[1]:
            best = (thr, ari, groups)

    thr, ari, groups = best
    print("   best: identity %.2f, agreement %.3f" % (thr, ari))

    print("4. the cleanly isolated groups, in the order he would peel them")
    allidx = set(range(len(names)))
    rows = []
    for g in groups:
        if len(g) < 5:
            continue
        others = sorted(allidx - set(g))
        iso = isolation(M, g, others)
        lab = Counter(key[names[i]] for i in g).most_common(1)[0]
        rows.append((iso, len(g), lab[0], lab[1] / float(len(g))))
    rows.sort(reverse=True)
    for iso, n, lab, pur in rows[:12]:
        print("   isolation %+.3f  %3d chunks  mostly %-4s (%.0f %% of the cluster)"
              % (iso, n, lab, 100 * pur))

    json.dump({"threshold": thr, "agreement": ari,
               "clusters": [[names[i] for i in g] for g in groups if len(g) >= 5],
               "truth": {n: key[n] for n in names}},
              open(os.path.join(work, "step1.json"), "w"), indent=1)
    print("wrote %s/step1.json" % work)


if __name__ == "__main__":
    main()
