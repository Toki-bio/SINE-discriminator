#!/usr/bin/env python3
"""Mosaicism by partition congruence, testing Sergei's k-mer proposal.

The current statistic (rank-1 residual on the copy x window identity matrix)
fails on a real mosaic: a copy that is 40 % subfamily A and 60 % B still matches
the A consensus at a level the rank-1 model reads as "an older copy".

Sergei proposed k-mers, 4 to 10, as separate graphs. Two things that could mean,
and they are not equally useful:

  (a) k-mer SPECTRUM per copy - composition without alignment. Useful for
      low-complexity and internal duplication, but it collapses position, and
      mosaicism is entirely about position. Tested below and expected to fail.

  (b) k-mers used to define a LOCAL PARTITION of the copies in each window,
      then asking whether that partition is the SAME from window to window.
      This is the discriminating quantity, and it is what the spec means by
      "permits legitimate subfamily structure - which is global and congruent
      across windows - while flagging block-swapped composites". Subfamily
      structure splits the copies the same way everywhere; a mosaic splits them
      differently in different blocks.

Congruence is measured with the adjusted Rand index between the two-group
partitions of consecutive windows.
"""
import json, glob, os, sys
from collections import defaultdict
import numpy as np
import measure_c as M

WIN, STEP = 40, 20


def kmer_vec(seq, k):
    d = defaultdict(int)
    for i in range(len(seq) - k + 1):
        s = seq[i:i + k]
        if "-" not in s and "N" not in s:
            d[s] += 1
    return d


def partition(block):
    """Split copies in one window into two groups by sequence similarity.
    Leading eigenvector of the copy x copy identity matrix - no k required, and
    equivalent to clustering on shared k-mers for k around the window size."""
    ok = (block != M.GAP)
    same = np.zeros((block.shape[0],) * 2, np.float32)
    for b in range(4):
        m = ((block == b) & ok).astype(np.float32)
        same += m @ m.T
    both = ok.astype(np.float32) @ ok.astype(np.float32).T
    with np.errstate(invalid="ignore", divide="ignore"):
        S = np.where(both >= WIN * 0.4, same / np.maximum(both, 1), np.nan)
    if not np.isfinite(S).any():
        return None
    S = np.where(np.isfinite(S), S, np.nanmean(S))
    Sc = S - S.mean(axis=0, keepdims=True)
    try:
        _, v = np.linalg.eigh(Sc @ Sc.T)
    except np.linalg.LinAlgError:
        return None
    e = v[:, -1]
    return (e > np.median(e)).astype(np.int8)


def adj_rand(a, b):
    n = len(a)
    if n < 8:
        return np.nan
    cm = np.zeros((2, 2))
    for i in range(n):
        cm[a[i], b[i]] += 1
    def c2(x):
        return x * (x - 1) / 2
    sij = c2(cm).sum()
    si = c2(cm.sum(axis=1)).sum()
    sj = c2(cm.sum(axis=0)).sum()
    tot = c2(n)
    exp = si * sj / tot
    mx = (si + sj) / 2
    return float((sij - exp) / (mx - exp)) if mx != exp else np.nan


def analyse(path):
    names, A = M.read_aln(path)
    ci = [i for i, n in enumerate(names) if "CONSENSUS_" in n]
    if not ci:
        return None
    k = ci[0]
    cons = A[k]
    nz = np.where(cons != M.GAP)[0]
    if len(nz) < 120:
        return None
    # index the consensus columns directly: nz[0]:nz[-1] spans the STRETCHED
    # alignment (1.7x longer than the consensus), so a coverage threshold against
    # its width can never be met
    C = np.delete(A, k, axis=0)[:, nz]
    n, L = C.shape

    # keep copies that are reasonably complete, so absence is not read as split
    full = (C != M.GAP).sum(axis=1) > 0.7 * L
    C = C[full]
    if C.shape[0] < 20:
        return None

    parts, pos = [], []
    for s in range(0, L - WIN + 1, STEP):
        p = partition(C[:, s:s + WIN])
        if p is not None:
            parts.append(p)
            pos.append(s + WIN // 2)
    if len(parts) < 4:
        return None

    # (b) congruence between consecutive windows, and between all pairs
    adj = [adj_rand(parts[i], parts[i + 1]) for i in range(len(parts) - 1)]
    allp = [adj_rand(parts[i], parts[j])
            for i in range(len(parts)) for j in range(i + 1, len(parts))]
    adj = [x for x in adj if np.isfinite(x)]
    allp = [x for x in allp if np.isfinite(x)]

    # (a) the k-mer spectrum route, for comparison: how much does per-copy k-mer
    # composition vary between copies, at several k
    seqs = ["".join("ACGT-"[c] for c in row) for row in C[:60]]
    spec = {}
    for kk in (4, 6, 8, 10):
        vecs = [kmer_vec(s.replace("-", ""), kk) for s in seqs]
        keys = sorted({x for v in vecs for x in v})
        if len(keys) < 4:
            continue
        Mx = np.array([[v.get(x, 0) for x in keys] for v in vecs], float)
        Mx = Mx / np.maximum(Mx.sum(axis=1, keepdims=True), 1)
        cmat = np.corrcoef(Mx)
        spec["kmer%d_mean_corr" % kk] = float(np.nanmean(cmat[np.triu_indices(len(vecs), 1)]))
    out = {"set": os.path.basename(path).replace(".aln.fa", ""),
           "n_used": int(C.shape[0]), "n_windows": len(parts),
           "congruence_adjacent": round(float(np.mean(adj)), 4) if adj else None,
           "congruence_allpairs": round(float(np.mean(allp)), 4) if allp else None,
           "congruence_min": round(float(np.min(allp)), 4) if allp else None,
           "window_pos": pos}
    out.update({k2: round(v, 4) for k2, v in spec.items()})
    return out


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "aln_v2"
    files = sorted(glob.glob(os.path.join(d, "*.aln.fa")))
    files = [f for f in files if not os.path.basename(f).startswith("SQ__")]
    res = {}
    for f in files:
        try:
            r = analyse(f)
        except Exception:
            r = None
        if r:
            res[r["set"]] = r
    json.dump(res, open("mosaic_kmer.json", "w"), separators=(",", ":"))

    byc = defaultdict(list)
    for s, v in res.items():
        byc[s.split("__")[0]].append(v)
    print("PARTITION CONGRUENCE - does the copy split stay the same along the element?")
    print("1.0 = identical split everywhere (one family, or clean subfamilies)")
    print("0.0 = the split is unrelated between windows (block-swapped mosaic)\n")
    print("%-12s %5s %12s %12s %10s | %s" % ("class", "sets", "adjacent", "all pairs",
                                             "minimum", "k-mer spectrum corr (4/6/8/10)"))
    print("-" * 100)
    for c in ("POS", "MIXSUBFAM", "NEGMOSAIC", "NEGSPLICE", "NEGCHIM", "MIXED30", "NEGRAND"):
        v = byc.get(c, [])
        if not v:
            continue
        g = lambda k: np.nanmean([x[k] for x in v if x.get(k) is not None])
        ks = " / ".join("%.3f" % g("kmer%d_mean_corr" % kk) if any(
            x.get("kmer%d_mean_corr" % kk) is not None for x in v) else "-"
            for kk in (4, 6, 8, 10))
        print("%-12s %5d %12.3f %12.3f %10.3f | %s"
              % (c, len(v), g("congruence_adjacent"), g("congruence_allpairs"),
                 g("congruence_min"), ks))


if __name__ == "__main__":
    main()
