#!/usr/bin/env python3
"""A weighted verdict with named sub-cases, instead of a binary call.

Two things the earlier builds got wrong, both raised by Sergei:

1. A single accept/reject throws away most of what was measured. The verdict is
   now a weighted score over the evidence groups, with every component visible so
   a disagreement can be traced to the variable that caused it.

2. "Negative" was one bucket. It is not one thing. A set can fail because there
   is no element at all, or because there IS an element and something else is
   also true - copies drawn from several structurally different subfamilies, or
   copies sitting inside other repeats. Those are different pieces of future
   work, so they get different names.

**A core of real-looking copies is always reported, however contaminated the set
is.** Twenty or thirty genuine copies inside a heavily polluted candidate are
worth knowing about - they are where the next clean family comes from.

The weights are stated, not fitted. Fitting now would fit the synthetic
negatives, which is the thing this project has been careful not to do.
"""
import json, glob, os, sys
import numpy as np
import measure_c as M

W = {                       # evidence group -> weight
    "element": 0.45,        # is there an element supported by the copies
    "homogeneity": 0.25,    # is it ONE element, uniform across copies
    "uniqueness": 0.20,     # are the copies in unique genomic locations
    "insertion": 0.10,      # one-sided bonuses: TSDs and the insertion site
}
SUPPORT = 0.55              # per-copy identity above which a copy supports the consensus
FLANK_SHARE = 0.55          # ungapped flank identity above which two copies share a flank
BG = 0.25                   # unrelated DNA


def sat(x, lo, hi):
    return float(min(1.0, max(0.0, (x - lo) / (hi - lo))))


def parts(path):
    names, A = M.read_aln(path)
    ci = [i for i, n in enumerate(names) if "CONSENSUS_" in n]
    if not ci:
        return None
    k = ci[0]
    cons = A[k]
    nz = np.where(cons != M.GAP)[0]
    if len(nz) < 60:
        return None
    lo, hi = int(nz[0]), int(nz[-1])
    idx = [i for i in range(len(names)) if i != k]
    C = A[idx]
    el = C[:, lo:hi + 1]
    pres = el != M.GAP
    agree = (el == cons[lo:hi + 1][None, :]) & pres
    with np.errstate(invalid="ignore"):
        ident = agree.sum(axis=1) / np.maximum(pres.sum(axis=1), 1)
    ident[pres.sum(axis=1) < 30] = 0.0
    lefts, rights = [], []
    for r in C:
        l = r[:lo]
        lefts.append(l[l != M.GAP][::-1])
        rr = r[hi + 1:]
        rights.append(rr[rr != M.GAP])
    return {"names": [names[i] for i in idx], "el": el, "pres": pres, "ident": ident,
            "lefts": lefts, "rights": rights, "cons_bp": len(nz), "n": len(idx)}


def flank_sharing(seqs, minlen=35):
    """Copy x copy flank identity, ungapped, anchored at each copy's own edge.
    Above background means two copies sit in the same genomic context - nested in
    a host repeat, inside a duplication, or in a satellite array."""
    n = len(seqs)
    L = min(70, max((len(s) for s in seqs), default=0))
    if L < minlen or n < 4:
        return np.zeros((n, n)), np.zeros(n, bool)
    Mx = np.full((n, L), 4, dtype=np.int8)
    for i, s in enumerate(seqs):
        m = min(L, len(s))
        Mx[i, :m] = s[:m]
    ok = Mx != 4
    same = np.zeros((n, n))
    for b in range(4):
        m = (Mx == b).astype(np.float32)
        same += m @ m.T
    both = ok.astype(np.float32) @ ok.astype(np.float32).T
    with np.errstate(invalid="ignore", divide="ignore"):
        F = np.where(both >= minlen, same / np.maximum(both, 1), 0.0)
    np.fill_diagonal(F, 0.0)
    return F, (F > FLANK_SHARE).any(axis=1)


def subfamily_split(el, pres, ident, min_grp=12):
    """Do the supported copies fall into structurally different groups?

    Copy x copy identity inside the element, split on the leading eigenvector.
    A real family varies smoothly with age; two subfamilies with different
    structure separate into blocks whose between-group identity is clearly below
    their within-group identity."""
    sel = np.where(ident >= SUPPORT)[0]
    if len(sel) < 2 * min_grp:
        return None
    E = el[sel]
    P = pres[sel]
    ok = P.astype(np.float32)
    same = np.zeros((len(sel), len(sel)), np.float32)
    for b in range(4):
        m = ((E == b) & P).astype(np.float32)
        same += m @ m.T
    both = ok @ ok.T
    with np.errstate(invalid="ignore", divide="ignore"):
        S = np.where(both >= 40, same / np.maximum(both, 1), np.nan)
    if not np.isfinite(S).any():
        return None
    S = np.where(np.isfinite(S), S, np.nanmean(S))
    Sc = S - S.mean(axis=0, keepdims=True)
    try:
        w, v = np.linalg.eigh(Sc @ Sc.T)
    except np.linalg.LinAlgError:
        return None
    g = v[:, -1] > np.median(v[:, -1])
    if g.sum() < min_grp or (~g).sum() < min_grp:
        return None
    within = np.nanmean([S[np.ix_(g, g)].mean(), S[np.ix_(~g, ~g)].mean()])
    between = S[np.ix_(g, ~g)].mean()
    return {"gap": float(within - between), "sizes": [int(g.sum()), int((~g).sum())],
            "within": float(within), "between": float(between),
            "members": [int(sel[i]) for i in np.where(g)[0]]}


def verdict(path):
    p = parts(path)
    if p is None:
        return None
    n, ident = p["n"], p["ident"]
    supported = ident >= SUPPORT
    n_sup = int(supported.sum())

    FL, sharedL = flank_sharing(p["lefts"])
    FR, sharedR = flank_sharing(p["rights"])
    shared = sharedL | sharedR
    n_shared = int(shared.sum())
    global_elev = float(np.median(np.concatenate([FL[FL > 0], FR[FR > 0]]))) \
        if (FL > 0).any() or (FR > 0).any() else 0.0

    # a "core" copy: supports the consensus AND sits in unique genomic sequence
    core = supported & ~shared
    n_core = int(core.sum())

    sub = subfamily_split(p["el"], p["pres"], ident)

    elen = (p["el"] != M.GAP).sum(axis=1).astype(float)
    # length spread over the SUPPORTED copies, not the core: excluding copies for
    # sharing flanks would hide truncation that the set really has
    cv_core = float(np.std(elen[supported]) / max(1e-9, np.mean(elen[supported])))         if n_sup >= 8 else 1.0
    id_core = float(np.median(ident[core])) if n_core else 0.0

    fl_id = []
    for seqs in (p["lefts"], p["rights"]):
        F, _ = flank_sharing(seqs)
        v = F[np.triu_indices(len(seqs), 1)]
        v = v[v > 0]
        if len(v):
            fl_id.append(float(np.median(v)))
    flank_bg = float(np.mean(fl_id)) if fl_id else BG

    tsd = [M.find_tsd(p["lefts"][i][::-1], p["rights"][i])
           for i in np.where(core)[0][:120]]
    tsd_frac = float(np.mean([t > 0 for t in tsd])) if tsd else 0.0

    # ---- evidence groups, each 0-1 --------------------------------------
    # The cliff is measured over ALL copies, not over the core. Measuring it on
    # the copies that were selected for matching the consensus is circular - it
    # is how random loci first scored 0.67 here, and the same trap that made the
    # gap-based pruning rule manufacture families out of noise.
    # mean, not median: with 30 % contamination the median is still a real copy,
    # so a median cliff cannot see contamination at all. The score must describe
    # the set as submitted; recoverability is carried by n_core and the flags.
    id_all = float(np.mean(ident))
    g_elem = (sat(id_all - flank_bg, 0.08, 0.55) ** 0.5) *              (sat(n_sup / max(1, n), 0.15, 0.98) ** 0.5)
    g_homog = np.mean([1 - sat(cv_core, 0.05, 0.25),
                       1 - (sat(sub["gap"], 0.02, 0.14) if sub else 0.0)])
    g_uniq = np.mean([1 - sat(n_shared / max(1, n), 0.05, 0.60),
                      1 - sat(global_elev - BG, 0.05, 0.30)])
    g_ins = sat(tsd_frac, 0.10, 0.60)

    # Geometric weighting, not a weighted sum: a set with no element is not a
    # family however unique its flanks are, and a sum lets one strong group mask
    # a fatal weakness in another. The insertion evidence is a one-sided bonus,
    # as the spec requires - it can raise a borderline score, never reject.
    eps = 1e-6
    base = ((g_elem + eps) ** W["element"] * (g_homog + eps) ** W["homogeneity"] *
            (g_uniq + eps) ** W["uniqueness"]) ** (1.0 / (W["element"] +
            W["homogeneity"] + W["uniqueness"]))
    score = 100 * min(1.0, base * (1 + W["insertion"] * g_ins))

    # ---- sub-cases, not one "negative" bucket ---------------------------
    flags = []
    if g_elem < 0.25:
        flags.append({"code": "NO_ELEMENT",
                      "text": "No element: too few copies support the consensus above background."})
    if n_core >= 8:
        if (n - n_sup) / max(1, n) > 0.08:
            flags.append({"code": "CONTAMINATED", "n": int(n - n_sup),
                          "text": "%d of %d copies do not support the consensus; prune them and "
                                  "re-score (see prune.py)." % (n - n_sup, n)})
        if sub and sub["gap"] > 0.03:
            flags.append({"code": "SUBFAMILY_STRUCTURE", "n": sub["sizes"],
                          "text": "Supported copies split into structurally different groups of "
                                  "%d and %d (within-group identity %.2f, between %.2f). Separate "
                                  "them and treat each as its own family."
                                  % (sub["sizes"][0], sub["sizes"][1], sub["within"], sub["between"])})
        if n_shared / max(1, n) > 0.15:
            kind = ("satellite or segmental duplication" if global_elev - BG > 0.15
                    else "another repeat or a duplication")
            flags.append({"code": "NESTED_COPIES", "n": n_shared,
                          "text": "%d of %d copies share flanking sequence with each other — they "
                                  "sit inside %s. Set them aside; they do not disqualify the "
                                  "family." % (n_shared, n, kind)})
    if n_core >= 20 and score < 55:
        flags.append({"code": "RECOVERABLE_CORE", "n": n_core,
                      "text": "Despite the problems above, %d copies look like genuine elements. "
                              "That is enough to build a clean family from — worth following up "
                              "even though the set as submitted scores poorly." % n_core})

    return {"set": os.path.basename(path).replace(".aln.fa", ""),
            "score": round(float(score), 1),
            "groups": {"element": round(float(g_elem), 3),
                       "homogeneity": round(float(g_homog), 3),
                       "uniqueness": round(float(g_uniq), 3),
                       "insertion": round(float(g_ins), 3)},
            "weights": W,
            "n": n, "n_supported": n_sup, "n_core": n_core, "n_shared": n_shared,
            "core_identity": round(id_core, 3), "all_identity": round(float(np.median(ident)), 3), "flank_bg": round(flank_bg, 3),
            "core_len_cv": round(cv_core, 3), "tsd_frac": round(tsd_frac, 3),
            "subfamily": sub and {k: v for k, v in sub.items() if k != "members"},
            "flags": flags}


def _one(f):
    try:
        return verdict(f)
    except Exception as exc:
        return {"set": os.path.basename(f).replace(".aln.fa", ""), "error": str(exc)}


def main():
    import multiprocessing as mp
    d = sys.argv[1] if len(sys.argv) > 1 else "aln_v2"
    out = sys.argv[2] if len(sys.argv) > 2 else "verdicts.json"
    files = sorted(glob.glob(os.path.join(d, "*.aln.fa")))
    files = [f for f in files if not os.path.basename(f).startswith("SQ__")]
    res = {}
    with mp.Pool(max(1, mp.cpu_count() - 2)) as pool:
        for v in pool.imap_unordered(_one, files, chunksize=4):
            if v:
                res[v["set"]] = v
    json.dump(res, open(out, "w"), separators=(",", ":"))

    from collections import defaultdict
    byc = defaultdict(list)
    for s, v in res.items():
        if "error" not in v:
            byc[s.split("__")[0]].append(v)
    print("%-12s %5s %7s   %-38s %s" % ("class", "sets", "score", "evidence groups", "n_core"))
    print("-" * 86)
    for c in ("POS", "MIXED10", "MIXED30", "NEGJITTER", "NEGSPLICE", "NEGTRUNC5", "NEGRAND"):
        v = byc.get(c, [])
        if not v:
            continue
        g = lambda k: np.mean([x["groups"][k] for x in v])
        print("%-12s %5d %7.1f   el %.2f  hom %.2f  uniq %.2f  ins %.2f   %5.1f"
              % (c, len(v), np.mean([x["score"] for x in v]),
                 g("element"), g("homogeneity"), g("uniqueness"), g("insertion"),
                 np.mean([x["n_core"] for x in v])))
    print("\nflag frequency by class")
    for c in ("POS", "MIXED10", "MIXED30", "NEGJITTER", "NEGSPLICE", "NEGTRUNC5", "NEGRAND"):
        v = byc.get(c, [])
        if not v:
            continue
        cnt = defaultdict(int)
        for x in v:
            for f in x["flags"]:
                cnt[f["code"]] += 1
        print("  %-12s %s" % (c, dict(cnt) or "clean"))


if __name__ == "__main__":
    main()
