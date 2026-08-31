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
SUPPORT = 0.55              # legacy constant, kept for reference only


def turbulence(el, pres):
    """How unsettled is the element along its length?

    Sergei on hedgehog e2-3: "can be SINE, but very turbulent at many places,
    should be regarded as SINE but with caution and need manual reinspection."
    Measured as the fraction of element columns whose identity falls more than
    0.25 below the element's own median. A clean family sits at 0.06-0.08;
    e2-3 is 0.18.

    This does not change the SINE verdict. It marks a set as one a person should
    look at, which is a different thing from a set being wrong.
    """
    idc = []
    for c in range(el.shape[1]):
        b = el[:, c][pres[:, c]]
        if len(b) >= 8:
            cnt = np.bincount(b, minlength=4)[:4].astype(float)
            idc.append((cnt * (cnt - 1)).sum() / (len(b) * (len(b) - 1)))
        else:
            idc.append(np.nan)
    idc = np.array(idc, float)
    if not np.isfinite(idc).sum() > 40:
        return None
    med = float(np.nanmedian(idc))
    return float(np.nanmean(idc < med - 0.25))


def contamination_split(ident, bg, min_gap=0.10, min_frac=0.05):
    """Contamination is a SEPARATE MODE, not a low tail.

    Asking "how many copies fall below a line" penalises any diverged family:
    a clean family at 30 % divergence has copies at identity 0.47-0.77 and every
    absolute or proportional cut clips its lower half. But that distribution is
    unimodal - largest internal gap 0.036 - while a genuinely contaminated set is
    bimodal, largest gap 0.247. The shape is divergence-independent; the location
    is not.

    Returns the identity cut separating the two modes, or None if the
    distribution is unimodal, i.e. no contamination.
    """
    d = np.sort(ident)
    n = len(d)
    if n < 20:
        return None
    lo, hi = max(1, int(n * 0.03)), min(n - 1, int(n * 0.97))
    if hi - lo < 3:
        return None
    gaps = np.diff(d[lo:hi + 1])
    j = int(np.argmax(gaps))
    gap = float(gaps[j])
    cut = float((d[lo + j] + d[lo + j + 1]) / 2)
    below = int(np.sum(ident < cut))
    if gap < min_gap or below < max(3, min_frac * n):
        return None                       # unimodal: a diverged family, not a mixture
    if np.median(d[d >= cut]) - bg < 0.20:
        return None                       # the upper mode is not an element either
    return cut


def support_threshold(ident, bg):
    """A copy supports the consensus if its identity is well above background
    FOR THIS FAMILY. A constant 0.55 flags an old but clean family as
    contaminated - SIMCLEAN at age 0.30 scored 80 with CONTAMINATED - because
    genuine old copies fall below any absolute cut. Scale it instead between the
    measured flank background and the family's own upper quartile.

    The guard matters: in a set with no element q75 is barely above background,
    and a proportional threshold would then pass almost everything. Below a
    minimum contrast there is nothing to support, so nothing is supported."""
    q75 = float(np.percentile(ident, 75))
    contrast = q75 - bg
    if contrast < 0.15:
        return None                      # no element: support is undefined
    return max(bg + 0.10, bg + 0.45 * contrast)
FLANK_SHARE = 0.55          # ungapped flank identity above which two copies share a flank
BG = 0.25                   # unrelated DNA


# Flank-decay readings, computed on 400 bp flanks by flankdecay.py. When a set
# has one, it REPLACES the within-set flank-sharing test: that test compared
# copies over 70 bp and called adjacent similarity "nesting", which produced
# false NESTED_COPIES on hedgehog e2-3 and e2-4. Decay distance distinguishes
# 50 bp of shared context from a satellite that never ends.
try:
    _DECAY = json.load(open("flankdecay.json"))
except Exception:
    _DECAY = {}


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


def subfamily_split(el, pres, ident, min_grp=12, thr=None):
    """Do the supported copies fall into structurally different groups?

    Copy x copy identity inside the element, split on the leading eigenvector.
    A real family varies smoothly with age; two subfamilies with different
    structure separate into blocks whose between-group identity is clearly below
    their within-group identity."""
    sel = np.where(ident >= (SUPPORT if thr is None else thr))[0]
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
    gap = float(within - between)

    # Null: split on a RANDOM vector instead of the leading eigenvector. Any
    # split of a homogeneous family yields a positive gap whose size grows with
    # divergence and shrinks with n - which is why the raw gap tracked age at
    # +0.99 and copy number at -0.97 in families with no subfamily structure at
    # all. Subtracting the null removes both dependencies at their source.
    rng = np.random.default_rng(0)
    null = []
    for _ in range(24):
        r = rng.normal(size=len(sel))
        gr = r > np.median(r)
        if gr.sum() < 4 or (~gr).sum() < 4:
            continue
        w = np.nanmean([S[np.ix_(gr, gr)].mean(), S[np.ix_(~gr, ~gr)].mean()])
        null.append(w - S[np.ix_(gr, ~gr)].mean())
    nullmean = float(np.mean(null)) if null else 0.0
    return {"gap": gap, "gap_excess": float(gap - nullmean), "gap_null": nullmean,
            "sizes": [int(g.sum()), int((~g).sum())],
            "within": float(within), "between": float(between),
            "members": [int(sel[i]) for i in np.where(g)[0]]}


def verdict(path):
    p = parts(path)
    if p is None:
        return None
    n, ident = p["n"], p["ident"]
    dec = _DECAY.get(os.path.basename(path).replace(".aln.fa", ""))

    # background first, because the support threshold is defined against it
    fl_pre = []
    for seqs in (p["lefts"], p["rights"]):
        F0, _ = flank_sharing(seqs)
        vv = F0[np.triu_indices(len(seqs), 1)]
        vv = vv[vv > 0]
        if len(vv):
            fl_pre.append(float(np.median(vv)))
    bg_meas = float(np.mean(fl_pre)) if fl_pre else BG
    # Two separate questions, previously conflated into one threshold:
    #   is there an element at all      -> support_threshold
    #   is this set a MIXTURE           -> contamination_split
    thr = support_threshold(ident, bg_meas)
    cut = contamination_split(ident, bg_meas)
    if thr is None:
        supported = np.zeros(n, bool)     # no element: nothing supports it
    elif cut is not None:
        supported = ident >= cut          # a real mixture: split at the mode boundary
    else:
        supported = np.ones(n, bool)      # unimodal: every copy belongs, however diverged
    n_sup = int(supported.sum())
    v_cut = cut

    FL, sharedL = flank_sharing(p["lefts"])
    FR, sharedR = flank_sharing(p["rights"])
    shared = sharedL | sharedR
    n_shared = int(shared.sum())
    global_elev = float(np.median(np.concatenate([FL[FL > 0], FR[FR > 0]]))) \
        if (FL > 0).any() or (FR > 0).any() else 0.0

    # a "core" copy: supports the consensus AND sits in unique genomic sequence
    core = supported & ~shared
    n_core = int(core.sum())

    # Structure is tested on the SUPPORTED copies only. Contamination masks a
    # scramble - seg_diag against scrambling went from -0.61 on the joint grid
    # to -0.92 once contamination was low - so pruning has to come first.
    sub = subfamily_split(p["el"], p["pres"], ident, thr=thr)
    turb = turbulence(p["el"], p["pres"])

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
    # "Is there an element" is a BINARY question, so the evidence must saturate
    # once it is answered. A cliff of 0.28 means copies match the consensus at
    # 0.57 while unrelated DNA matches at 0.29 - already unambiguous. Scoring on
    # to 0.55 made the statistic report youth, not existence, and cost a clean
    # 30 %-divergent family a third of its element score.
    g_elem = (sat(id_all - flank_bg, 0.10, 0.30) ** 0.5) *              (sat(n_sup / max(1, n), 0.15, 0.98) ** 0.5)
    # Subfamily structure is deliberately NOT scored. Sergei: "subfamily
    # ambiguity doesn't matter, we answer sine/not-sine, not this sine or that
    # sine." Two subfamilies in one set is still a set of SINEs. It stays as a
    # note so the copies can be split later, but it must not cost score.
    # Length spread is softened deliberately. Truncated copies are still copies
    # of that family - Sergei on the truncated sets: "more like a SINE", "looks
    # more like a LINE", never "not an element". Scored on a 0.10-0.50 range it
    # notes heterogeneity without destroying the verdict; at 0.05-0.25 it drove
    # the truncated class to a mean of 9.7.
    g_homog = 1 - sat(cv_core, 0.10, 0.50)
    if dec:
        # BOTH parts of the decay curve matter. Distance alone rates a LINE
        # fragment as isolated, because its shared flank runs only ~50-75 bp
        # before the copies truncate - yet identity right at the edge is 0.89,
        # which says the element does not end there at all.
        g_uniq = (1 - sat(dec["decay_max"], 75, 300)) *                  (1 - sat(dec.get("edge_max") or 0.0, 0.45, 0.85))
    else:
        g_uniq = np.mean([1 - sat(n_shared / max(1, n), 0.05, 0.60),
                          1 - sat(global_elev - BG, 0.05, 0.30)])
    g_ins = sat(tsd_frac, 0.10, 0.60)

    # The element group GATES the score rather than voting in it. As a co-equal
    # geometric term, a set with no element (g_elem 0.31) still scored 56 because
    # its flanks were unique and its lengths uniform - both true and both
    # irrelevant when there is nothing there. Uniqueness and homogeneity only
    # mean something once an element exists.
    rest = (max(g_homog, 1e-6) ** 0.55) * (max(g_uniq, 1e-6) ** 0.45)
    score = 100 * min(1.0, g_elem * rest * (1 + W["insertion"] * g_ins))

    # ---- sub-cases, not one "negative" bucket ---------------------------
    flags = []
    if dec:
        call = dec.get("call")
        if call == "SATELLITE_OR_DUPLICATION":
            flags.append({"code": "NOT_ISOLATED", "n": dec["decay_max"],
                          "text": "Similarity between copies continues %d bp beyond the "
                                  "element and does not reach background — a satellite "
                                  "array or a duplicated region, not independent "
                                  "insertions." % dec["decay_max"]})
        elif call == "ELEMENT_CONTINUES":
            flags.append({"code": "FRAGMENT_OF_LONGER", "n": dec["decay_max"],
                          "text": "The sequence flanking these copies is itself shared "
                                  "(identity %.2f at the edge, reaching background only "
                                  "%d bp out). They look like fragments of a longer "
                                  "element — a LINE 3' end or similar — rather than a "
                                  "short element with its own boundaries."
                                  % (dec["edge_max"], dec["decay_max"])})
        elif call == "ADJACENT_SIMILARITY":
            flags.append({"code": "ADJACENT_SIMILARITY", "n": dec["decay_max"],
                          "text": "Copies share about %d bp of sequence just outside the "
                                  "element, then reach background. Worth noting, but they "
                                  "are in independent locations."
                                  % dec["decay_max"]})
    # Flank sharing is reported whatever the core count. Previously this sat
    # behind "n_core >= 8", so a set whose copies ALL share flanking sequence
    # scored near zero with no explanation at all - hedgehog e1-4 and e2-2 both
    # did exactly that, at flank identity 0.68 against 0.28 for a normal family.
    if (not dec) and global_elev - BG > 0.15 and n >= 15:
        flags.append({"code": "SHARED_FLANKS", "n": n_shared,
                      "text": "Flanking sequence is %.2f identical between copies against "
                              "%.2f for unrelated DNA — these loci are not in independent "
                              "genomic contexts. A satellite array, a duplicated region, or "
                              "copies inside one host repeat." % (global_elev, BG)})
    if g_elem < 0.25:
        flags.append({"code": "NO_ELEMENT",
                      "text": "No element: too few copies support the consensus above background."})
    if n_core >= 8:
        if cut is not None and (n - n_sup) / max(1, n) > 0.03:
            flags.append({"code": "CONTAMINATED", "n": int(n - n_sup),
                          "text": "%d of %d copies do not support the consensus; prune them and "
                                  "re-score (see prune.py)." % (n - n_sup, n)})
        if sub and sub.get("gap_excess", sub["gap"]) > 0.03:
            flags.append({"code": "SUBFAMILY_NOTE", "n": sub["sizes"],
                          "text": "Supported copies split into structurally different groups of "
                                  "%d and %d (within-group identity %.2f, between %.2f, excess "
                                  "over a random split %.3f). Both are SINEs - this does not "
                                  "affect the verdict, but the copies can be split if the "
                                  "subfamilies are wanted separately."
                                  % (sub["sizes"][0], sub["sizes"][1], sub["within"],
                                     sub["between"], sub.get("gap_excess", 0))})
        if (not dec) and n_shared / max(1, n) > 0.15:
            kind = ("satellite or segmental duplication" if global_elev - BG > 0.15
                    else "another repeat or a duplication")
            flags.append({"code": "NESTED_COPIES", "n": n_shared,
                          "text": "%d of %d copies share flanking sequence with each other — they "
                                  "sit inside %s. Set them aside; they do not disqualify the "
                                  "family." % (n_shared, n, kind)})
    if cv_core > 0.18 and n_core >= 15:
        flags.append({"code": "TRUNCATED_COPIES", "n": round(cv_core, 3),
                      "text": "Copy length varies widely (CV %.2f against 0.05 for a clean "
                              "family). Many copies are 5'-truncated. Still the same family, "
                              "but the truncated copies may be worth setting aside."
                              % cv_core})
    if turb is not None and turb > 0.12 and n_core >= 15:
        flags.append({"code": "NEEDS_REVIEW", "n": round(turb, 3),
                      "text": "The element is unsettled along its length: %.0f %% of "
                              "positions drop well below its own median identity, against "
                              "6-8 %% in a clean family. Still a SINE, but one to look at "
                              "by eye before using." % (100 * turb)})
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
            "support_threshold": None if thr is None else round(float(thr), 3),
            "contamination_cut": None if v_cut is None else round(float(v_cut), 3),
            "turbulence": None if turb is None else round(turb, 3),
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
