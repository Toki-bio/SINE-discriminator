#!/usr/bin/env python3
"""Inspect the MEASURE() feature table: per-feature separation and the
subsampling side quest. No model is fitted - the point is to see whether the
statistics behave before anything is trained on them.
"""
import sys, json, math
from collections import defaultdict
import numpy as np

CLASSES = ["POS", "MIXED10", "MIXED30", "NEGJITTER", "NEGTRUNC5", "NEGSPLICE", "NEGRAND"]
SKIP = {"set", "error", "anchorL", "anchorR", "aln_len", "n_copies", "n_edged",
        "n_unique", "core_cols"}


def load(path):
    rows = [json.loads(l) for l in open(path)]
    return rows


def cls_of(s):
    return s.split("__")[0]


def auc(pos, neg):
    """P(random positive > random negative), ties counted as 0.5."""
    pos = np.asarray([x for x in pos if np.isfinite(x)], float)
    neg = np.asarray([x for x in neg if np.isfinite(x)], float)
    if len(pos) < 3 or len(neg) < 3:
        return float("nan")
    allv = np.concatenate([pos, neg])
    r = allv.argsort().argsort().astype(float) + 1
    # average ranks for ties
    order = np.argsort(allv)
    sv = allv[order]
    i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        if j > i:
            r[order[i:j + 1]] = (i + j + 2) / 2.0
        i = j + 1
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def fmt(x, w=7, p=3):
    return ("%*.*f" % (w, p, x)) if np.isfinite(x) else "%*s" % (w, "-")


def main():
    rows = load(sys.argv[1] if len(sys.argv) > 1 else "features.jsonl")
    ok = [r for r in rows if "error" not in r]
    bad = [r for r in rows if "error" in r]

    print("=" * 78)
    print("MEASURE() ran on %d sets: %d complete, %d aborted" % (len(rows), len(ok), len(bad)))
    ec = defaultdict(lambda: defaultdict(int))
    for r in bad:
        ec[cls_of(r["set"])][r["error"]] += 1
    for c in sorted(ec):
        print("  aborted %-12s %s" % (c, dict(ec[c])))

    byc = defaultdict(list)
    for r in ok:
        byc[cls_of(r["set"])].append(r)
    tot = defaultdict(int)
    for r in rows:
        tot[cls_of(r["set"])] += 1

    print("\n" + "=" * 78)
    print("TIER-2 GATE: does a consensus core survive at all?")
    print("An abort is a verdict, not a failure - MEASURE could not find a block")
    print("of >=75 alignment columns where copies are both present and agreeing.")
    print("\n%-12s %6s %8s %8s   %s" % ("class", "sets", "measured", "aborted", "abort rate"))
    print("-" * 60)
    for c in CLASSES + ["SQ"]:
        if not tot[c]:
            continue
        a = tot[c] - len(byc[c])
        bar = "#" * int(round(20.0 * a / tot[c]))
        print("%-12s %6d %8d %8d   %5.0f%% %s"
              % (c, tot[c], len(byc[c]), a, 100.0 * a / tot[c], bar))
    print("\nEverything below is conditioned on surviving this gate, so it is")
    print("measured on the negatives that most resemble real families.")

    feats = sorted({k for r in ok for k in r if k not in SKIP and
                    isinstance(r.get(k), (int, float))})
    pos = byc.get("POS", [])
    negs = [c for c in CLASSES[1:] if c in byc]

    # ---------------- per-feature separation, POS vs EACH negative --------
    print("\n" + "=" * 78)
    print("PER-FEATURE SEPARATION  (AUC: 1.0 = feature is high in real families,")
    print("0.0 = high in that negative, 0.5 = blind.  |AUC-0.5| is the signal.)")
    hdr = "%-18s %8s %8s | " % ("feature", "POSmed", "POSiqr") + " ".join("%9s" % n[:9] for n in negs)
    print("\n" + hdr)
    print("-" * len(hdr))
    table = {}
    for f in feats:
        pv = [r[f] for r in pos if isinstance(r.get(f), (int, float))]
        pv = [x for x in pv if np.isfinite(x)]
        if len(pv) < 3:
            continue
        med = float(np.median(pv))
        iqr = float(np.subtract(*np.percentile(pv, [75, 25])))
        aucs = []
        for nc in negs:
            nv = [r.get(f) for r in byc[nc] if isinstance(r.get(f), (int, float))]
            aucs.append(auc(pv, nv))
        table[f] = aucs
        print("%-18s %s %s | " % (f, fmt(med, 8, 3), fmt(iqr, 8, 3))
              + " ".join(fmt(a, 9, 3) for a in aucs))

    print("\nBest single feature against each negative class:")
    for i, nc in enumerate(negs):
        best = max((f for f in table if np.isfinite(table[f][i])),
                   key=lambda f: abs(table[f][i] - 0.5), default=None)
        if best:
            print("  %-11s %-18s AUC=%.3f" % (nc, best, table[best][i]))
    dead = [f for f in table if all(not np.isfinite(a) or abs(a - .5) < .12 for a in table[f])]
    print("\nFeatures that separate nothing (|AUC-0.5| < 0.12 everywhere): %s"
          % (", ".join(sorted(dead)) or "none"))

    # ---------------- redundancy ------------------------------------------
    print("\n" + "=" * 78)
    print("REDUNDANT PAIRS (|r| > 0.9 across all sets) - keep one of each")
    use = [f for f in table]
    Xall = []
    for f in use:
        Xall.append([r.get(f, np.nan) if isinstance(r.get(f), (int, float)) else np.nan
                     for r in ok])
    X = np.array(Xall, float)
    n_pairs = 0
    for i in range(len(use)):
        for j in range(i + 1, len(use)):
            m = np.isfinite(X[i]) & np.isfinite(X[j])
            if m.sum() < 20 or X[i][m].std() == 0 or X[j][m].std() == 0:
                continue
            r = np.corrcoef(X[i][m], X[j][m])[0, 1]
            if abs(r) > 0.9:
                print("  %-18s %-18s r=%+.2f" % (use[i], use[j], r))
                n_pairs += 1
    if not n_pairs:
        print("  none")

    # ---------------- side quest -------------------------------------------
    sq = [r for r in ok if r["set"].startswith("SQ__")]
    if sq:
        print("\n" + "=" * 78)
        print("SIDE QUEST - subsampling stability (spec section 3)")
        g = defaultdict(list)
        for r in sq:
            p = r["set"].split("__")
            g[(p[1] + "/" + p[2], int(p[3][1:]))].append(r)
        keyf = ["core_cols", "cliff", "cliff_z", "rank1_excess", "len_cv",
                "resL_spread", "resR_spread"]
        print("\n%-22s %4s %4s  %-14s %s" % ("family", "n", "K", "core_cols",
              "  ".join("%-14s" % k for k in keyf[1:])))
        print("-" * 118)
        for fam in sorted({k[0] for k in g}):
            for n in sorted({k[1] for k in g if k[0] == fam}):
                rs = g[(fam, n)]
                cells = []
                for k in keyf:
                    v = np.array([x.get(k, np.nan) for x in rs], float)
                    v = v[np.isfinite(v)]
                    cells.append("%7.2f+-%-5.2f" % (v.mean(), v.std()) if len(v) else "     -      ")
                print("%-22s %4d %4d  %s" % (fam, n, len(rs), "  ".join(cells)))
        print("\nRead: 'core_cols' spread across K independent subsets of the same")
        print("family IS the boundary rule's reproducibility. Where it stops shrinking")
        print("with n is n_min. The +- on the other columns is the sampling-only error")
        print("bar every threshold must clear.")

    # ---------------- headline contrasts -----------------------------------
    print("\n" + "=" * 78)
    print("HEADLINE: medians by class")
    show = ["core_cols", "cliff", "cliff_z", "core_identity", "flank_bg",
            "rank1_frac", "rank1_excess", "rank2_frac", "len_cv",
            "resL_spread", "resR_spread", "res_asymmetry", "tsd_frac",
            "nested_frac_L", "flank_med_L", "arich_score", "tail_score"]
    hdr = "%-16s " % "class" + " ".join("%9s" % s[:9] for s in show)
    print(hdr)
    print("-" * len(hdr))
    for c in CLASSES:
        if c not in byc:
            continue
        line = "%-16s " % ("%s(%d)" % (c, len(byc[c])))
        for s in show:
            v = [r.get(s) for r in byc[c] if isinstance(r.get(s), (int, float))]
            v = [x for x in v if np.isfinite(x)]
            line += fmt(np.median(v) if v else float("nan"), 9, 3) + " "
        print(line)


if __name__ == "__main__":
    main()
