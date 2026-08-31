#!/usr/bin/env python3
"""The refinement loop from spec section 4.3, and its validation.

A mixed set is not a negative. The spec asks for two distinct verdicts - "not a
family at all" versus "a family plus contaminants" - and for the second the
required output is the family, the refined bounds, and a cleaned copy list.
Scoring MIXED as a negative class, as earlier builds did, answers the wrong
question.

Loop: measure per-copy identity to the known consensus, drop the copies that do
not support it, re-measure on what remains, repeat up to three passes or until
nothing more is dropped.

Two cut rules are compared, because the choice is not obvious:
  mad    median - k*MAD, a robust outlier rule that assumes one population
  gap    the largest gap in the sorted identities, i.e. an explicit two-population
         split, searched only in a plausible middle band

Ground truth is available: planted contaminants carry "contam" in their name, so
precision and recall of the removal are measurable rather than assumed. The
control that matters is POS - a clean family should lose almost nothing.
"""
import json, glob, os, sys
import numpy as np
import measure_c as M


def per_copy_identity(A, names):
    ci = [i for i, n in enumerate(names) if "CONSENSUS_" in n]
    if not ci:
        return None
    k = ci[0]
    cons = A[k]
    nz = np.where(cons != M.GAP)[0]
    lo, hi = int(nz[0]), int(nz[-1])
    idx = [i for i in range(len(names)) if i != k]
    C = A[idx][:, lo:hi + 1]
    pres = C != M.GAP
    agree = (C == cons[lo:hi + 1][None, :]) & pres
    with np.errstate(invalid="ignore"):
        ident = agree.sum(axis=1) / np.maximum(pres.sum(axis=1), 1)
    ident[pres.sum(axis=1) < 30] = 0.0
    return ident, [names[i] for i in idx]


def cut_mad(x, k=3.0):
    med = np.median(x)
    mad = np.median(np.abs(x - med)) * 1.4826
    return med - k * max(mad, 0.01)


def cut_gap(x, lo_q=0.05, hi_q=0.95):
    """Largest gap in the sorted values, searched inside a middle band so a
    single stray copy at either extreme cannot define the split."""
    s = np.sort(x)
    a, b = int(len(s) * lo_q), int(len(s) * hi_q)
    if b - a < 3:
        return cut_mad(x)
    d = np.diff(s[a:b])
    j = int(np.argmax(d))
    return float((s[a + j] + s[a + j + 1]) / 2)


def run(path, rule="gap", passes=3):
    names, A = M.read_aln(path)
    r = per_copy_identity(A, names)
    if r is None:
        return None
    ident0, cnames = r
    keep = np.ones(len(ident0), bool)
    hist = []
    for p in range(passes):
        x = ident0[keep]
        if len(x) < 10:
            break
        thr = cut_gap(x) if rule == "gap" else cut_mad(x)
        # never cut above the level that would remove a real family wholesale
        thr = min(thr, 0.6)
        new = keep & (ident0 >= thr)
        hist.append({"pass": p + 1, "thr": round(float(thr), 3),
                     "kept": int(new.sum()), "dropped": int(keep.sum() - new.sum())})
        if new.sum() == keep.sum():
            break
        keep = new
    truth = np.array(["contam" in n for n in cnames])
    removed = ~keep
    tp = int((removed & truth).sum())
    fp = int((removed & ~truth).sum())
    fn = int((~removed & truth).sum())
    return {"set": os.path.basename(path).replace(".aln.fa", ""),
            "n": len(cnames), "n_contam": int(truth.sum()),
            "removed": int(removed.sum()), "tp": tp, "fp": fp, "fn": fn,
            "precision": round(tp / max(1, tp + fp), 3) if removed.sum() else None,
            "recall": round(tp / max(1, tp + fn), 3) if truth.sum() else None,
            "ident_kept_med": round(float(np.median(ident0[keep])), 3) if keep.sum() else None,
            "hist": hist, "keep_mask": keep.tolist()}


def main():
    rule = sys.argv[1] if len(sys.argv) > 1 else "gap"
    files = sorted(glob.glob("aln_v2/*.aln.fa"))
    files = [f for f in files if os.path.basename(f).startswith(("MIXED", "POS", "NEGRAND"))]
    out = {}
    for f in files:
        r = run(f, rule)
        if r:
            out[r["set"]] = r
    json.dump(out, open("prune_%s.json" % rule, "w"), separators=(",", ":"))

    print("REFINEMENT LOOP (%s rule) - contaminant removal against ground truth\n" % rule)
    print("%-12s %5s %8s %10s %10s %10s" % ("class", "sets", "removed", "precision", "recall", "kept ident"))
    print("-" * 62)
    for cls in ("POS", "MIXED10", "MIXED30", "NEGRAND"):
        rs = [v for k, v in out.items() if k.startswith(cls + "__")]
        if not rs:
            continue
        prec = [v["precision"] for v in rs if v["precision"] is not None]
        rec = [v["recall"] for v in rs if v["recall"] is not None]
        print("%-12s %5d %8.1f %10s %10s %10.3f"
              % (cls, len(rs), np.mean([v["removed"] for v in rs]),
                 ("%.3f" % np.mean(prec)) if prec else "-",
                 ("%.3f" % np.mean(rec)) if rec else "-",
                 np.mean([v["ident_kept_med"] for v in rs if v["ident_kept_med"]])))
    print("\nPOS is the control: a clean family should lose almost nothing.")


if __name__ == "__main__":
    main()
