#!/usr/bin/env python3
"""Exploratory edge refinement: treat the consensus edge as a starting guess.

The consensus used for the preliminary search may be mis-trimmed, so the
element boundary is searched over offsets around it rather than taken as given.

Objective. At a candidate edge, the boundary quality is the step in mean
pairwise identity between copies across it:

    step(d) = mean pair_id over [d, d+W]  -  mean pair_id over [d-W, d]

for the left edge (sign reversed on the right). The best edge is the offset
maximising that step.

Circularity. Optimising the same quantity used to discriminate would inflate it
for negatives too, so the refined step value is NOT a discriminating variable
here. What is reported instead is the SHAPE of the landscape:

    d_best      how far the edge should move
    gain        step at the optimum minus step at the consensus edge
    prominence  step at the optimum minus the median step over the scan
    width       offsets within 90 % of the peak - a sharp optimum is narrow

A real family should show a sharp, well-localised optimum. Junk should have a
flat landscape with no preferred edge, whatever its best step happens to be.

Coordinates. Inside the element, positions are consensus columns and copies are
aligned. Outside, positions are ungapped offsets from each copy's own edge, so
no aligner is involved. Both measure the same thing - mean pairwise identity at
a homologous position - but the method changes at the edge; noted, not hidden.
"""
import json, os, glob, sys
import numpy as np
import measure_c as M

W = 25            # half-window either side of a candidate edge
SEARCH = 60       # offsets scanned outward; limited by the flank in the alignment


def pid(col):
    b = col[col != M.GAP]
    n = len(b)
    if n < 6:
        return np.nan
    c = np.bincount(b, minlength=4)[:4].astype(float)
    return float((c * (c - 1)).sum() / (n * (n - 1)))


def edge_profiles(path):
    """Mean pairwise identity on a single axis spanning flank-element-flank."""
    names, A = M.read_aln(path)
    ci = [i for i, n in enumerate(names) if "CONSENSUS_" in n]
    if not ci:
        return None
    k = ci[0]
    cons = A[k]
    nzc = np.where(cons != M.GAP)[0]
    if len(nzc) < 80:
        return None
    lo, hi = int(nzc[0]), int(nzc[-1])
    C = np.delete(A, k, axis=0)
    n = C.shape[0]
    L = len(nzc)

    inside = np.array([pid(C[:, c]) for c in nzc])            # element positions

    lefts, rights = [], []
    for i in range(n):
        l = C[i, :lo]
        lefts.append(l[l != M.GAP][::-1])
        r = C[i, hi + 1:]
        rights.append(r[r != M.GAP])

    def flank_prof(seqs, k_max):
        out = []
        for off in range(1, k_max + 1):
            col = np.array([s[off - 1] if len(s) >= off else M.GAP for s in seqs],
                           dtype=np.int8)
            out.append(pid(col))
        return np.array(out)

    lf = flank_prof(lefts, SEARCH + W + 5)      # index 0 = adjacent to element
    rf = flank_prof(rights, SEARCH + W + 5)
    return {"inside": inside, "lf": lf, "rf": rf, "L": L, "n": n}


def scan_edge(inside, flank, search=SEARCH, w=W):
    """step(d) for d in -search..+search. d<0 extends into the flank, d>0 trims
    into the element. Returns the scan and its summary."""
    # One continuous axis: flank outermost-first, then the element. Indexing the
    # flank and the element separately means a candidate edge inside the element
    # reads its "outside" window off the end of the flank array and averages
    # NaN, which manufactures a spurious optimum.
    axis = np.concatenate([flank[::-1], inside])
    F = len(flank)                        # element position p is axis[F + p]
    ds, steps = [], []
    for d in range(-search, search + 1):
        q = F + d                         # candidate edge in axis coordinates
        if q - w < 0 or q + w > len(axis):
            ds.append(d)
            steps.append(np.nan)
            continue
        o_m = np.nanmean(axis[q - w:q])
        i_m = np.nanmean(axis[q:q + w])
        ds.append(d)
        steps.append(i_m - o_m if np.isfinite(o_m) and np.isfinite(i_m) else np.nan)
    ds, steps = np.array(ds), np.array(steps)
    ok = np.isfinite(steps)
    if ok.sum() < 20:
        return None
    best = int(ds[ok][np.nanargmax(steps[ok])])
    smax = float(np.nanmax(steps))
    at0 = float(steps[ds == 0][0]) if np.isfinite(steps[ds == 0][0]) else np.nan
    med = float(np.nanmedian(steps))
    thresh = med + 0.9 * (smax - med)
    width = int(np.sum(steps[ok] >= thresh))
    return {"d_best": best, "step_best": smax, "step_at_cons": at0,
            "gain": smax - at0 if np.isfinite(at0) else None,
            "prominence": smax - med, "width": width,
            "scan_d": ds.tolist(),
            "scan": [None if not np.isfinite(x) else round(float(x), 4) for x in steps]}


def refine(path):
    p = edge_profiles(path)
    if p is None:
        return None
    left = scan_edge(p["inside"], p["lf"])
    right = scan_edge(p["inside"][::-1], p["rf"])
    out = {"set": os.path.basename(path).replace(".aln.fa", ""),
           "L": p["L"], "n": p["n"]}
    for side, r in (("L", left), ("R", right)):
        if r is None:
            continue
        for key in ("d_best", "step_best", "step_at_cons", "gain", "prominence", "width"):
            out[key + "_" + side] = r[key]
        out["scan_" + side] = r["scan"]
        out["scan_d"] = r["scan_d"]
    if "d_best_L" in out and "d_best_R" in out:
        out["len_change"] = -(out["d_best_L"]) - (out["d_best_R"])
        out["prominence_mean"] = float(np.mean([out["prominence_L"], out["prominence_R"]]))
        out["width_mean"] = float(np.mean([out["width_L"], out["width_R"]]))
    return out


def _one(f):
    try:
        return refine(f)
    except Exception as exc:
        return {"set": os.path.basename(f).replace(".aln.fa", ""), "error": str(exc)}


def main():
    import multiprocessing as mp
    d = sys.argv[1] if len(sys.argv) > 1 else "aln_c"
    out = sys.argv[2] if len(sys.argv) > 2 else "edges.jsonl"
    files = sorted(glob.glob(os.path.join(d, "*.aln.fa")))
    files = [f for f in files if not os.path.basename(f).startswith("SQ__")]
    print("refining %d" % len(files), flush=True)
    with mp.Pool(max(1, mp.cpu_count() - 2)) as pool, open(out, "w") as fh:
        for i, v in enumerate(pool.imap_unordered(_one, files, chunksize=4)):
            if v:
                fh.write(json.dumps(v) + "\n")
            if i % 100 == 0:
                print("  %d" % i, flush=True)
    print("wrote", out)


if __name__ == "__main__":
    main()
