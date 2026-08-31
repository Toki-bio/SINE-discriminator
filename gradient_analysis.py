#!/usr/bin/env python3
"""Which statistic responds to which perturbation, and to nothing else?

For each gradient (one parameter swept over 100 sets) and each statistic,
Spearman rho between the swept parameter and the statistic. That answers two
questions a pass/fail comparison cannot:

  sensitivity  does the statistic move monotonically with the thing it claims
               to measure?
  specificity  does it stay flat when something else is varied? A statistic
               that responds to everything is a general quality score, not a
               diagnostic, and cannot support a named sub-case verdict.

Also reports where each statistic saturates - the parameter value beyond which
it stops changing - because that is the working range, and a threshold set
outside it is meaningless.
"""
import json, glob, os, sys
from collections import defaultdict
import numpy as np
import measure_c as M


def spearman(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 12 or np.std(x[m]) == 0 or np.std(y[m]) == 0:
        return np.nan
    rx = np.argsort(np.argsort(x[m])).astype(float)
    ry = np.argsort(np.argsort(y[m])).astype(float)
    return float(np.corrcoef(rx, ry)[0, 1])


def stats_for(path):
    """Every statistic this project has, on one alignment."""
    import measure_c
    out = {}
    try:
        v = measure_c.measure(path)
        if "error" not in v:
            out.update({k: x for k, x in v.items() if isinstance(x, float)})
    except Exception:
        pass
    try:
        import verdict
        w = verdict.verdict(path)
        if w and "error" not in w:
            out["score"] = w["score"]
            out["n_core_frac"] = w["n_core"] / max(1, w["n"])
            for g, val in w["groups"].items():
                out["g_" + g] = val
            if w.get("subfamily"):
                out["sub_gap"] = w["subfamily"]["gap"]
    except Exception:
        pass
    try:
        import segmap
        sm = segmap.analyse(path)
        if sm:
            out["seg_diag"] = sm["diag_frac"]
    except Exception:
        pass
    # long internal gaps
    try:
        names, A = M.read_aln(path)
        ci = [i for i, n in enumerate(names) if "CONSENSUS_" in n][0]
        nz = np.where(A[ci] != M.GAP)[0]
        C = np.delete(A, ci, axis=0)[:, nz]
        L = C.shape[1]
        keep = (C != M.GAP).sum(axis=1) > 0.35 * L
        C = C[keep]
        if C.shape[0] >= 10:
            wg, starts = 0, []
            for i in range(C.shape[0]):
                g = (C[i] == M.GAP).astype(np.int8)
                d = np.diff(np.concatenate([[0], g, [0]]))
                st, en = np.where(d == 1)[0], np.where(d == -1)[0]
                h = [(int(a), int(b - a)) for a, b in zip(st, en)
                     if b - a >= 20 and a > 5 and b < L - 5]
                if h:
                    wg += 1
                    starts += [a for a, _ in h]
            out["long_gap_frac"] = wg / C.shape[0]
            if len(starts) >= 4:
                hh, _ = np.histogram(starts, bins=max(4, L // 20), range=(0, L))
                out["gap_concentration"] = float(hh.max() / max(1, hh.sum()))
    except Exception:
        pass
    return out


def _one(f):
    return os.path.basename(f).replace(".aln.fa", ""), stats_for(f)


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "grad_aln"
    truth = json.load(open(sys.argv[2] if len(sys.argv) > 2 else "grad_truth.json"))
    files = sorted(glob.glob(os.path.join(d, "GRAD_*.aln.fa")))
    print("measuring %d gradient alignments" % len(files), flush=True)
    import multiprocessing as mp
    res = {}
    with mp.Pool(max(1, mp.cpu_count() - 2)) as pool:
        for i, (sid, st) in enumerate(pool.imap_unordered(_one, files, chunksize=8)):
            res[sid] = st
            if i % 200 == 0:
                print("  %d" % i, flush=True)
    json.dump(res, open("grad_stats.json", "w"), separators=(",", ":"))

    grids = defaultdict(list)
    for sid, t in truth.items():
        if sid in res:
            grids[t["grid"]].append((t["value"], res[sid]))
    keys = sorted({k for v in res.values() for k in v})

    table = {}
    for g, rows in grids.items():
        x = np.array([r[0] for r in rows], float)
        table[g] = {}
        for k in keys:
            y = np.array([r[1].get(k, np.nan) for r in rows], float)
            table[g][k] = spearman(x, y)
    json.dump(table, open("grad_response.json", "w"), indent=1)

    gnames = [g for g, _, _, _ in
              [("AGE",)*4, ("NCOPIES",)*4] if False] or sorted(grids)
    interesting = [k for k in keys
                   if max(abs(table[g].get(k, 0) or 0) for g in gnames) > 0.5]
    print("\nDOSE-RESPONSE  |rho| Spearman, statistic vs the swept parameter")
    print("bold |rho|>0.8 = strong.  A row responding to many columns is a")
    print("general quality score, not a diagnostic.\n")
    hdr = "%-20s" % "statistic" + "".join("%9s" % g[:8] for g in gnames)
    print(hdr)
    print("-" * len(hdr))
    for k in sorted(interesting,
                    key=lambda k: -max(abs(table[g].get(k, 0) or 0) for g in gnames)):
        row = "%-20s" % k[:20]
        for g in gnames:
            r = table[g].get(k)
            row += "%9s" % ("  -  " if r is None or not np.isfinite(r)
                            else ("%+.2f" % r))
        n_resp = sum(1 for g in gnames
                     if table[g].get(k) is not None and np.isfinite(table[g][k])
                     and abs(table[g][k]) > 0.6)
        print(row + "   <- responds to %d" % n_resp)


if __name__ == "__main__":
    main()
