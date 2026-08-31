#!/usr/bin/env python3
"""Positional profiles: every statistic that varies along the element, plotted
against nucleotide position rather than collapsed to one number.

The x axis has three zones stitched together, so the Tier-2 cliff is visible as
a curve rather than inferred from a scalar:

    -70 .. -1     left flank, offset outward from the element edge
      0 .. L-1    consensus positions of the element itself
     +1 .. +70    right flank, offset outward from the element edge

Flank positions are indexed from each copy's OWN edge on ungapped sequence, so
no aligner is involved out there and the background is the true one.

Tracks
  pair_id   mean pairwise identity between copies - the one metric defined
            identically in all three zones, so the cliff is directly readable
  cover     fraction of copies present - shows truncation as a ramp
  cons_id   identity to the known consensus (element only)
  mosaic    per-window residual after a rank-1 fit of the copy x window
            identity matrix. Rank 1 means "every copy is just older or younger";
            what does not fit that is positional structure specific to a subset
            of copies, i.e. the spec's mosaicism, localised.
  at        A+T fraction - the A-rich insertion signature and the poly-A tail
"""
import json, sys, os, glob
import numpy as np
import measure_c as M

FL = 70          # flank offsets profiled either side
WIN, STEP = 20, 4


def pair_identity_col(col):
    """Mean pairwise identity within one column of symbols (gaps excluded)."""
    b = col[col != M.GAP]
    n = len(b)
    if n < 4:
        return np.nan, n
    cnt = np.bincount(b, minlength=4)[:4].astype(float)
    same = float((cnt * (cnt - 1)).sum())
    return same / (n * (n - 1)), n


def profile(path):
    names, A = M.read_aln(path)
    ci = [i for i, n in enumerate(names) if "CONSENSUS_" in n]
    if not ci:
        return None
    k = ci[0]
    cons = A[k]
    nzc = np.where(cons != M.GAP)[0]           # consensus position -> aln column
    lo, hi = int(nzc[0]), int(nzc[-1])
    C = np.delete(A, k, axis=0)
    n = C.shape[0]
    L = len(nzc)

    xs, pair, cover, consid, at = [], [], [], [], []

    # ---- left flank, ungapped, indexed outward from each copy's own edge ----
    lefts = []
    for i in range(n):
        r = C[i, :lo]
        r = r[r != M.GAP]
        lefts.append(r[::-1])
    for off in range(FL, 0, -1):
        col = np.array([l[off - 1] if len(l) >= off else M.GAP for l in lefts], dtype=np.int8)
        p, m = pair_identity_col(col)
        xs.append(-off)
        pair.append(p)
        cover.append(m / n)
        consid.append(np.nan)
        b = col[col != M.GAP]
        at.append(float(np.mean(np.isin(b, [0, 3]))) if len(b) > 3 else np.nan)

    # ---- element, one point per consensus position -------------------------
    for p_i, c in enumerate(nzc):
        col = C[:, c]
        p, m = pair_identity_col(col)
        xs.append(p_i)
        pair.append(p)
        cover.append(m / n)
        present = col != M.GAP
        consid.append(float(np.mean(col[present] == cons[c])) if present.sum() >= 4 else np.nan)
        b = col[present]
        at.append(float(np.mean(np.isin(b, [0, 3]))) if len(b) > 3 else np.nan)

    # ---- right flank -------------------------------------------------------
    rights = []
    for i in range(n):
        r = C[i, hi + 1:]
        rights.append(r[r != M.GAP])
    for off in range(1, FL + 1):
        col = np.array([r[off - 1] if len(r) >= off else M.GAP for r in rights], dtype=np.int8)
        p, m = pair_identity_col(col)
        xs.append(L - 1 + off)
        pair.append(p)
        cover.append(m / n)
        consid.append(np.nan)
        b = col[col != M.GAP]
        at.append(float(np.mean(np.isin(b, [0, 3]))) if len(b) > 3 else np.nan)

    # ---- mosaicism: rank-1 residual per window over consensus positions -----
    El = C[:, nzc]
    pres = El != M.GAP
    agree = (El == cons[nzc][None, :]) & pres
    starts = np.arange(0, max(1, L - WIN + 1), STEP)
    Mw = np.full((n, len(starts)), np.nan)
    for j, s in enumerate(starts):
        e = min(s + WIN, L)
        den = pres[:, s:e].sum(axis=1)
        with np.errstate(invalid="ignore"):
            Mw[:, j] = np.where(den >= 6, agree[:, s:e].sum(axis=1) / np.maximum(den, 1), np.nan)
    keep = np.isfinite(Mw).mean(axis=1) > 0.8
    Mk = Mw[keep]
    colmed = np.nanmedian(Mk, axis=0)
    Mk = np.where(np.isfinite(Mk), Mk, colmed[None, :])
    mosaic = np.full(len(starts), np.nan)
    if Mk.shape[0] >= 8:
        U, S, Vt = np.linalg.svd(Mk - Mk.mean(), full_matrices=False)
        fit = np.outer(U[:, 0] * S[0], Vt[0])
        resid = (Mk - Mk.mean()) - fit
        denom = np.linalg.norm(Mk - Mk.mean(), axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            mosaic = np.where(denom > 1e-9, np.linalg.norm(resid, axis=0) / denom, np.nan)
    mos_x = (starts + WIN // 2).tolist()

    # ---- motif profiles: score the motif at EVERY position, not just its best
    # hit. A solid box says only where the best match is; the profile shows how
    # well every other part of the element matches too, so a duplicated A box, a
    # decayed second copy of the tRNA head, or a box that is barely better than
    # background all become visible.
    import re as _re
    IU = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]",
          "W": "[AT]", "S": "[GC]", "K": "[GT]", "M": "[AC]", "N": "[ACGT]"}
    cseq = "".join("ACGT"[c] for c in cons[nzc])

    def motif_track(motif):
        pat = [IU[c] for c in motif]
        k = len(motif)
        out = []
        for i in range(L):
            if i + k > L:
                out.append(None)
                continue
            hit = sum(1 for j in range(k) if _re.match(pat[j], cseq[i + j]))
            out.append(round(hit / k, 3))
        return out

    abox_p = motif_track("TRGCNNARYGG")
    bbox_p = motif_track("GWTCRANNC")

    # self-similarity to the tRNA-derived head, to expose internal duplications
    head_lo, head_hi = 5, min(80, L)
    head = cseq[head_lo:head_hi]
    hl = len(head)
    self_p = []
    for i in range(L):
        if i + hl > L:
            self_p.append(None)
            continue
        seg = cseq[i:i + hl]
        self_p.append(round(sum(1 for a, b in zip(seg, head) if a == b) / hl, 3))

    f = lambda a: [None if (x is None or not np.isfinite(x)) else round(float(x), 4) for x in a]
    return {"x": xs, "elem_len": L, "n": int(n),
            "pair_id": f(pair), "cover": f(cover), "cons_id": f(consid), "at": f(at),
            "mosaic_x": mos_x, "mosaic": f(mosaic),
            "motif_x": list(range(L)),
            "abox_p": abox_p, "bbox_p": bbox_p, "selfsim_p": self_p}


def main():
    sets = [os.path.basename(p)[:-7] for p in sorted(glob.glob("aln_c/*.aln.fa"))]
    want = sys.argv[1:] if len(sys.argv) > 1 else None
    if want:
        sets = [s for s in sets if s in want]
    out = {}
    for s in sets:
        p = profile(os.path.join("aln_c", s + ".aln.fa"))
        if p:
            out[s] = p
    json.dump(out, open("profiles.json", "w"), separators=(",", ":"))
    print("profiled %d sets, %.0f KB" % (len(out), os.path.getsize("profiles.json") / 1024))


if __name__ == "__main__":
    main()
