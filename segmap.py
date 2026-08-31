#!/usr/bin/env python3
"""Segment mapping: which part of the consensus does each slot of each copy
actually resemble?

Sergei's definition of mosaicism, which is not the two-parent recombination the
earlier code assumed:

    consensus  a b c d
    copy 1     a b b d      slot 3 carries segment b
    copy 2     a d c d      slot 2 carries segment d

The element's own segments are reshuffled or duplicated between slots, and which
slot is affected varies from copy to copy. Neither the rank-1 residual nor
partition congruence can see this: every copy is still built entirely of this
element's own parts, in nearly the right proportions.

Method - a per-copy dotplot against the consensus, reduced to one number:
for each copy and each slot j, compare that slot's sequence to EVERY consensus
segment k and take the best match. A clean copy maps j -> j for every slot. A
scrambled copy maps some j -> k, k != j.

  diag_frac      fraction of (copy, slot) pairs mapping to their own segment
  offdiag_frac   fraction mapping elsewhere, with a clear margin
  slot_entropy   per slot, the entropy of which segment the copies match - a
                 slot that is scrambled in only some copies has high entropy
"""
import json, glob, os, sys
from collections import defaultdict, Counter
import numpy as np
import measure_c as M

NSEG = 8            # consensus segments; ~30 bp each for a 250 bp element
MARGIN = 0.06       # a call must beat the runner-up by this much to count


def analyse(path, nseg=NSEG):
    names, A = M.read_aln(path)
    ci = [i for i, n in enumerate(names) if "CONSENSUS_" in n]
    if not ci:
        return None
    k = ci[0]
    cons = A[k]
    nz = np.where(cons != M.GAP)[0]
    L = len(nz)
    if L < nseg * 20:
        return None
    C = np.delete(A, k, axis=0)[:, nz]          # copies in consensus coordinates
    cseq = cons[nz]
    full = (C != M.GAP).sum(axis=1) > 0.75 * L
    C = C[full]
    if C.shape[0] < 15:
        return None

    b = [round(i * L / nseg) for i in range(nseg + 1)]
    segs = [cseq[b[i]:b[i + 1]] for i in range(nseg)]

    # UNALIGNED slots. On aligned columns this statistic is blind: MAFFT places a
    # duplicated or swapped segment back at the segment it came from, so a slot
    # substitution becomes an indel pair. Measured across a full scramble
    # gradient, the aligned version moved by 0.019; the unaligned version by 0.40.
    raw = []
    for i in range(C.shape[0]):
        r = C[i]
        raw.append(r[r != M.GAP])
    Cn = np.full((C.shape[0], L), M.GAP, dtype=np.int8)
    for i, r in enumerate(raw):
        rb = [round(j * len(r) / nseg) for j in range(nseg + 1)]
        for j in range(nseg):
            piece = r[rb[j]:rb[j + 1]]
            w = min(len(piece), b[j + 1] - b[j])
            Cn[i, b[j]:b[j] + w] = piece[:w]
    C = Cn

    best = np.full((C.shape[0], nseg), -1, np.int8)
    for j in range(nseg):
        slot = C[:, b[j]:b[j + 1]]
        sc = np.full((C.shape[0], nseg), np.nan)
        for kk in range(nseg):
            ref = segs[kk]
            m = min(slot.shape[1], len(ref))
            if m < 12:
                continue
            sub = slot[:, :m]
            ok = sub != M.GAP
            agree = (sub == ref[None, :m]) & ok
            with np.errstate(invalid="ignore"):
                sc[:, kk] = agree.sum(axis=1) / np.maximum(ok.sum(axis=1), 1)
        for i in range(C.shape[0]):
            row = sc[i]
            if not np.isfinite(row).any():
                continue
            order = np.argsort(np.where(np.isfinite(row), row, -1))[::-1]
            if row[order[0]] - row[order[1]] >= MARGIN:
                best[i, j] = order[0]

    called = best >= 0
    if called.sum() < nseg * 5:
        return None
    diag = (best == np.arange(nseg)[None, :]) & called
    diag_frac = float(diag.sum() / called.sum())

    ent, mism = [], []
    for j in range(nseg):
        col = best[:, j][called[:, j]]
        if len(col) < 8:
            ent.append(None)
            mism.append(None)
            continue
        cnt = Counter(col.tolist())
        p = np.array(list(cnt.values()), float)
        p /= p.sum()
        ent.append(round(float(-(p * np.log2(p)).sum()), 3))
        mism.append(round(float(np.mean(col != j)), 3))

    # which off-diagonal moves actually occur, for the report
    moves = Counter()
    for i in range(C.shape[0]):
        for j in range(nseg):
            if called[i, j] and best[i, j] != j:
                moves[(j, int(best[i, j]))] += 1
    top = [{"slot": a, "carries": c, "n": n} for (a, c), n in moves.most_common(5)]

    return {"set": os.path.basename(path).replace(".aln.fa", ""),
            "n_used": int(C.shape[0]), "nseg": nseg,
            "diag_frac": round(diag_frac, 4),
            "offdiag_frac": round(1 - diag_frac, 4),
            "called_frac": round(float(called.mean()), 3),
            "slot_entropy": ent, "slot_mismatch": mism,
            "max_slot_mismatch": round(float(max([m for m in mism if m is not None], default=0)), 3),
            "top_moves": top}


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
    json.dump(res, open("segmap.json", "w"), separators=(",", ":"))

    byc = defaultdict(list)
    for s, v in res.items():
        byc[s.split("__")[0]].append(v)
    print("SEGMENT MAPPING - does each slot of each copy match its OWN consensus segment?\n")
    print("%-14s %5s %11s %14s %12s" % ("class", "sets", "diag_frac", "worst slot", "called"))
    print("-" * 62)
    order = ["SIMCLEAN", "SIMSCRAM", "POS", "SIMSUBFAM", "SIMMOSAIC", "MIXSUBFAM",
             "NEGMOSAIC", "NEGCHIM", "NEGRAND"]
    for c in order:
        v = byc.get(c, [])
        if not v:
            continue
        print("%-14s %5d %11.3f %14.3f %12.2f"
              % (c, len(v), np.mean([x["diag_frac"] for x in v]),
                 np.mean([x["max_slot_mismatch"] for x in v]),
                 np.mean([x["called_frac"] for x in v])))
    print("\nper-set, the simulated grid (known truth):")
    for s in sorted(res):
        if not s.startswith(("SIMSCRAM", "SIMCLEAN")):
            continue
        v = res[s]
        mv = ", ".join("slot%d<-seg%d x%d" % (m["slot"], m["carries"], m["n"])
                       for m in v["top_moves"][:3])
        print("  %-26s diag %.3f  worst slot %.2f   %s"
              % (s, v["diag_frac"], v["max_slot_mismatch"], mv or "-"))


if __name__ == "__main__":
    main()
