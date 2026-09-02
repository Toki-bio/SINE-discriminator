#!/usr/bin/env python3
"""Flank decay: how FAR outward the similarity between copies persists.

Sergei kept saying the same thing about different sets - "extend it and the
picture may change" - and he was right. A single averaged flank identity blurs
two completely different situations together:

  a real insertion   similarity is at background from the first base outside
  everything else    similarity continues outward, and HOW FAR separates the
                     classes

Measured at 400 bp (which costs nothing now that flanks are justified, never
aligned):

  set                 at the edge   50 bp out   300 bp out
  real SINE               0.28        0.25         0.27
  hedgehog e2-3           0.71        0.25         0.26     <- 50 bp of adjacent
  LINE fragment           0.92        0.31         0.26        similarity only
  satellite               0.87        0.73         0.48
  segmental duplication   0.90        0.88         0.92

Two numbers capture it: identity right at the edge, and the distance at which it
reaches background. A SINE ends; a satellite and a duplication do not.
"""
import json, glob, os, sys
import numpy as np
import measure_c as M

BG = 0.25
STEP, MAXOFF = 25, 400


def pid(col):
    b = col[col != M.GAP]
    if len(b) < 8:
        return np.nan
    c = np.bincount(b, minlength=4)[:4].astype(float)
    return float((c * (c - 1)).sum() / (len(b) * (len(b) - 1)))


def decay(path):
    names, A = M.read_aln(path)
    ci = [i for i, n in enumerate(names) if "CONSENSUS_" in n]
    if not ci:
        return None
    k = ci[0]
    nz = np.where(A[k] != M.GAP)[0]
    lo, hi = int(nz[0]), int(nz[-1])
    C = np.delete(A, k, axis=0)
    out = {}
    for side, block, rev in (("L", C[:, :lo], True), ("R", C[:, hi + 1:], False)):
        seqs = []
        for r in block:
            x = r[r != M.GAP]
            seqs.append(x[::-1] if rev else x)
        prof, offs = [], list(range(0, MAXOFF, STEP))
        for off in offs:
            col = np.array([s[off] if len(s) > off else M.GAP for s in seqs], dtype=np.int8)
            prof.append(pid(col))
        prof = np.array(prof, float)
        edge = float(np.nanmean(prof[:2])) if np.isfinite(prof[:2]).any() else np.nan
        # first offset at which identity has reached background
        dist = None
        for i, v in enumerate(prof):
            if np.isfinite(v) and v < BG + 0.06:
                dist = offs[i]
                break
        if dist is None:
            # Never saw a drop to background - but "never saw" is not "never
            # happened". These copies usually do not carry the full 400 bp, so
            # the profile goes unmeasurable after the first two or three
            # offsets, and defaulting to MAXOFF claimed similarity persisted
            # over 400 bp that nobody looked at. That collapsed the uniqueness
            # term to zero and rejected POS__saq__s4_60seqs, whose right-flank
            # profile is measurable only to 50 bp.
            #
            # Claim only as far as the profile was actually readable.
            last_ok = 0
            for i, v in enumerate(prof):
                if np.isfinite(v):
                    last_ok = offs[i]
            dist = last_ok if (np.isfinite(edge) and edge > BG + 0.06) else 0
        out["edge_" + side] = round(edge, 3) if np.isfinite(edge) else None
        out["decay_" + side] = dist
        out["prof_" + side] = [None if not np.isfinite(x) else round(float(x), 3) for x in prof]
    e = [out.get("edge_L"), out.get("edge_R")]
    e = [x for x in e if x is not None]
    out["edge_max"] = round(max(e), 3) if e else None
    out["decay_max"] = max(out.get("decay_L", 0), out.get("decay_R", 0))
    return out


def classify(d):
    """The reading the profiles support, stated as a rule so it can be tested."""
    if d.get("edge_max") is None:
        return "unknown"
    # The distance alone is not enough, and reading it alone was a real bug: a
    # set whose flank identity is 0.26 at the edge - pure background - was called
    # a satellite because the profile never *confidently* dropped below the
    # threshold, so `dist` fell through to its 400 default. That default is an
    # unmeasured value being scored as guilt, and it rejected 10 sets outright,
    # POS__saq__s4_60seqs among them.
    #
    # A satellite or a duplication is similar THROUGHOUT: the real ones here run
    # 0.73-0.95 at the edge, while every one of those 10 sat between 0.26 and
    # 0.55. Distance means nothing until there is elevated similarity to decay
    # from.
    if d["decay_max"] >= 300 and d["edge_max"] > 0.60:
        return "SATELLITE_OR_DUPLICATION"     # similarity never ends
    if d["edge_max"] > 0.55 and d["decay_max"] >= 50:
        return "ELEMENT_CONTINUES"            # a fragment of something longer
    if d["edge_max"] > 0.55:
        return "ADJACENT_SIMILARITY"          # only the first ~50 bp
    return "ISOLATED_INSERTION"               # ends at a real boundary


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "aln_ext"
    res = {}
    for f in sorted(glob.glob(os.path.join(d, "*.aln.fa"))):
        s = os.path.basename(f)[:-7]
        try:
            r = decay(f)
        except Exception:
            r = None
        if r:
            r["call"] = classify(r)
            res[s] = r
    json.dump(res, open("flankdecay.json", "w"), separators=(",", ":"))
    print("%-26s %9s %9s %8s %8s   %s"
          % ("set", "edge 5'", "edge 3'", "decay5'", "decay3'", "reading"))
    print("-" * 92)
    for s in sorted(res):
        r = res[s]
        print("%-26s %9s %9s %8s %8s   %s"
              % (s[:26], r.get("edge_L"), r.get("edge_R"),
                 r.get("decay_L"), r.get("decay_R"), r["call"]))


if __name__ == "__main__":
    main()
