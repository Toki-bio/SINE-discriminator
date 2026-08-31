#!/usr/bin/env python3
"""Pick the sets worth looking at by eye and pack them for the report page.

Selection is deliberate: clean examples of each class, plus the borderline
cases - the real family that scores worst and the negatives that score most
like real families. Those are where a human eye is actually needed.
"""
import json, os, glob
import numpy as np
import measure_c as M

ROWS = 36          # copies shown per alignment
PAD = 70           # flank columns kept around the element


def load(path):
    names, A = M.read_aln(path)
    ci = [i for i, n in enumerate(names) if "CONSENSUS_" in n][0]
    nz = np.where(A[ci] != M.GAP)[0]
    lo, hi = int(nz[0]), int(nz[-1])
    keep = [ci] + [i for i in range(len(names)) if i != ci][:ROWS]
    a, b = max(0, lo - PAD), min(A.shape[1], hi + PAD + 1)
    sub = A[keep, a:b]
    txt = ["".join("ACGT-"[c] for c in row) for row in sub]
    # keep the genomic coordinates - they are what makes a manual check possible
    nm = []
    for i in keep:
        h = names[i].replace("_R_", "")
        nm.append(h.split("::")[-1] if "::" in h else h)
    bits = [names[i].split("|")[1].split("::")[0] for i in keep[1:]
            if "|" in names[i] and names[i].split("|")[1].split("::")[0].isdigit()]
    return {"names": nm, "rows": txt, "lo": lo - a, "hi": hi - a,
            "n_total": len(names) - 1,
            "bitscore_med": int(np.median([int(b) for b in bits])) if bits else None}


def main():
    feats = {}
    for line in open("features_c.jsonl"):
        r = json.loads(line)
        feats[r["set"]] = r
    ok = [r for r in feats.values() if "error" not in r]
    bycls = {}
    for r in ok:
        bycls.setdefault(r["set"].split("__")[0], []).append(r)

    picks = []

    def add(r, why, tag):
        picks.append({"set": r["set"], "why": why, "tag": tag})

    # clean, typical members of each class - one per species where possible
    for sp in ("saq", "ccr", "teu", "dmo"):
        c = [r for r in bycls["POS"] if r["set"].split("__")[1] == sp]
        if c:
            best = max(c, key=lambda r: r.get("cliff", 0))
            add(best, "typical real family, %s" % sp, "good")
    for cls, why in (("NEGRAND", "random genomic loci - the pure null"),
                     ("NEGTRUNC5", "5'-truncated, LINE-fragment mimic"),
                     ("MIXED30", "30 % random-locus contamination"),
                     ("MIXED10", "10 % random-locus contamination"),
                     ("NEGSPLICE", "chimera of two families"),
                     ("NEGJITTER", "edges displaced 20-100 bp")):
        c = bycls.get(cls, [])
        if c:
            add(sorted(c, key=lambda r: r.get("cliff", 0))[len(c) // 2], why, "bad")

    # the borderline cases - where the numbers and the eye may disagree
    pos = bycls["POS"]
    add(min(pos, key=lambda r: r.get("cliff", 9)),
        "REAL family with the weakest cliff - does it still look real?", "edge")
    add(max(pos, key=lambda r: r.get("elem_len_cv", 0)),
        "REAL family with the most variable copy length", "edge")
    add(min(pos, key=lambda r: r.get("tsd_frac", 9)),
        "REAL family with the fewest detected TSDs", "edge")
    for cls in ("MIXED10", "NEGTRUNC5", "NEGSPLICE"):
        c = bycls.get(cls, [])
        if c:
            add(max(c, key=lambda r: r.get("cliff", 0)),
                "%s that most resembles a real family" % cls, "edge")
    # subsampling stability: same family, smallest and largest n
    sq = bycls.get("SQ", [])
    for n in (25, 200):
        c = [r for r in sq if "__n%d__" % n in r["set"] and "s1_30seqs" in r["set"]]
        if c:
            add(c[0], "side quest: saq/s1 at n=%d - same family, resampled" % n, "sq")

    out = []
    seen = set()
    for p in picks:
        if p["set"] in seen:
            continue
        seen.add(p["set"])
        f = os.path.join("aln_c", p["set"] + ".aln.fa")
        if not os.path.exists(f):
            continue
        d = load(f)
        d.update(p)
        d["feat"] = {k: v for k, v in feats[p["set"]].items()
                     if isinstance(v, (int, float))}
        out.append(d)
    json.dump(out, open("viewer_data.json", "w"))
    print("packed %d alignments, %.0f KB"
          % (len(out), os.path.getsize("viewer_data.json") / 1024))


if __name__ == "__main__":
    main()
