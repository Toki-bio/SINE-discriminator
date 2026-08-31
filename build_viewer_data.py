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
    for line in open("features_all.jsonl"):
        r = json.loads(line)
        feats[r["set"]] = r
    V = json.load(open("verdicts.json"))

    # Deliberately weighted toward negatives and the grey zone, which is where
    # judgement is actually needed - a page of clean positives proves nothing.
    PICKS = [
        ("POS__saq__s5_5seqs", "clear SINE, a curated family", "good"),
        ("POS__teu__t2_75seqs", "clear SINE, a second genome", "good"),
        ("SIMCLEAN__age005", "clear SINE, simulated, young", "good"),
        ("SIMCLEAN__age030", "simulated CLEAN family at high age - scored down, and "
                             "flagged contaminated, purely because it is old", "edge"),
        ("POS__ccr__a_ccr", "a curated family that lands in the grey zone", "edge"),
        ("POS__ccr__g7_3seqs", "curated, grey zone, flagged contaminated", "edge"),
        ("NEGTRUNC5__saq__s5_5seqs", "grey zone: 5'-truncated copies", "edge"),
        ("NEGTRUNC5__teu__t1_45seqs", "grey zone: truncation plus other flags", "edge"),
        ("MIXED30__saq__s3_43seqs", "a real family under 30 % contamination - "
                                    "recoverable, and the pruning loop recovers it", "edge"),
        ("MIXED10__ccr__g6_58seqs", "10 % contamination, grey zone", "edge"),
        ("NEGRAND__dmo__r00", "clear non-SINE: random genomic loci", "bad"),
        ("NEGRAND__saq__r00", "clear non-SINE, second genome", "bad"),
        ("NEGCHIM__saq__s3_43seqs", "clear non-SINE: half element, half foreign DNA", "bad"),
        ("NEGCHIM__ccr__g5_7seqs", "grey zone chimera - scores far higher than it should", "edge"),
        ("MIXSUBFAM__saq__s1_30seqs_s8_225seqs", "two subfamilies mixed - a SINE, but "
                                                 "one that should be split", "edge"),
        ("NEGMOSAIC__saq__s1_30seqs_s8_225seqs", "per-copy recombinants - scores as a "
                                                 "clean family, undetected", "edge"),
        ("SIMSCRAM__swap__f100", "segments swapped between slots - the mosaic Sergei means", "edge"),
        ("SIMDEL__one__f050", "half the copies missing the same 50 bp block", "edge"),
        ("SIMNEST__f050", "half the copies sharing a host flank", "edge"),
        ("SIMTRUNC__r060", "simulated heavy truncation", "bad"),
        ("NEGLINE__teu__r00", "REAL LINE 3' ends — a natural non-SINE that scores 99 "
                              "and is not caught by anything here", "bad"),
        ("NEGLINE__teu__r03", "LINE 3' ends, second sample", "bad"),
        ("ERI__eri__e1-1", "a second SINE family: hedgehog, clean", "good"),
        ("ERI__eri__e1-4", "hedgehog subfamily whose copies all share flanking "
                           "sequence — not independent insertions", "bad"),
        ("ERI__eri__e2-3", "hedgehog, partly nested", "edge"),
    ]
    out, seen = [], set()
    for st, why, tag in PICKS:
        if st in seen:
            continue
        f = os.path.join("aln_v2", st + ".aln.fa")
        if not os.path.exists(f) or st not in feats:
            print("  missing:", st)
            continue
        seen.add(st)
        d = load(f)
        d.update({"set": st, "why": why, "tag": tag})
        d["feat"] = {k: v for k, v in feats[st].items() if isinstance(v, (int, float))}
        out.append(d)
    json.dump(out, open("viewer_data.json", "w"))
    print("packed %d alignments, %.0f KB"
          % (len(out), os.path.getsize("viewer_data.json") / 1024))


if __name__ == "__main__":
    main()
