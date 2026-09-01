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


# Sets from the human and Timema benchmarks, so those legs get profile graphs
# too. Chosen for what they demonstrate, not for looking good.
EXTRA_PICKS = [
    ("TIM__top100__t1_1", "Timema, manually curated by Sergei and boundary-confirmed - "
                          "what a known-good answer looks like", "good", "tim_bench"),
    ("TIM__top100__t2", "a second curated Timema subfamily", "good", "tim_bench"),
    ("TIM__subfam__t3-1", "curated, but element-only: no genomic flanks at all. Scored 0.0 "
                          "until the flank background stopped being computed from a few "
                          "stub bases", "edge", "tim_bench"),
    ("TIMB__timb__SINE_17", "AnnoSINE candidate Sergei judged real; scores 100 clean", "good", "bench"),
    ("TIMB__timb__SINE_47", "judged real, scored 44 at 50 bp and 94 at 400 bp - the flank "
                            "width was the whole disagreement", "edge", "bench"),
    ("TIMB__timb__SINE_25", "judged real; at 400 bp the element continues past its "
                            "annotated boundary", "edge", "bench"),
    ("TIMB__timb__SINE_43", "labelled matched, but its copies sit at identity 0.38 against "
                            "a 0.25 baseline - consensus resembles a family, loci do not "
                            "support it", "bad", "bench"),
    ("HUM__hum__AluY", "human, a young Alu - the easiest possible positive", "good", "bench"),
    ("HUM__hum__MIR1_Amn", "human, ancient: consensus longer than the element it describes. "
                           "Scored 0.0 NO_ELEMENT until CONSENSUS_OVEREXTENDED", "edge", "bench"),
    ("HUM__hum__AluYk11", "human, a genuine mixture - verdict deferred rather than scored", "edge", "bench"),
    ("HUM__hum__LmeSINE1c", "human, open false positive: 69.2 here, 7.1 % RepeatMasker "
                            "overlap, and no internal statistic catches it", "bad", "bench"),
]


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
        ("SIMCLEAN__age030", "simulated CLEAN family at high age - Sergei: pure clear SINE with high divergence", "good"),
        ("POS__ccr__a_ccr", "a curated family that lands in the grey zone", "edge"),
        ("POS__ccr__g7_3seqs", "curated, grey zone, flagged contaminated", "edge"),
        ("NEGTRUNC5__saq__s5_5seqs", "grey zone: 5'-truncated copies", "edge"),
        ("NEGTRUNC5__teu__t1_45seqs", "grey zone: truncation plus other flags", "edge"),
        ("MIXED30__saq__s3_43seqs", "a real family under 30 % contamination - "
                                    "recoverable, and the pruning loop recovers it", "edge"),
        ("MIXED10__ccr__g6_58seqs", "10 % contamination, grey zone", "edge"),
        ("NEGRAND__dmo__r00", "clear non-SINE: random genomic loci", "bad"),
        ("NEGRAND__saq__r00", "clear non-SINE, second genome", "bad"),
        ("NEGCHIM__saq__s3_43seqs", "built as half element, half foreign DNA - but Sergei read it as a clear SINE with a wobbly right end, and the tool now agrees", "edge"),
        ("NEGCHIM__ccr__g5_7seqs", "chimera by construction; scores high because the element in it is real. CONSENSUS_OVEREXTENDED now names the join", "edge"),
        ("MIXSUBFAM__saq__s1_30seqs_s8_225seqs", "two subfamilies mixed - a SINE, but "
                                                 "one that should be split", "edge"),
        ("NEGMOSAIC__saq__s1_30seqs_s8_225seqs", "per-copy recombinants; Sergei called this an absolutely clear SINE, so scoring it high is correct", "good"),
        ("SIMSCRAM__swap__f100", "segments swapped between slots - the mosaic Sergei means", "edge"),
        ("SIMDEL__one__f050", "half the copies missing the same 50 bp block", "edge"),
        ("SIMNEST__f050", "half the copies sharing a host flank", "edge"),
        ("SIMTRUNC__r060", "simulated heavy truncation", "bad"),
        ("NEGLINEORF__teu__r00", "a REAL LINE, queried from its interior — scores 54, "
                                 "no TSDs at all, flanks shared", "bad"),
        ("NEGLINE__teu__r00", "queried with a LINE 3-prime terminus, but 66 % of these loci "
                              "are known Tal SINEs — the shared 3-prime end. Scoring it 99 "
                              "is correct", "edge"),
        ("ERI__eri__e1-1", "a second SINE family: hedgehog, clean", "good"),
        ("ERI__eri__e1-4", "hedgehog subfamily whose copies all share flanking "
                           "sequence — not independent insertions", "bad"),
        ("ERI__eri__e2-3", "hedgehog, partly nested", "edge"),
    ]
    # Human and Timema picks live in other corpora, and their features are not
    # in features_all.jsonl (which only covers aln_v2), so compute them here.
    PICKS = PICKS + EXTRA_PICKS
    out, seen = [], set()
    for item in PICKS:
        st, why, tag = item[0], item[1], item[2]
        srcdir = item[3] if len(item) > 3 else "aln_v2"
        if st in seen:
            continue
        f = os.path.join(srcdir, st + ".aln.fa")
        if not os.path.exists(f):
            print("  missing:", st)
            continue
        seen.add(st)
        d = load(f)
        d.update({"set": st, "why": why, "tag": tag})
        if st in feats:
            d["feat"] = {k: v for k, v in feats[st].items() if isinstance(v, (int, float))}
        else:
            import measure_c as _M
            try:
                mv = _M.measure(f)
                d["feat"] = {k: v for k, v in mv.items() if isinstance(v, (int, float))}
            except Exception:
                d["feat"] = {}
        out.append(d)
    json.dump(out, open("viewer_data.json", "w"))
    print("packed %d alignments, %.0f KB"
          % (len(out), os.path.getsize("viewer_data.json") / 1024))


if __name__ == "__main__":
    main()
