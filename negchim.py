#!/usr/bin/env python3
"""The negative NEGSPLICE was supposed to be.

Every NEGSPLICE set joins two Tal SUBFAMILIES, so it is a within-family
recombinant - a real biological entity that correctly scores as a SINE. It was
never a negative. This builds the actual case: half a real element joined to
non-element sequence, with the breakpoint varying per copy so the result is a
mosaic rather than a coherent new element.

  NEGCHIM   3' half replaced with random genomic sequence, per-copy breakpoint
  NEGMOSAIC two subfamilies, but the breakpoint differs per copy, so different
            copies carry different blocks - the spec's Tier-3 case
"""
import os, random, sys
OUT = "/data/W/toki/SINE_disc"
FL = 70
rng = random.Random(11)


def read(p):
    d, n = {}, None
    for line in open(p):
        line = line.rstrip()
        if line.startswith(">"):
            n = line[1:]
            d[n] = []
        elif n:
            d[n].append(line)
    return {k: "".join(v) for k, v in d.items()}


SP = ["saq", "ccr", "teu", "dmo"]
FAMS = {"saq": ["s1_30seqs", "s3_43seqs", "s5_5seqs", "s8_225seqs"],
        "ccr": ["g1_180seqs", "g3_71seqs", "g5_7seqs", "g6_58seqs"],
        "teu": ["t1_45seqs", "t2_75seqs", "t5_31seqs", "t6_324seqs"],
        "dmo": ["d1_16seqs", "d3_38seqs", "d4_266seqs", "d5_268seqs"]}
made = 0
for sp in SP:
    rand = read("%s/sets_c/NEGRAND__%s__r00.fa" % (OUT, sp))
    rseq = [v for k, v in rand.items() if not k.startswith("CONSENSUS_")]
    for fam in FAMS[sp]:
        f = "%s/sets_c/POS__%s__%s.fa" % (OUT, sp, fam)
        if not os.path.exists(f) or not rseq:
            continue
        A = read(f)
        cons = [(k, v) for k, v in A.items() if k.startswith("CONSENSUS_")]
        cop = [(k, v) for k, v in A.items() if not k.startswith("CONSENSUS_")]
        if not cons or len(cop) < 60:
            continue

        # NEGCHIM: keep the 5' flank + 5' part of the element, then foreign DNA
        rec = [cons[0]]
        for i, (k, v) in enumerate(rng.sample(cop, 100)):
            el = len(v) - 2 * FL
            cut = FL + int(el * rng.uniform(0.35, 0.65))
            donor = rng.choice(rseq)
            off = rng.randint(0, max(0, len(donor) - (len(v) - cut) - 1))
            rec.append((k + "|chim", v[:cut] + donor[off:off + (len(v) - cut)]))
        with open("%s/sets_c/NEGCHIM__%s__%s.fa" % (OUT, sp, fam), "w") as fh:
            for k, v in rec:
                fh.write(">%s\n%s\n" % (k, v))
        made += 1

    # NEGMOSAIC: two subfamilies, breakpoint varying per copy
    for a, b in [(FAMS[sp][0], FAMS[sp][3]), (FAMS[sp][1], FAMS[sp][2])]:
        fa = "%s/sets_c/POS__%s__%s.fa" % (OUT, sp, a)
        fb = "%s/sets_c/POS__%s__%s.fa" % (OUT, sp, b)
        if not (os.path.exists(fa) and os.path.exists(fb)):
            continue
        A, B = read(fa), read(fb)
        ca = [(k, v) for k, v in A.items() if not k.startswith("CONSENSUS_")]
        cb = [(k, v) for k, v in B.items() if not k.startswith("CONSENSUS_")]
        consA = [(k, v) for k, v in A.items() if k.startswith("CONSENSUS_")]
        if not consA or len(ca) < 100 or len(cb) < 100:
            continue
        rec = [consA[0]]
        sa, sb = rng.sample(ca, 100), rng.sample(cb, 100)
        for i in range(100):
            va, vb = sa[i][1], sb[i][1]
            fa_ = rng.uniform(0.2, 0.8)          # breakpoint varies per copy
            ca_ = FL + int((len(va) - 2 * FL) * fa_)
            cb_ = FL + int((len(vb) - 2 * FL) * fa_)
            rec.append(("mos%03d|%s|%s" % (i, a, b), va[:ca_] + vb[cb_:]))
        with open("%s/sets_c/NEGMOSAIC__%s__%s_%s.fa" % (OUT, sp, a, b), "w") as fh:
            for k, v in rec:
                fh.write(">%s\n%s\n" % (k, v))
        made += 1
print("wrote %d sets" % made)
