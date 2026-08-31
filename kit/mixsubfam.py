#!/usr/bin/env python3
"""Build the missing negative: one set drawn from TWO structurally different
subfamilies of the same species, searched with only the first one's consensus.

This is Sergei's case "looks like a SINE but with subfamilies of different
structure which need to be separated". Nothing in the corpus tested it, so the
SUBFAMILY_STRUCTURE flag had never had the chance to fire or to be wrong.
"""
import os, random, glob, sys
OUT = "/data/W/toki/SINE_disc"
rng = random.Random(7)


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


PAIRS = [("saq", "s1_30seqs", "s8_225seqs"), ("saq", "s3_43seqs", "s5_5seqs"),
         ("ccr", "g1_180seqs", "g6_58seqs"), ("ccr", "g2_103seqs", "g5_7seqs"),
         ("teu", "t1_45seqs", "t5_31seqs"), ("teu", "t3_27seqs", "t6_324seqs"),
         ("dmo", "d1_16seqs", "d4_266seqs"), ("dmo", "d3_38seqs", "d5_268seqs")]

os.makedirs(OUT + "/sets_c", exist_ok=True)
made = 0
for sp, a, b in PAIRS:
    fa = "%s/sets_c/POS__%s__%s.fa" % (OUT, sp, a)
    fb = "%s/sets_c/POS__%s__%s.fa" % (OUT, sp, b)
    if not (os.path.exists(fa) and os.path.exists(fb)):
        sys.stderr.write("missing %s or %s\n" % (fa, fb))
        continue
    A, B = read(fa), read(fb)
    consA = [(k, v) for k, v in A.items() if k.startswith("CONSENSUS_")]
    ca = [(k, v) for k, v in A.items() if not k.startswith("CONSENSUS_")]
    cb = [(k, v) for k, v in B.items() if not k.startswith("CONSENSUS_")]
    if not consA or len(ca) < 50 or len(cb) < 50:
        continue
    rec = [consA[0]] + rng.sample(ca, 50) + rng.sample(cb, 50)
    with open("%s/sets_c/MIXSUBFAM__%s__%s_%s.fa" % (OUT, sp, a, b), "w") as fh:
        for k, v in rec:
            fh.write(">%s\n%s\n" % (k, v))
    made += 1
print("wrote %d mixed-subfamily sets" % made)
