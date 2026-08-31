#!/usr/bin/env python3
"""A satellite negative, sampled the way a satellite actually presents.

The blast route failed because a 250 bp window spans several monomers and comes
back as many short hits. But the defining property of a satellite does not need
a search at all: the copies are ADJACENT, in one array, not dispersed insertions
anywhere in the genome. So sample windows along a single long tandem array.

This is the case the spec calls "satellites, cheaply excluded by neighbour_hit +
coordinate clustering" - and the point of building it is to check whether they
actually are.
"""
import os, sys, subprocess, random
B = "/data/W/toki/Genomes/Mammalia/Eulipotyphla"
OUT = "/data/W/toki/SINE_disc"
FL = 70
rng = random.Random(808)
G = B + "/Erniacidae/GCF_950295315.1_mEriEur2.1_genomic.fna"
TRF = B + "/Erniacidae/run_20260820_221537/trf_out/tandem_repeats.merged.bed"


def sh(c):
    r = subprocess.run(c, shell=True, capture_output=True, text=True)
    if r.returncode:
        sys.stderr.write("FAIL %s\n" % r.stderr[-300:])
    return r.stdout


def read_fa(p):
    d, n = {}, None
    for line in open(p):
        line = line.rstrip()
        if line.startswith(">"):
            n = line[1:].split()[0]
            d[n] = []
        elif n:
            d[n].append(line)
    return dict((k, "".join(v)) for k, v in d.items())


ivs = []
with open(TRF) as fh:
    for line in fh:
        f = line.split("\t")
        if len(f) >= 3 and int(f[2]) - int(f[1]) > 30000:
            ivs.append((f[0].replace("@U@", "_"), int(f[1]), int(f[2])))  # TRF bed encodes NC@U@080185.1
print("arrays longer than 30 kb: %d" % len(ivs))
made = 0
for c, s, e in ivs[:40]:
    n = min(100, (e - s) // 300)
    if n < 60:
        continue
    bd = OUT + "/tmp/sat%d.bed" % made
    with open(bd, "w") as fh:
        for i in range(n):
            p = s + 50 + i * ((e - s - 400) // n)
            fh.write("%s\t%d\t%d\tSAT\t0\t+\n" % (c, p, p + 250))
    sh("cut -f1,2 %s.fai > %s/tmp/eri.gsize" % (G, OUT))
    sh("bedtools slop -i %s -g %s/tmp/eri.gsize -b %d > %s.fl" % (bd, OUT, FL, bd))
    sh("bedtools getfasta -fi %s -bed %s.fl -s -name+ -fo %s.fa" % (G, bd, bd))
    seqs = [(k, v) for k, v in read_fa(bd + ".fa").items()
            if len(v) >= 300 and v.upper().count("N") < 0.1 * len(v)]
    if len(seqs) < 60:
        continue
    cons = seqs[0][1][FL:FL + 250]
    with open("%s/sets_c/NEGSAT__eri__r%02d.fa" % (OUT, made), "w") as fh:
        fh.write(">CONSENSUS_SAT\n%s\n" % cons)
        for k, v in seqs[1:101]:
            fh.write(">%s\n%s\n" % (k, v))
    print("  NEGSAT r%02d  %s:%d-%d  %d windows" % (made, c, s, e, len(seqs)))
    made += 1
    if made >= 3:
        break
print("satellite sets:", made)
