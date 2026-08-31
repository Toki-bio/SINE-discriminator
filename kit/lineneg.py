#!/usr/bin/env python3
"""LINE 3' ends as a natural negative, found the way they really are found.

Only 40 full-length LINE loci exist in this genome, but a search with the LINE's
3' terminus finds thousands of fragments - which is precisely the situation the
spec warns about, where a SINE-like query pulls in LINE 3' ends. That is the
honest negative: a real, abundant, ancient repeat family that is not a SINE.
"""
import os, sys, subprocess, random, re
B = "/data/W/toki/Genomes/Mammalia/Eulipotyphla"
OUT = "/data/W/toki/SINE_disc"
FL = 70
rng = random.Random(77)


def sh(c):
    r = subprocess.run(c, shell=True, capture_output=True, text=True)
    if r.returncode:
        sys.stderr.write("FAIL %s\n%s\n" % (c, r.stderr[-600:]))
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


g = B + "/teu/GCA_964194135.1_mTalEur1.hap1.1_genomic.fna"
lin = list(read_fa(B + "/teu/line/LINE_candidates.fa").values())[0]
tail = lin[-260:]
os.makedirs(OUT + "/tmp", exist_ok=True)
q = OUT + "/tmp/line_tail.fa"
open(q, "w").write(">LINE3p\n%s\n" % tail)
if not os.path.exists(g + ".nsq"):
    sh("makeblastdb -in %s -dbtype nucl" % g)
sh("blastn -query %s -db %s -outfmt '6 sseqid sstart send pident length bitscore' "
   "-evalue 1e-10 -num_threads 16 -max_target_seqs 20000 > %s/tmp/line_hits.tsv"
   % (q, g, OUT))
rows = [l.split("\t") for l in open(OUT + "/tmp/line_hits.tsv") if l.strip()]
rows = [r for r in rows if int(r[4]) >= 150]
print("LINE 3' hits with >=150 bp: %d" % len(rows))
sh("cut -f1,2 %s.fai > %s/tmp/teu.gsize" % (g, OUT))
if len(rows) > 600:
    rows = rng.sample(rows, 600)
bd = OUT + "/tmp/lineneg.bed"
with open(bd, "w") as fh:
    for sid, s, e, pid, ln, bs in rows:
        s, e = int(s), int(e)
        st = "+" if s < e else "-"
        a, b = (s - 1, e) if s < e else (e - 1, s)
        fh.write("%s\t%d\t%d\tLINE\t0\t%s\n" % (sid, max(0, a), b, st))
sh("bedtools slop -i %s -g %s/tmp/teu.gsize -b %d > %s.fl" % (bd, OUT, FL, bd))
sh("bedtools getfasta -fi %s -bed %s.fl -s -name+ -fo %s.fa" % (g, bd, bd))
seqs = [(k, v) for k, v in read_fa(bd + ".fa").items()
        if len(v) >= 200 and v.upper().count("N") < 0.1 * len(v)]
print("usable LINE loci: %d" % len(seqs))
made = 0
for i in range(5):
    if len(seqs) < 60:
        break
    with open("%s/sets_c/NEGLINE__teu__r%02d.fa" % (OUT, i), "w") as fh:
        fh.write(">CONSENSUS_LINE3p\n%s\n" % tail)
        for k, v in rng.sample(seqs, min(100, len(seqs))):
            fh.write(">%s\n%s\n" % (k, v))
    made += 1
print("NEGLINE sets: %d" % made)
