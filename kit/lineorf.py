#!/usr/bin/env python3
"""A genuine LINE negative, from a region SINEs cannot share.

The previous attempt queried the LINE's 3' TERMINUS - which is precisely the
segment a SINE mobilised by that LINE also carries. 66 % of the hits overlapped
known Tal SINE loci, so the set was two-thirds real SINEs and scoring it 98 was
correct behaviour, not a failure.

This queries the LINE interior instead (~7 kb in, ORF territory), which no SINE
carries. Overlap with the known SINE annotation is reported, so the set can be
trusted or discarded on evidence rather than assumption.
"""
import os, sys, subprocess, random
B = "/data/W/toki/Genomes/Mammalia/Eulipotyphla"
OUT = "/data/W/toki/SINE_disc"
FL = 70
rng = random.Random(99)


def sh(c):
    r = subprocess.run(c, shell=True, capture_output=True, text=True)
    if r.returncode:
        sys.stderr.write("FAIL %s\n%s\n" % (c, r.stderr[-500:]))
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
print("LINE length %d" % len(lin))
mid = lin[7000:7260]                      # interior: ORF2 region, not shared with SINEs
q = OUT + "/tmp/line_orf.fa"
open(q, "w").write(">LINE_ORF\n%s\n" % mid)
sh("blastn -query %s -db %s -outfmt '6 sseqid sstart send pident length bitscore' "
   "-evalue 1e-5 -num_threads 16 -max_target_seqs 20000 > %s/tmp/orf_hits.tsv" % (q, g, OUT))
rows = [l.split("\t") for l in open(OUT + "/tmp/orf_hits.tsv") if l.strip()]
rows = [r for r in rows if int(r[4]) >= 60]
print("LINE-interior hits >=60 bp: %d" % len(rows))
if len(rows) < 60:
    sys.exit("too few interior hits")
if len(rows) > 600:
    rows = rng.sample(rows, 600)
bd = OUT + "/tmp/lineorf.bed"
with open(bd, "w") as fh:
    for sid, s, e, pid, ln, bs in rows:
        s, e = int(s), int(e)
        st = "+" if s < e else "-"
        a, b = (s - 1, e) if s < e else (e - 1, s)
        fh.write("%s\t%d\t%d\tLINEORF\t0\t%s\n" % (sid, max(0, a), b, st))
sh("sort -k1,1 -k2,2n %s > %s.s" % (bd, bd))
ov = sh("bedtools intersect -a %s.s -b %s/tmp/teu_sines.bed -u -f 0.5 | wc -l" % (bd, OUT)).strip()
print("of %d interior hits, overlapping a known Tal SINE: %s" % (len(rows), ov))
sh("bedtools slop -i %s.s -g %s/tmp/teu.gsize -b %d > %s.fl" % (bd, OUT, FL, bd))
sh("bedtools getfasta -fi %s -bed %s.fl -s -name+ -fo %s.fa" % (g, bd, bd))
seqs = [(k, v) for k, v in read_fa(bd + ".fa").items()
        if len(v) >= 2 * FL + 40 and v.upper().count("N") < 0.1 * len(v)]
print("usable: %d" % len(seqs))
for i in range(4):
    if len(seqs) < 60:
        break
    with open("%s/sets_c/NEGLINEORF__teu__r%02d.fa" % (OUT, i), "w") as fh:
        fh.write(">CONSENSUS_LINEORF\n%s\n" % mid)
        for k, v in rng.sample(seqs, min(100, len(seqs))):
            fh.write(">%s\n%s\n" % (k, v))
print("NEGLINEORF sets written")
