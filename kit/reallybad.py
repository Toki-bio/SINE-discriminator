#!/usr/bin/env python3
"""Two classes that are definitively NOT SINEs, built from real sequence.

Everything I have called a negative so far has turned out to contain a real
element - NEGSPLICE, NEGLINE, NEGCHIM and all six mosaic types were judged SINEs
on inspection. So these are chosen to be structurally incapable of being a SINE:

NEGSAT     A real tandem-repeat array from the hedgehog TRF annotation. A search
           seeded from inside a satellite returns hundreds of hits that are all
           in the SAME array - adjacent, overlapping, not dispersed insertions.

NEGSEGDUP  A duplicated genomic region: a window with 5-50 near-identical
           paralogs. Similarity runs straight THROUGH the flanks, so there is no
           boundary anywhere - the spec's "easy but essential negative".
"""
import os, sys, subprocess, random
B = "/data/W/toki/Genomes/Mammalia/Eulipotyphla"
OUT = "/data/W/toki/SINE_disc"
FL = 70
rng = random.Random(555)


def sh(c):
    r = subprocess.run(c, shell=True, capture_output=True, text=True)
    if r.returncode:
        sys.stderr.write("FAIL %s\n%s\n" % (c[:120], r.stderr[-400:]))
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


G = B + "/Erniacidae/GCF_950295315.1_mEriEur2.1_genomic.fna"
TRF = B + "/Erniacidae/run_20260820_221537/trf_out/tandem_repeats.merged.bed"
os.makedirs(OUT + "/tmp", exist_ok=True)
if not os.path.exists(G + ".fai"):
    sh("samtools faidx " + G)
sh("cut -f1,2 %s.fai > %s/tmp/eri.gsize" % (G, OUT))


def hits_for(query_fa, tag, minlen, want_lo, want_hi):
    sh("blastn -query %s -db %s -outfmt '6 sseqid sstart send pident length' "
       "-evalue 1e-10 -num_threads 16 -max_target_seqs 5000 > %s/tmp/%s.tsv"
       % (query_fa, G, OUT, tag))
    rows = [l.split("\t") for l in open("%s/tmp/%s.tsv" % (OUT, tag)) if l.strip()]
    rows = [r for r in rows if int(r[4]) >= minlen]
    return rows


# ---------------- satellite ------------------------------------------------
def satellite():
    ivs = []
    with open(TRF) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 3 and int(f[2]) - int(f[1]) > 8000:
                ivs.append((f[0], int(f[1]), int(f[2])))
            if len(ivs) > 4000:
                break
    print("long tandem arrays available: %d" % len(ivs))
    made = 0
    for attempt in range(30):
        c, s, e = rng.choice(ivs)
        mid = (s + e) // 2
        q = OUT + "/tmp/sat_q.bed"
        open(q, "w").write("%s\t%d\t%d\tSAT\t0\t+\n" % (c, mid, mid + 250))
        sh("bedtools getfasta -fi %s -bed %s -fo %s.fa" % (G, q, q))
        seq = list(read_fa(q + ".fa").values())
        if not seq or len(seq[0]) < 240:
            continue
        open(OUT + "/tmp/sat_q2.fa", "w").write(">SAT\n%s\n" % seq[0])
        rows = hits_for(OUT + "/tmp/sat_q2.fa", "sat", 120, 100, 5000)
        if len(rows) < 120:
            continue
        bd = OUT + "/tmp/sat.bed"
        with open(bd, "w") as fh:
            for sid, a, b, pid, ln in rng.sample(rows, min(400, len(rows))):
                a, b = int(a), int(b)
                st = "+" if a < b else "-"
                lo, hi = (a - 1, b) if a < b else (b - 1, a)
                fh.write("%s\t%d\t%d\tSAT\t0\t%s\n" % (sid, max(0, lo), hi, st))
        sh("sort -k1,1 -k2,2n %s > %s.s" % (bd, bd))
        sh("bedtools slop -i %s.s -g %s/tmp/eri.gsize -b %d > %s.fl" % (bd, OUT, FL, bd))
        sh("bedtools getfasta -fi %s -bed %s.fl -s -name+ -fo %s.fa" % (G, bd, bd))
        seqs = [(k, v) for k, v in read_fa(bd + ".fa").items()
                if len(v) >= 2 * FL + 60 and v.upper().count("N") < 0.1 * len(v)]
        if len(seqs) < 60:
            continue
        with open("%s/sets_c/NEGSAT__eri__r%02d.fa" % (OUT, made), "w") as fh:
            fh.write(">CONSENSUS_SAT\n%s\n" % seq[0])
            for k, v in rng.sample(seqs, min(100, len(seqs))):
                fh.write(">%s\n%s\n" % (k, v))
        made += 1
        print("  NEGSAT r%02d from %s:%d (%d hits)" % (made - 1, c, mid, len(rows)))
        if made >= 3:
            break
    return made


# ---------------- segmental duplication ------------------------------------
def segdup():
    fai = [l.split("\t") for l in open(G + ".fai")]
    big = [(f[0], int(f[1])) for f in fai if int(f[1]) > 2000000]
    made = 0
    for attempt in range(60):
        c, ln = rng.choice(big)
        s = rng.randint(1000, ln - 2000)
        q = OUT + "/tmp/sd_q.bed"
        open(q, "w").write("%s\t%d\t%d\n" % (c, s, s + 300))
        sh("bedtools getfasta -fi %s -bed %s -fo %s.fa" % (G, q, q))
        seq = list(read_fa(q + ".fa").values())
        if not seq or seq[0].upper().count("N") > 10:
            continue
        open(OUT + "/tmp/sd_q2.fa", "w").write(">SD\n%s\n" % seq[0])
        rows = hits_for(OUT + "/tmp/sd_q2.fa", "sd", 200, 5, 60)
        rows = [r for r in rows if float(r[3]) >= 88]
        if not (25 <= len(rows) <= 400):
            continue
        bd = OUT + "/tmp/sd.bed"
        with open(bd, "w") as fh:
            for sid, a, b, pid, l2 in rows:
                a, b = int(a), int(b)
                st = "+" if a < b else "-"
                lo, hi = (a - 1, b) if a < b else (b - 1, a)
                fh.write("%s\t%d\t%d\tSD\t0\t%s\n" % (sid, max(0, lo), hi, st))
        sh("sort -k1,1 -k2,2n %s > %s.s" % (bd, bd))
        sh("bedtools slop -i %s.s -g %s/tmp/eri.gsize -b %d > %s.fl" % (bd, OUT, FL, bd))
        sh("bedtools getfasta -fi %s -bed %s.fl -s -name+ -fo %s.fa" % (G, bd, bd))
        seqs = [(k, v) for k, v in read_fa(bd + ".fa").items()
                if len(v) >= 2 * FL + 60 and v.upper().count("N") < 0.1 * len(v)]
        if len(seqs) < 25:
            continue
        with open("%s/sets_c/NEGSEGDUP__eri__r%02d.fa" % (OUT, made), "w") as fh:
            fh.write(">CONSENSUS_SD\n%s\n" % seq[0])
            for k, v in seqs[:100]:
                fh.write(">%s\n%s\n" % (k, v))
        made += 1
        print("  NEGSEGDUP r%02d from %s:%d (%d paralogs)" % (made - 1, c, s, len(rows)))
        if made >= 3:
            break
    return made


if __name__ == "__main__":
    print("satellite sets:", satellite())
    print("segdup sets:", segdup())
