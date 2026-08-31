#!/usr/bin/env python3
"""Two things the corpus has been missing: a NATURAL negative, and a SECOND
SINE family.

NEGLINE  Real LINE 3' ends from Talpa europaea. The spec names LINE 3' ends as a
         confusable class, and until now every negative except random loci was
         synthetic or carved from SINE data. A LINE is a real, old, abundant
         repeat family that is genuinely NOT a SINE - the hardest honest test.

ERI      Erinaceus europaeus (hedgehog), 8 curated subfamilies, run through the
         identical pipeline. Every biological number so far comes from one SINE
         family type; this is the check on which of them generalise.
"""
import os, re, sys, subprocess, random
from collections import defaultdict

B = "/data/W/toki/Genomes/Mammalia/Eulipotyphla"
OUT = "/data/W/toki/SINE_disc"
FL = 70
rng = random.Random(31)
NAME = re.compile(r"^(.+):(\d+)-(\d+)\(([+-])\)$")


def sh(c):
    r = subprocess.run(c, shell=True, capture_output=True, text=True)
    if r.returncode:
        sys.stderr.write("FAIL %s\n%s\n" % (c, r.stderr[-800:]))
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


def emit(name, cons, recs):
    with open("%s/sets_c/%s.fa" % (OUT, name), "w") as fh:
        fh.write(">CONSENSUS_%s\n%s\n" % (name.split("__")[-1], cons))
        for k, v in recs:
            fh.write(">%s\n%s\n" % (k, v))


def line_negative():
    g = B + "/teu/GCA_964194135.1_mTalEur1.hap1.1_genomic.fna"
    bed = B + "/teu/line/best_loci/best_loci_main.bed6"
    cand = read_fa(B + "/teu/line/LINE_candidates.fa")
    lin = list(cand.values())[0]
    tail = lin[-250:]                     # the 3' end, which is what a search finds
    if not os.path.exists(g + ".fai"):
        sh("samtools faidx " + g)
    sh("cut -f1,2 %s.fai > %s/tmp/teu.gsize" % (g, OUT))
    rows = [l.split("\t") for l in open(bed) if l.strip()]
    rows = [r for r in rows if len(r) >= 6 and int(r[2]) - int(r[1]) > 400]
    if len(rows) > 400:
        rows = rng.sample(rows, 400)
    tb = OUT + "/tmp/line.bed"
    with open(tb, "w") as fh:
        for r in rows:
            st, en, sd = int(r[1]), int(r[2]), r[5]
            # the 3'-most 250 bp of the element, in its own orientation
            a, b = (en - 250, en) if sd == "+" else (st, st + 250)
            fh.write("%s\t%d\t%d\tLINE\t0\t%s\n" % (r[0], max(0, a), b, sd))
    sh("bedtools slop -i %s -g %s/tmp/teu.gsize -b %d > %s.fl" % (tb, OUT, FL, tb))
    sh("bedtools getfasta -fi %s -bed %s.fl -s -name+ -fo %s.fa" % (g, tb, tb))
    seqs = [(k, v) for k, v in read_fa(tb + ".fa").items()
            if len(v) >= 250 and v.upper().count("N") < 0.1 * len(v)]
    print("LINE loci extracted: %d" % len(seqs))
    for i in range(4):
        if len(seqs) < 60:
            break
        emit("NEGLINE__teu__r%02d" % i, tail, rng.sample(seqs, min(100, len(seqs))))
    return len(seqs)


def eri():
    R = B + "/Erniacidae/run_20260820_221537"
    g = B + "/Erniacidae/GCF_950295315.1_mEriEur2.1_genomic.fna"
    if not os.path.exists(g + ".fai"):
        sh("samtools faidx " + g)
    sh("cut -f1,2 %s.fai > %s/tmp/eri.gsize" % (g, OUT))
    cons = read_fa(R + "/results/consensuses.fa")
    bysub = defaultdict(list)
    with open(R + "/results/assignment_full.tsv") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[4] != "assigned":
                continue
            m = NAME.match(f[0].replace("@U@", "_"))
            if m:
                bysub[f[1]].append((m.group(1), int(m.group(2)), int(m.group(3)),
                                    m.group(4), int(f[2])))
    made = 0
    for sub, rows in sorted(bysub.items()):
        if sub not in cons or len(rows) < 150:
            continue
        rows = sorted(rows, key=lambda r: r[4])
        pick = [rows[len(rows) * i // 10 + j] for i in range(10)
                for j in range(min(12, len(rows) // 10))][:120]
        bd = "%s/tmp/eri.%s.bed" % (OUT, sub)
        with open(bd, "w") as fh:
            for c, s, e, st, bs in pick:
                fh.write("%s\t%d\t%d\t%s\t0\t%s\n" % (c, s, e, sub, st))
        sh("bedtools slop -i %s -g %s/tmp/eri.gsize -b %d > %s.fl" % (bd, OUT, FL, bd))
        sh("bedtools getfasta -fi %s -bed %s.fl -s -name+ -fo %s.fa" % (g, bd, bd))
        seqs = [(k, v) for k, v in read_fa(bd + ".fa").items()
                if len(v) >= 2 * FL + 60 and v.upper().count("N") < 0.1 * len(v)]
        if len(seqs) < 40:
            continue
        emit("ERI__eri__%s" % sub, cons[sub], seqs[:100])
        made += 1
    print("eri subfamilies emitted: %d" % made)
    return made


if __name__ == "__main__":
    os.makedirs(OUT + "/tmp", exist_ok=True)
    line_negative()
    eri()
