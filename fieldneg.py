#!/usr/bin/env python3
"""Negatives and grey cases from a real run, rather than constructed by me.

Every negative I have built has needed Sergei to check whether it was actually a
SINE, and several were. The curated SINEderella runs already contain loci the
pipeline itself rejected, with a stated reason:

  NEGLOWBIT   rejected_low_bitscore - search score ~223 against ~1770 for an
              assigned copy. Real loci a real search returned and a curated
              pipeline threw away.
  AMBIG       no_unanimous_votes - the copy matched, but the subfamily vote was
              split. Probably real SINEs of uncertain subfamily, so a
              field-realistic GREY case rather than a negative.

Extracted with 400 bp flanks, body aligned with the consensus, flanks attached
unaligned - the same geometry as everything else, so the numbers are comparable.
"""
import os, sys, subprocess, random
B = "/data/W/toki/Genomes/Mammalia/Eulipotyphla"
OUT = "/data/W/toki/SINE_disc"
BIG = 400
rng = random.Random(1234)
RUNS = {"saq": ("saq/run_20260425_182219", "saq/GCA_004024925.1_ScaAqu_v1_BIUU_genomic.fna"),
        "teu": ("teu/run_20260427_130055", "teu/GCA_964194135.1_mTalEur1.hap1.1_genomic.fna")}


def sh(c):
    r = subprocess.run(c, shell=True, capture_output=True, text=True)
    if r.returncode:
        sys.stderr.write("FAIL %s\n" % r.stderr[-250:])
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


def build(sp, reason, tag, nsets):
    run, gen = RUNS[sp]
    R, g = B + "/" + run, B + "/" + gen
    cons = read_fa(R + "/results/consensuses.fa")
    rows = []
    with open(R + "/results/unassigned.tsv") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ir = hdr.index("Reason")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > ir and f[ir] == reason:
                rows.append((f[1], int(f[2]), int(f[3]), f[4], f[5]))
    print("  %s %s: %d loci" % (sp, reason, len(rows)))
    if len(rows) < 120:
        return 0
    sh("cut -f1,2 %s.fai > %s/tmp/%s.gsize" % (g, OUT, sp))
    made = 0
    for i in range(nsets):
        pick = rng.sample(rows, 100)
        sub = max(set(p[4] for p in pick), key=[p[4] for p in pick].count)
        if sub not in cons:
            continue
        bd = "%s/tmp/%s_%s_%d.bed" % (OUT, tag, sp, i)
        with open(bd, "w") as fh:
            for c, s, e, st, _ in pick:
                fh.write("%s\t%d\t%d\tX\t0\t%s\n" % (c, s, e, st))
        sh("bedtools getfasta -fi %s -bed %s -s -name+ -fo %s.body.fa" % (g, bd, bd))
        sh("bedtools slop -i %s -g %s/tmp/%s.gsize -b %d > %s.fl" % (bd, OUT, sp, BIG, bd))
        sh("bedtools getfasta -fi %s -bed %s.fl -s -name+ -fo %s.full.fa" % (g, bd, bd))
        body = list(read_fa(bd + ".body.fa").items())
        full = list(read_fa(bd + ".full.fa").items())
        if len(body) != len(full) or len(body) < 60:
            continue
        tmp = "%s/tmp/%s_%s_%d.in.fa" % (OUT, tag, sp, i)
        with open(tmp, "w") as fh:
            fh.write(">CONSENSUS_%s\n%s\n" % (sub, cons[sub]))
            for j, (k, v) in enumerate(body):
                fh.write(">b%04d\n%s\n" % (j, v))
        sh("mafft --retree 2 --maxiterate 0 --adjustdirection --quiet --thread 12 %s > %s.aln"
           % (tmp, tmp))
        al = read_fa(tmp + ".aln")
        lf, rf = {}, {}
        for j, (k, v) in enumerate(full):
            b = body[j][1]
            p = v.find(b)
            lf[j], rf[j] = (v[:p], v[p + len(b):]) if p >= 0 else ("", "")
        wL = max((len(x) for x in lf.values()), default=0)
        wR = max((len(x) for x in rf.values()), default=0)
        os.makedirs(OUT + "/aln_ext", exist_ok=True)
        name = "%s__%s__r%02d" % (tag, sp, i)
        with open("%s/aln_ext/%s.aln.fa" % (OUT, name), "w") as fh:
            ck = [k for k in al if "CONSENSUS" in k.upper()][0]
            fh.write(">CONSENSUS_%s\n%s%s%s\n" % (sub, "-" * wL, al[ck], "-" * wR))
            for j, (k, v) in enumerate(body):
                key = "b%04d" % j
                key = key if key in al else "_R_" + key
                if key in al:
                    fh.write(">%s\n%s%s%s\n" % (k, lf[j].rjust(wL, "-"), al[key],
                                                rf[j].ljust(wR, "-")))
        print("    wrote %s (consensus %s)" % (name, sub))
        made += 1
    return made


if __name__ == "__main__":
    for sp in ("saq", "teu"):
        build(sp, "rejected_low_bitscore", "NEGLOWBIT", 2)
        build(sp, "no_unanimous_votes", "AMBIG", 2)
