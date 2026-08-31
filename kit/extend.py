#!/usr/bin/env python3
"""Re-extract with 400 bp flanks, body-aligned only.

Sergei's repeated point: several sets look bad at an edge, but the flank is too
short to tell whether the edge is real or just where the extraction stopped. The
70 bp limit existed only to stop MAFFT smearing unrelated flank across the
alignment - and since flanks are now JUSTIFIED rather than aligned, that
constraint is gone.

Architecture (the same one the Tal repo uses): align the element body with the
consensus, then attach the flanks as unaligned justified blocks. Flank length
then costs nothing at all.
"""
import os, re, sys, subprocess
B = "/data/W/toki/Genomes/Mammalia/Eulipotyphla"
OUT = "/data/W/toki/SINE_disc"
BIG = 400
GEN = {"eri": B + "/Erniacidae/GCF_950295315.1_mEriEur2.1_genomic.fna",
       "teu": B + "/teu/GCA_964194135.1_mTalEur1.hap1.1_genomic.fna",
       "ccr": B + "/ccr/GCF_000260355.1_ConCri1.0_genomic.fna",
       "saq": B + "/saq/GCA_004024925.1_ScaAqu_v1_BIUU_genomic.fna"}
LOC = re.compile(r"([^:]+):(\d+)-(\d+)\(([+-])\)")


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


def go(setname, sp):
    src = "%s/sets_c/%s.fa" % (OUT, setname)
    if not os.path.exists(src):
        print("  missing", setname)
        return
    old = read_fa(src)
    cons = [v for k, v in old.items() if k.startswith("CONSENSUS_")][0]
    g = GEN[sp]
    if not os.path.exists(g + ".fai"):
        sh("samtools faidx " + g)
    sh("cut -f1,2 %s.fai > %s/tmp/%s.gsize" % (g, OUT, sp))
    rows = []
    for k in old:
        if k.startswith("CONSENSUS_"):
            continue
        m = LOC.search(k)
        if m:
            rows.append((m.group(1), int(m.group(2)), int(m.group(3)), m.group(4), k))
    if len(rows) < 20:
        print("  %s: only %d parsable names" % (setname, len(rows)))
        return
    bd = "%s/tmp/ext_%s.bed" % (OUT, setname)
    with open(bd, "w") as fh:
        for c, s, e, st, k in rows:
            fh.write("%s\t%d\t%d\tE\t0\t%s\n" % (c, s, e, st))
    # body only, for the alignment
    sh("bedtools getfasta -fi %s -bed %s -s -name+ -fo %s.body.fa" % (g, bd, bd))
    # body plus 400 bp either side, for the flanks
    sh("bedtools slop -i %s -g %s/tmp/%s.gsize -b %d > %s.fl" % (bd, OUT, sp, BIG, bd))
    sh("bedtools getfasta -fi %s -bed %s.fl -s -name+ -fo %s.full.fa" % (g, bd, bd))
    body = list(read_fa(bd + ".body.fa").items())
    full = list(read_fa(bd + ".full.fa").items())
    if len(body) != len(full):
        print("  %s: body/full mismatch" % setname)
        return
    # align the body with the consensus - no flank in the aligner at all
    tmp = "%s/tmp/ext_%s.in.fa" % (OUT, setname)
    with open(tmp, "w") as fh:
        fh.write(">CONSENSUS_ext\n%s\n" % cons)
        for i, (k, v) in enumerate(body):
            fh.write(">b%04d\n%s\n" % (i, v))
    sh("mafft --retree 2 --maxiterate 0 --adjustdirection --quiet --thread 12 %s > %s.aln"
       % (tmp, tmp))
    al = read_fa(tmp + ".aln")
    out = "%s/aln_ext/%s.aln.fa" % (OUT, setname)
    os.makedirs(OUT + "/aln_ext", exist_ok=True)
    lf, rf = {}, {}
    for i, (k, v) in enumerate(full):
        b = body[i][1]
        j = v.find(b)
        if j < 0:
            lf[i], rf[i] = "", ""
        else:
            lf[i], rf[i] = v[:j], v[j + len(b):]
    wL = max((len(x) for x in lf.values()), default=0)
    wR = max((len(x) for x in rf.values()), default=0)
    with open(out, "w") as fh:
        ck = [k for k in al if "CONSENSUS" in k.upper()]
        fh.write(">CONSENSUS_ext\n%s%s%s\n" % ("-" * wL, al[ck[0]], "-" * wR))
        for i, (k, v) in enumerate(body):
            key = "b%04d" % i
            key = key if key in al else "_R_" + key
            if key not in al:
                continue
            fh.write(">%s\n%s%s%s\n" % (k, lf[i].rjust(wL, "-"), al[key],
                                        rf[i].ljust(wR, "-")))
    print("  wrote %s  (%d copies, flanks up to %d/%d bp)" % (setname, len(body), wL, wR))


if __name__ == "__main__":
    import glob
    # every set whose copies are REAL genomic loci; synthetic sets have
    # synthetic flanks and nothing to extend into
    todo = []
    for p in sorted(glob.glob(OUT + "/sets_c/*.fa")):
        b = os.path.basename(p)[:-3]
        pref = b.split("__")[0]
        if pref not in ("POS", "ERI", "NEGRAND", "NEGSAT", "NEGSEGDUP", "NEGLINEORF"):
            continue
        parts = b.split("__")
        sp = parts[1] if len(parts) > 1 else None
        if sp in GEN:
            todo.append((b, sp))
    print("extending %d sets" % len(todo))
    for s, sp in todo:
        go(s, sp)
