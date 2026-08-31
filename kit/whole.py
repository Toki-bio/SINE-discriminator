#!/usr/bin/env python3
"""Align the WHOLE sequence where the flanks are not really flanks.

Sergei: "you provided badly aligned files, they should be aligned wholly because
there are no true flanks."

He is right, and it exposes an assumption baked into the v2 justify step: that
flanking sequence is non-homologous and so must not be aligned. That holds for a
SINE, whose neighbours are unrelated genomic DNA. It is exactly false for a
satellite, a segmental duplication, or a LINE fragment - there the "flank" is
more of the same element, and justifying it splits a continuous alignment in two.

Rule: when flank decay says NOT_ISOLATED or FRAGMENT_OF_LONGER, the flanks are
part of the element, so align end to end.
"""
import os, sys, subprocess, json
OUT = "/data/W/toki/SINE_disc"


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


SETS = ["NEGSAT__eri__r00", "NEGSAT__eri__r01", "NEGSAT__eri__r02",
        "NEGSEGDUP__eri__r00", "NEGSEGDUP__eri__r01", "NEGSEGDUP__eri__r02",
        "NEGLINEORF__teu__r00", "NEGLINEORF__teu__r01",
        "NEGLINEORF__teu__r02", "NEGLINEORF__teu__r03",
        "ERI__eri__e1-4", "ERI__eri__e2-2", "ERI__eri__e2-3",
        "POS__saq__s5_5seqs"]

os.makedirs(OUT + "/aln_whole", exist_ok=True)
for s in SETS:
    src = "%s/aln_ext/%s.aln.fa" % (OUT, s)
    if not os.path.exists(src):
        print("  no extended alignment for", s)
        continue
    d = read_fa(src)
    tmp = "%s/tmp/whole_%s.fa" % (OUT, s)
    with open(tmp, "w") as fh:
        for k, v in d.items():
            u = v.replace("-", "")
            if len(u) > 40:
                fh.write(">%s\n%s\n" % (k, u))
    out = "%s/aln_whole/%s.aln.fa" % (OUT, s)
    sh("mafft --retree 2 --maxiterate 0 --adjustdirection --quiet --thread 12 %s > %s"
       % (tmp, out))
    n = sum(1 for l in open(out) if l.startswith(">"))
    print("  %-26s %d sequences aligned end to end" % (s, n))
