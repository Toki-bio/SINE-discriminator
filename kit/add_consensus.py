#!/usr/bin/env python3
"""Prepend the known subfamily consensus to every corpus set, then realign.

The boundary is then read directly off the consensus row: the first and last
alignment column at which it carries a nucleotide. No threshold, no smoothing,
no free parameter - and nothing that changes when the copies are resampled.
"""
import os, sys, glob, subprocess

BASE = "/data/W/toki/Genomes/Mammalia/Eulipotyphla"
OUT = "/data/W/toki/SINE_disc"
KEEP = 70                 # bp of flank retained in the alignment
RUNS = {"saq": "saq/run_20260425_182219", "ccr": "ccr/run_20260517_205955",
        "teu": "teu/run_20260427_130055", "dmo": "dmo/run_20260427_103405"}
# stand-in query for the pure-null sets: a real search is always run WITH some
# consensus, so the null case is that consensus finding nothing that supports it
FALLBACK = {"saq": "s7g_172seqs", "ccr": "g3_71seqs",
            "teu": "t5_31seqs", "dmo": "d5_268seqs"}


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    return dict((k, "".join(v)) for k, v in seqs.items())


def main():
    cons = {}
    for sp, run in RUNS.items():
        cons[sp] = read_fasta(os.path.join(BASE, run, "results/consensuses.fa"))
        sys.stderr.write("%s: %d consensuses\n" % (sp, len(cons[sp])))

    src = os.path.join(OUT, "sets")
    dst = os.path.join(OUT, "sets_c")
    os.makedirs(dst, exist_ok=True)
    missing = 0
    for f in sorted(glob.glob(os.path.join(src, "*.fa"))):
        base = os.path.basename(f)[:-3]
        p = base.split("__")
        cls = p[0]
        if cls == "NEGSPLICE":            # "saq_s1_30seqs" -> species, subfamily
            sp, sub = p[1].split("_", 1)
        elif cls == "NEGRAND":
            sp, sub = p[1], FALLBACK[p[1]]
        else:
            sp, sub = p[1], p[2]
        seq = cons.get(sp, {}).get(sub)
        if not seq:
            sys.stderr.write("no consensus for %s (%s/%s)\n" % (base, sp, sub))
            missing += 1
            continue
        # Trim the stored 300 bp flanks back to KEEP bp. 300 bp of unrelated
        # flanking DNA per side forces MAFFT to invent thousands of junk
        # columns: measured, a ~250 bp element ends up smeared over 1200-1475
        # alignment columns, and the dilution grows with copy number. The
        # background identity does not need to come from the alignment - it is
        # measured separately from random genomic windows.
        cut = 300 - KEEP
        with open(os.path.join(dst, base + ".fa"), "w") as out:
            out.write(">CONSENSUS_%s\n%s\n" % (sub, seq))
            nm, buf = None, []
            for line in open(f):
                line = line.rstrip()
                if line.startswith(">"):
                    if nm:
                        s = "".join(buf)
                        out.write(">%s\n%s\n" % (nm, s[cut:len(s) - cut]))
                    nm, buf = line[1:], []
                else:
                    buf.append(line)
            if nm:
                s = "".join(buf)
                out.write(">%s\n%s\n" % (nm, s[cut:len(s) - cut]))
    sys.stderr.write("wrote %d sets, %d missing\n"
                     % (len(glob.glob(os.path.join(dst, "*.fa"))), missing))


if __name__ == "__main__":
    main()
