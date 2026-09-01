#!/usr/bin/env python3
"""candidate consensus + genome  ->  ~100-locus alignments with flanks.

This is his steps 1-4, which have never existed as one command:

  1. start from a preliminary candidate consensus
  2. search the genome with it
  3. inspect the hits
  4. as several ~100-locus alignments with flanks - "random samples, best top
     hits, everything which works"

Produces, per candidate, the same three views the Tal data has:
  top100    the 100 highest-scoring hits
  rand100   a random 100
  all       everything, when there are fewer than 100

The consensus is written first in every alignment, because that is what the
scoring code anchors the element on.
"""
import os
import random
import subprocess
import sys

FLANK = 400
MAXCOPIES = 100
MIN_HITS = 5


def sh(c):
    return subprocess.run(c, shell=True, capture_output=True, text=True)


def read_fa(p):
    names, seqs, cur, buf = [], [], None, []
    for l in open(p):
        l = l.rstrip("\n\r")
        if l.startswith(">"):
            if cur is not None:
                seqs.append("".join(buf))
            cur = l[1:]
            names.append(cur)
            buf = []
        else:
            buf.append(l.strip())
    if cur is not None:
        seqs.append("".join(buf))
    return names, seqs


def one(cons_name, cons_seq, genome, out, tag):
    os.makedirs(out, exist_ok=True)
    q = os.path.join(out, "%s.q.fa" % tag)
    with open(q, "w") as fh:
        fh.write(">%s\n%s\n" % (cons_name, cons_seq))

    # 2. search the genome
    b6 = os.path.join(out, "%s.hits.tsv" % tag)
    r = sh("blastn -query %s -subject %s -evalue 1e-10 -word_size 11 -dust no "
           "-outfmt '6 sseqid sstart send pident length bitscore sstrand' "
           "-max_target_seqs 100000 > %s" % (q, genome, b6))
    if r.returncode or not os.path.exists(b6):
        return None, "blastn failed: %s" % r.stderr[:120]

    hits = []
    for l in open(b6):
        f = l.rstrip("\n").split("\t")
        if len(f) < 7:
            continue
        c, s, e = f[0], int(f[1]), int(f[2])
        strand = "-" if s > e else "+"
        lo, hi = (e, s) if s > e else (s, e)
        if hi - lo < 40:
            continue
        hits.append((c, lo - 1, hi, float(f[5]), strand))
    if len(hits) < MIN_HITS:
        return None, "only %d hits" % len(hits)

    # dedupe overlapping hits, keep the best-scoring
    hits.sort(key=lambda h: (h[0], h[1]))
    keep, last = [], None
    for h in hits:
        if last and h[0] == last[0] and h[1] < last[2]:
            if h[3] > last[3]:
                keep[-1] = h
                last = h
            continue
        keep.append(h)
        last = h
    hits = keep

    gsize = os.path.join(out, "genome.sizes")
    if not os.path.exists(gsize):
        sh("cut -f1,2 %s.fai > %s" % (genome, gsize))

    views = {}
    byscore = sorted(hits, key=lambda h: -h[3])
    views["top100"] = byscore[:MAXCOPIES]
    if len(hits) > MAXCOPIES:
        rng = random.Random(0)
        views["rand100"] = rng.sample(hits, MAXCOPIES)
    else:
        views["all"] = hits

    made = []
    for view, sel in views.items():
        bed = os.path.join(out, "%s.%s.bed" % (tag, view))
        with open(bed, "w") as fh:
            for i, (c, s, e, sc, st) in enumerate(sel):
                fh.write("%s\t%d\t%d\t%s_%d\t0\t%s\n" % (c, s, e, view, i, st))
        sl = bed + ".slop"
        sh("bedtools slop -i %s -g %s -b %d > %s" % (bed, gsize, FLANK, sl))
        fa = os.path.join(out, "%s.%s.fa" % (tag, view))
        sh("bedtools getfasta -fi %s -bed %s -s -name+ -fo %s" % (genome, sl, fa))
        if not os.path.exists(fa):
            continue
        fn, fs = read_fa(fa)
        inp = os.path.join(out, "%s.%s.in.fa" % (tag, view))
        with open(inp, "w") as fh:
            fh.write(">CONSENSUS_%s\n%s\n" % (tag, cons_seq))
            for n, s in zip(fn, fs):
                fh.write(">%s\n%s\n" % (n, s))
        aln = os.path.join(out, "%s__%s.aln.fa" % (tag, view))
        sh("mafft --retree 2 --maxiterate 0 --adjustdirection --quiet --thread 8 %s > %s"
           % (inp, aln))
        if os.path.exists(aln) and os.path.getsize(aln) > 0:
            made.append((view, len(sel), aln))
    return made, None


def main():
    cons_fa, genome, out, prefix = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    names, seqs = read_fa(cons_fa)
    print("candidates in %s: %d" % (cons_fa, len(names)))
    for n, s in zip(names, seqs):
        cid = n.split()[0]
        tag = "%s_%s" % (prefix, cid)
        made, err = one(cid, s, genome, out, tag)
        if err:
            print("  %-18s %s" % (tag, err))
            continue
        for view, n_hits, path in made:
            print("  %-18s %-8s %4d copies -> %s" % (tag, view, n_hits, os.path.basename(path)))


if __name__ == "__main__":
    main()
