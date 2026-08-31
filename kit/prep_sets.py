#!/usr/bin/env python3
"""Build the SINE-discriminator corpus: flanked, unaligned FASTA sets on KIT.

For each species/subfamily, extract a pool of copies with +-FLANK bp of real
genomic flank, then emit labelled sets (positives, synthetic negatives,
subsampling side-quest replicates) as plain FASTA ready for MAFFT.
"""
import os, sys, random, subprocess, re, json
from collections import defaultdict

BASE = "/data/W/toki/Genomes/Mammalia/Eulipotyphla"
OUT = "/data/W/toki/SINE_disc"
FLANK = 300
POOL = 2000          # copies extracted per subfamily
RANDPOOL = 6000      # random genomic loci per species
SEED = 20260831

SPECIES = {
    "saq": ("saq/run_20260425_182219", "saq/GCA_004024925.1_ScaAqu_v1_BIUU_genomic.fna"),
    "ccr": ("ccr/run_20260517_205955", "ccr/GCF_000260355.1_ConCri1.0_genomic.fna"),
    "teu": ("teu/run_20260427_130055", "teu/GCA_964194135.1_mTalEur1.hap1.1_genomic.fna"),
    "dmo": ("dmo/run_20260427_103405", "dmo/GCA_051107935.1_ASM5110793v1_genomic.fna"),
}
# families used for the subsampling side quest (spread of age and copy number)
SIDEQUEST = [("saq", "s1_30seqs"), ("saq", "s7g_172seqs"), ("saq", "s5_5seqs"),
             ("ccr", "g3_71seqs"), ("teu", "t5_31seqs")]
SQ_N = [25, 50, 100, 200]
SQ_K = 20

NAME_RE = re.compile(r"^(.+):(\d+)-(\d+)\(([+-])\)$")


def sh(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode:
        sys.stderr.write("FAILED: %s\n%s\n" % (cmd, r.stderr[-2000:]))
        raise SystemExit(1)
    return r.stdout


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return dict((k, "".join(v)) for k, v in seqs.items())


def stratified(rows, n, rng):
    """rows = [(payload, bitscore)] -> n copies spread evenly over bitscore deciles."""
    rows = sorted(rows, key=lambda r: r[1])
    if len(rows) <= n:
        return list(rows)
    out, per = [], max(1, n // 10)
    for d in range(10):
        lo, hi = len(rows) * d // 10, len(rows) * (d + 1) // 10
        chunk = rows[lo:hi]
        out += rng.sample(chunk, min(per, len(chunk)))
    rng.shuffle(out)
    return out[:n]


def write_set(sets_dir, label, records, manifest):
    path = os.path.join(sets_dir, label + ".fa")
    with open(path, "w") as fh:
        for nm, sq in records:
            fh.write(">%s\n%s\n" % (nm, sq))
    manifest.append({"set": label, "n": len(records)})


def main():
    rng = random.Random(SEED)
    sets_dir = os.path.join(OUT, "sets")
    os.makedirs(sets_dir, exist_ok=True)
    tmp = os.path.join(OUT, "tmp")
    os.makedirs(tmp, exist_ok=True)
    manifest = []

    pools = {}        # (sp, subfam) -> [(name, seq)]
    randpools = {}    # sp -> [(name, seq)]

    for sp in sorted(SPECIES):
        rundir, genome = SPECIES[sp]
        gpath = os.path.join(BASE, genome)
        fai = gpath + ".fai"
        if not os.path.exists(fai):
            sh("samtools faidx %s" % gpath)
        gsize = os.path.join(tmp, sp + ".gsize")
        sh("cut -f1,2 %s > %s" % (fai, gsize))

        af = os.path.join(BASE, rundir, "results/assignment_full.tsv")
        bysub = defaultdict(list)
        lens = []
        with open(af) as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 5 or f[4] != "assigned":
                    continue
                m = NAME_RE.match(f[0])
                if not m:
                    continue
                ctg, s, e, st = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
                bysub[f[1]].append((ctg, s, e, st, int(f[2])))
                lens.append(e - s)
        medlen = sorted(lens)[len(lens) // 2]
        sys.stderr.write("[%s] %d subfams, median copy len %d\n" % (sp, len(bysub), medlen))

        # --- real copies, one pool per subfamily -------------------------
        for sub in sorted(bysub):
            rows = bysub[sub]
            keep = stratified([(r, r[4]) for r in rows], POOL, rng)
            bed = os.path.join(tmp, "%s.%s.bed" % (sp, sub))
            with open(bed, "w") as fh:
                for r, _bs in keep:
                    ctg, s, e, st, bs = r
                    fh.write("%s\t%d\t%d\t%s|%d\t0\t%s\n" % (ctg, s, e, sub, bs, st))
            fl = bed + ".flank"
            sh("bedtools slop -i %s -g %s -b %d > %s" % (bed, gsize, FLANK, fl))
            fa = bed + ".fa"
            sh("bedtools getfasta -fi %s -bed %s -s -name+ -fo %s" % (gpath, fl, fa))
            seqs = read_fasta(fa)
            # keep only full-flank extractions (contig edges truncate some)
            pool = [(k, v) for k, v in seqs.items()
                    if len(v) >= 2 * FLANK + 60 and v.upper().count("N") < 0.1 * len(v)]
            pools[(sp, sub)] = pool
            sys.stderr.write("  %s/%s pool=%d\n" % (sp, sub, len(pool)))

        # --- random genomic loci (negatives / diluent) -------------------
        rbed = os.path.join(tmp, sp + ".rand.bed")
        sh("bedtools random -l %d -n %d -seed %d -g %s | sort -k1,1 -k2,2n > %s"
           % (medlen, RANDPOOL, SEED, gsize, rbed))
        rfl = rbed + ".flank"
        sh("bedtools slop -i %s -g %s -b %d > %s" % (rbed, gsize, FLANK, rfl))
        rfa = rbed + ".fa"
        sh("bedtools getfasta -fi %s -bed %s -s -name+ -fo %s" % (gpath, rfl, rfa))
        rs = read_fasta(rfa)
        randpools[sp] = [(k, v) for k, v in rs.items()
                         if len(v) >= 2 * FLANK + 60 and v.upper().count("N") < 0.1 * len(v)]
        sys.stderr.write("  %s randpool=%d\n" % (sp, len(randpools[sp])))

    # ===================== emit the labelled sets =========================
    fams = sorted(pools)

    # A. positives
    for key in fams:
        sp, sub = key
        rec = rng.sample(pools[key], min(100, len(pools[key])))
        write_set(sets_dir, "POS__%s__%s" % (sp, sub), rec, manifest)

    # B. side quest: K independent subsets at each n
    for sp, sub in SIDEQUEST:
        pool = pools.get((sp, sub))
        if not pool:
            continue
        for n in SQ_N:
            if len(pool) < n:
                continue
            for k in range(SQ_K):
                rec = rng.sample(pool, n)
                write_set(sets_dir, "SQ__%s__%s__n%d__k%02d" % (sp, sub, n, k), rec, manifest)

    # C. random-locus negatives
    for sp in sorted(SPECIES):
        for k in range(5):
            rec = rng.sample(randpools[sp], 100)
            write_set(sets_dir, "NEGRAND__%s__r%02d" % (sp, k), rec, manifest)

    # D. per-copy boundary jitter: element edges are individually wrong
    for key in fams:
        sp, sub = key
        rec = []
        for nm, sq in rng.sample(pools[key], min(100, len(pools[key]))):
            jl = rng.randint(20, 100) * rng.choice([-1, 1])
            jr = rng.randint(20, 100) * rng.choice([-1, 1])
            a = max(0, FLANK + jl)
            b = min(len(sq), len(sq) - FLANK + jr)
            if b - a < 80:
                continue
            rec.append((nm + "|jit", sq[max(0, a - FLANK):min(len(sq), b + FLANK)]))
        write_set(sets_dir, "NEGJITTER__%s__%s" % (sp, sub), rec, manifest)

    # E. spliced chimeras: left half of family A, right half of family B
    for k in range(20):
        ka, kb = rng.sample(fams, 2)
        spa, sa = ka
        spb, sb = kb
        A = rng.sample(pools[ka], min(50, len(pools[ka])))
        B = rng.sample(pools[kb], min(50, len(pools[kb])))
        rec = []
        for i in range(len(A)):
            na, qa = A[i]
            nb, qb = B[i % len(B)]
            ma = FLANK + (len(qa) - 2 * FLANK) // 2
            mb = FLANK + (len(qb) - 2 * FLANK) // 2
            rec.append(("chim%03d|%s|%s" % (i, sa, sb), qa[:ma] + qb[mb:]))
        write_set(sets_dir, "NEGSPLICE__%s_%s__%s_%s__%02d" % (spa, sa, spb, sb, k), rec, manifest)

    # F. 5'-truncated, LINE-fragment-like: geometric truncation, foreign 5' flank
    for key in fams:
        sp, sub = key
        rec = []
        for nm, sq in rng.sample(pools[key], min(100, len(pools[key]))):
            elen = len(sq) - 2 * FLANK
            t = int(rng.expovariate(1.0 / (0.45 * elen)))
            t = max(0, min(elen - 40, t))
            donor = rng.choice(randpools[sp])[1]
            off = rng.randint(0, max(0, len(donor) - FLANK - 1))
            rec.append((nm + "|t%d" % t, donor[off:off + FLANK] + sq[FLANK + t:]))
        write_set(sets_dir, "NEGTRUNC5__%s__%s" % (sp, sub), rec, manifest)

    # G. MIXED: real family diluted with random loci
    for key in fams:
        sp, sub = key
        for frac in (0.10, 0.30):
            nr = int(100 * frac)
            rec = rng.sample(pools[key], min(100 - nr, len(pools[key]))) \
                + [(n + "|contam", s) for n, s in rng.sample(randpools[sp], nr)]
            rng.shuffle(rec)
            write_set(sets_dir, "MIXED%02d__%s__%s" % (int(frac * 100), sp, sub), rec, manifest)

    with open(os.path.join(OUT, "manifest.json"), "w") as fh:
        json.dump({"flank": FLANK, "seed": SEED, "sets": manifest}, fh, indent=1)
    sys.stderr.write("wrote %d sets\n" % len(manifest))


if __name__ == "__main__":
    main()
