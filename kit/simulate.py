#!/usr/bin/env python3
"""Synthetic families with known ground truth.

Every negative so far has been carved out of real data, which means its truth is
only as good as the carving - NEGSPLICE turned out not to be a negative at all.
A generator gives exact control and exact labels, so a statistic can be tested
against the thing it claims to measure.

The generative process is the one the spec describes: plant a consensus in real
genomic background, then vary age, truncation, copy number, TSD, poly-A, and the
structure of interest. Background is real genomic sequence from the same species,
so composition and repeat content are not idealised away.

Grids produced (label -> what it should be):
  SIMCLEAN     one family, varying age                       -> SINE
  SIMMOSAIC    fraction f of copies are A/B recombinants,
               breakpoint spread s                           -> mosaic
  SIMSUBFAM    two subfamilies at divergence d, no recomb.   -> subfamily structure
  SIMNEST      fraction f of copies share a host flank       -> nested
  SIMTRUNC     geometric 5' truncation at rate r             -> truncated
"""
import os, random, sys
OUT = "/data/W/toki/SINE_disc"
FL = 70
rng = random.Random(20260831)
BASES = "ACGT"


def read(p):
    d, n = {}, None
    for line in open(p):
        line = line.rstrip()
        if line.startswith(">"):
            n = line[1:]
            d[n] = []
        elif n:
            d[n].append(line)
    return {k: "".join(v) for k, v in d.items()}


def mutate(seq, rate, indel=0.08):
    out = []
    for c in seq:
        r = rng.random()
        if r < rate * indel:
            continue                                    # deletion
        if r < rate * indel * 2:
            out.append(rng.choice(BASES))               # insertion
            out.append(c)
        elif r < rate:
            out.append(rng.choice([b for b in BASES if b != c]))
        else:
            out.append(c)
    return "".join(out)


def main():
    sp = "saq"
    src = read("%s/sets_c/POS__%s__s3_43seqs.fa" % (OUT, sp))
    consA = [v for k, v in src.items() if k.startswith("CONSENSUS_")][0]
    src2 = read("%s/sets_c/POS__%s__s8_225seqs.fa" % (OUT, sp))
    consB = [v for k, v in src2.items() if k.startswith("CONSENSUS_")][0]
    rand = read("%s/sets_c/NEGRAND__%s__r00.fa" % (OUT, sp))
    bg = [v for k, v in rand.items() if not k.startswith("CONSENSUS_")]
    if not bg:
        sys.exit("no background")

    def flank(n):
        d = rng.choice(bg)
        o = rng.randint(0, max(0, len(d) - n - 1))
        return d[o:o + n]

    def one_copy(body, age, tsd_len=12, polya=14, trunc=0):
        b = mutate(body[trunc:], age)
        t = "".join(rng.choice(BASES) for _ in range(tsd_len)) if tsd_len else ""
        return flank(FL - len(t)) + t + b + "A" * rng.randint(polya - 4, polya + 4) + \
               t + flank(FL - len(t))

    def emit(name, recs, cons):
        with open("%s/sets_c/%s.fa" % (OUT, name), "w") as fh:
            fh.write(">CONSENSUS_sim\n%s\n" % cons)
            for i, r in enumerate(recs):
                fh.write(">sim%03d\n%s\n" % (i, r))

    made = []
    # ---- clean families across the age range -----------------------------
    for age in (0.05, 0.12, 0.20, 0.30):
        recs = [one_copy(consA, age) for _ in range(100)]
        n = "SIMCLEAN__age%03d" % int(age * 100)
        emit(n, recs, consA)
        made.append(n)

    # ---- mosaics: fraction f recombinant, breakpoint spread s ------------
    for f in (0.2, 0.5, 1.0):
        for spread in ("tight", "wide"):
            recs = []
            for i in range(100):
                if rng.random() < f:
                    bp = rng.uniform(0.45, 0.55) if spread == "tight" else rng.uniform(0.15, 0.85)
                    a = consA[:int(len(consA) * bp)]
                    b = consB[int(len(consB) * bp):]
                    recs.append(one_copy(a + b, 0.12))
                else:
                    recs.append(one_copy(consA, 0.12))
            n = "SIMMOSAIC__f%03d__%s" % (int(f * 100), spread)
            emit(n, recs, consA)
            made.append(n)

    # ---- two subfamilies at controlled divergence ------------------------
    for d in (0.05, 0.12, 0.25):
        sub = mutate(consA, d)
        recs = [one_copy(consA if i < 50 else sub, 0.10) for i in range(100)]
        n = "SIMSUBFAM__d%03d" % int(d * 100)
        emit(n, recs, consA)
        made.append(n)

    # ---- nested copies: a shared host flank on a fraction ------------------
    host = "".join(rng.choice(BASES) for _ in range(400))
    for f in (0.2, 0.5):
        recs = []
        for i in range(100):
            c = one_copy(consA, 0.12)
            if rng.random() < f:
                o = rng.randint(0, 200)
                c = mutate(host[o:o + FL], 0.05) + c[FL:len(c) - FL] + \
                    mutate(host[o + 200:o + 200 + FL], 0.05)
            recs.append(c)
        n = "SIMNEST__f%03d" % int(f * 100)
        emit(n, recs, consA)
        made.append(n)

    # ---- 5' truncation ----------------------------------------------------
    for r in (0.3, 0.6):
        recs = []
        for i in range(100):
            t = int(rng.expovariate(1.0 / (r * len(consA)))) if rng.random() < 0.8 else 0
            recs.append(one_copy(consA, 0.12, trunc=min(t, len(consA) - 60)))
        n = "SIMTRUNC__r%03d" % int(r * 100)
        emit(n, recs, consA)
        made.append(n)

    print("generated %d simulated sets" % len(made))
    for n in made:
        print("  " + n)


if __name__ == "__main__":
    main()


def scrambled():
    """The mosaic Sergei actually means: abcd -> abbd, adcd. A copy carries one
    of the element's OWN segments in a slot where the consensus has a different
    one. Not a two-parent recombinant - an internal rearrangement, and which
    slot is affected varies between copies."""
    sp = "saq"
    src = read("%s/sets_c/POS__%s__s3_43seqs.fa" % (OUT, sp))
    consA = [v for k, v in src.items() if k.startswith("CONSENSUS_")][0]
    rand = read("%s/sets_c/NEGRAND__%s__r00.fa" % (OUT, sp))
    bg = [v for k, v in rand.items() if not k.startswith("CONSENSUS_")]

    def flank(n):
        d = rng.choice(bg)
        o = rng.randint(0, max(0, len(d) - n - 1))
        return d[o:o + n]

    NSEG = 6
    L = len(consA)
    bounds = [round(i * L / NSEG) for i in range(NSEG + 1)]
    segs = [consA[bounds[i]:bounds[i + 1]] for i in range(NSEG)]
    made = []
    for f in (0.2, 0.5, 1.0):
        for mode in ("swap", "dup"):
            recs = []
            for _ in range(100):
                s = list(segs)
                if rng.random() < f:
                    i, j = rng.sample(range(NSEG), 2)
                    if mode == "dup":
                        s[j] = s[i]                 # segment i also occupies slot j
                    else:
                        s[i], s[j] = s[j], s[i]     # two slots exchange contents
                body = "".join(s)
                b = mutate(body, 0.10)
                t = "".join(rng.choice(BASES) for _ in range(12))
                recs.append(flank(58) + t + b + "A" * 14 + t + flank(58))
            n = "SIMSCRAM__%s__f%03d" % (mode, int(f * 100))
            with open("%s/sets_c/%s.fa" % (OUT, n), "w") as fh:
                fh.write(">CONSENSUS_sim\n%s\n" % consA)
                for i, r in enumerate(recs):
                    fh.write(">scr%03d\n%s\n" % (i, r))
            made.append(n)
    print("generated %d scrambled sets" % len(made))


if __name__ == "__main__" and "--scram" in sys.argv:
    scrambled()


def deleted():
    """Sergei: in non-SINE repeats the mosaic often includes a long stretch
    MISSING in some copies where others have sequence. Unlike a swap, a deletion
    survives alignment as a shared gap block, so it should be detectable."""
    sp = "saq"
    src = read("%s/sets_c/POS__%s__s3_43seqs.fa" % (OUT, sp))
    consA = [v for k, v in src.items() if k.startswith("CONSENSUS_")][0]
    rand = read("%s/sets_c/NEGRAND__%s__r00.fa" % (OUT, sp))
    bg = [v for k, v in rand.items() if not k.startswith("CONSENSUS_")]

    def flank(n):
        d = rng.choice(bg)
        o = rng.randint(0, max(0, len(d) - n - 1))
        return d[o:o + n]

    L = len(consA)
    for f in (0.2, 0.5):
        for mode in ("one", "many"):
            recs = []
            for _ in range(100):
                body = consA
                if rng.random() < f:
                    if mode == "one":
                        st = round(L * 0.45)
                        body = consA[:st] + consA[st + round(L * 0.2):]
                    else:                       # a different block per copy
                        st = rng.randint(20, L - 70)
                        ln = rng.randint(30, 70)
                        body = consA[:st] + consA[st + ln:]
                b = mutate(body, 0.10)
                t = "".join(rng.choice(BASES) for _ in range(12))
                recs.append(flank(58) + t + b + "A" * 14 + t + flank(58))
            n = "SIMDEL__%s__f%03d" % (mode, int(f * 100))
            with open("%s/sets_c/%s.fa" % (OUT, n), "w") as fh:
                fh.write(">CONSENSUS_sim\n%s\n" % consA)
                for i, r in enumerate(recs):
                    fh.write(">del%03d\n%s\n" % (i, r))
    print("generated 4 deletion sets")


if __name__ == "__main__" and "--del" in sys.argv:
    deleted()
