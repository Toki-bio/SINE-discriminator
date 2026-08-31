#!/usr/bin/env python3
"""One gradient per feature: 100 alignments each, parameter swept continuously.

Every statistic in this project was calibrated on a handful of hand-carved sets,
which is how NEGSPLICE spent three builds being scored as a negative when it was
not one. A continuous sweep with exact labels answers a sharper question than
"does it separate": does the statistic move MONOTONICALLY with the thing it
claims to measure, and does it stay flat for everything else?

Each grid varies ONE parameter across 100 sets; everything else is held at a
realistic default. Background and the seed consensus are real sequence from saq,
so composition and repeat content are not idealised away.

Usage:  python3.12 gradients.py [outdir]
"""
import os, sys, random, json

OUT = "/data/W/toki/SINE_disc"
FL = 70
BASES = "ACGT"
rng = random.Random(4242)

DEFAULTS = dict(age=0.12, n=100, trunc=0.0, contam=0.0, subdiv=0.0, subfrac=0.0,
                scram=0.0, dup=0.0, delshared=0.0, delcopy=0.0, nest=0.0,
                tsd=12, polya=14, recomb=0.0)

# name -> (parameter, low, high)   100 steps each
GRIDS = [
    ("AGE",       "age",       0.02, 0.35),
    ("NCOPIES",   "n",         20,   200),
    ("TRUNC",     "trunc",     0.0,  0.70),
    ("CONTAM",    "contam",    0.0,  0.60),
    ("SUBDIV",    "subdiv",    0.0,  0.30),
    ("SUBFRAC",   "subfrac",   0.0,  0.50),
    ("SCRAM",     "scram",     0.0,  1.00),
    ("DUP",       "dup",       0.0,  1.00),
    ("DELSHARED", "delshared", 0.0,  0.80),
    ("DELCOPY",   "delcopy",   0.0,  0.80),
    ("NEST",      "nest",      0.0,  0.80),
    ("TSD",       "tsd",       0,    22),
    ("POLYA",     "polya",     0,    40),
    ("RECOMB",    "recomb",    0.0,  1.00),
]
NSTEP = 100
NSEG = 6


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
            continue
        if r < rate * indel * 2:
            out.append(rng.choice(BASES))
            out.append(c)
        elif r < rate:
            out.append(rng.choice([b for b in BASES if b != c]))
        else:
            out.append(c)
    return "".join(out)


def build(p, consA, consB, bg, host):
    """One set of copies under parameter dict p."""
    L = len(consA)
    b = [round(i * L / NSEG) for i in range(NSEG + 1)]
    segs = [consA[b[i]:b[i + 1]] for i in range(NSEG)]
    subcons = mutate(consA, p["subdiv"]) if p["subdiv"] > 0 else consA

    def flank(n):
        d = rng.choice(bg)
        o = rng.randint(0, max(0, len(d) - n - 1))
        return d[o:o + n]

    recs = []
    for i in range(int(p["n"])):
        if p["contam"] and rng.random() < p["contam"]:
            d = rng.choice(bg)
            o = rng.randint(0, max(0, len(d) - (L + 2 * FL) - 1))
            recs.append(d[o:o + L + 2 * FL])
            continue

        base = subcons if (p["subfrac"] and rng.random() < p["subfrac"]) else consA
        s = list(segs) if base is consA else None

        if s is not None and (p["scram"] or p["dup"]):
            if p["scram"] and rng.random() < p["scram"]:
                i1, j1 = rng.sample(range(NSEG), 2)
                s[i1], s[j1] = s[j1], s[i1]
            if p["dup"] and rng.random() < p["dup"]:
                i1, j1 = rng.sample(range(NSEG), 2)
                s[j1] = s[i1]
            base = "".join(s)

        if p["recomb"] and rng.random() < p["recomb"]:
            bp = rng.uniform(0.2, 0.8)
            base = consA[:int(L * bp)] + consB[int(len(consB) * bp):]

        if p["delshared"] and rng.random() < p["delshared"]:
            st, ln = round(L * 0.45), round(L * 0.2)
            base = base[:st] + base[st + ln:]
        if p["delcopy"] and rng.random() < p["delcopy"]:
            st = rng.randint(20, max(21, len(base) - 70))
            base = base[:st] + base[st + rng.randint(30, 70):]

        tr = 0
        if p["trunc"] and rng.random() < 0.8:
            tr = min(int(rng.expovariate(1.0 / max(1e-6, p["trunc"] * L))), len(base) - 60)
        body = mutate(base[max(0, tr):], p["age"])

        tsd = "".join(rng.choice(BASES) for _ in range(int(p["tsd"])))
        tail = "A" * max(0, int(p["polya"]) + rng.randint(-3, 3))
        lf, rf = flank(FL - len(tsd)), flank(FL - len(tsd))
        if p["nest"] and rng.random() < p["nest"]:
            o = rng.randint(0, 200)
            lf = mutate(host[o:o + FL - len(tsd)], 0.05)
            rf = mutate(host[o + 220:o + 220 + FL - len(tsd)], 0.05)
        recs.append(lf + tsd + body + tail + tsd + rf)
    return recs


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else OUT + "/grad_sets"
    os.makedirs(outdir, exist_ok=True)
    src = read("%s/sets_c/POS__saq__s3_43seqs.fa" % OUT)
    consA = [v for k, v in src.items() if k.startswith("CONSENSUS_")][0]
    src2 = read("%s/sets_c/POS__saq__s8_225seqs.fa" % OUT)
    consB = [v for k, v in src2.items() if k.startswith("CONSENSUS_")][0]
    rand = read("%s/sets_c/NEGRAND__saq__r00.fa" % OUT)
    bg = [v for k, v in rand.items() if not k.startswith("CONSENSUS_")]
    host = "".join(rng.choice(BASES) for _ in range(600))

    truth = {}
    for name, par, lo, hi in GRIDS:
        for k in range(NSTEP):
            val = lo + (hi - lo) * k / (NSTEP - 1)
            if par == "n" or par in ("tsd", "polya"):
                val = int(round(val))
            p = dict(DEFAULTS)
            p[par] = val
            if par == "subfrac":
                p["subdiv"] = 0.15          # a fraction is meaningless at zero divergence
            recs = build(p, consA, consB, bg, host)
            sid = "GRAD_%s_%03d" % (name, k)
            with open(os.path.join(outdir, sid + ".fa"), "w") as fh:
                fh.write(">CONSENSUS_sim\n%s\n" % consA)
                for i, r in enumerate(recs):
                    fh.write(">c%03d\n%s\n" % (i, r))
            truth[sid] = {"grid": name, "param": par, "value": val, "step": k}
    json.dump(truth, open(os.path.join(outdir, "truth.json"), "w"), indent=0)
    print("wrote %d sets across %d grids" % (len(truth), len(GRIDS)))


if __name__ == "__main__":
    main()
