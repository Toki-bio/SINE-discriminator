#!/usr/bin/env python3
"""Six different things "mosaic" could mean, so Sergei can pick the closest.

Same real Tal consensus (saq s3), same low divergence (0.08) so the structure is
visible rather than buried in noise, same flanks, same TSDs. ONLY the
rearrangement differs. Segments are sixths of the element, ~42 bp each, labelled
a-f below.

  A  DUP        one segment replaces another:      a b c d e f -> a b b d e f
  B  SWAP       two segments exchange places:      a b c d e f -> a e c d b f
  C  DEL        a segment is missing:              a b c d e f -> a b _ d e f
  D  INS        foreign DNA inserted mid-element:  a b c d e f -> a b X c d e f
  E  TANDEM     a segment duplicated in place:     a b c d e f -> a b c c d e f
  F  KALEID     every copy rearranged differently - the "different copies carry
                different blocks" case, mixing all of the above per copy

In A-E the SAME rearrangement is applied to every affected copy, so the set has
one consistent variant plus wild type. In F no two copies agree.
"""
import os, sys, random
OUT = "/data/W/toki/SINE_disc"
FL = 70
BASES = "ACGT"
rng = random.Random(2026)
NSEG = 6
AGE = 0.08
FRAC = 0.5                     # half the copies carry the rearrangement


def read(p):
    d, n = {}, None
    for line in open(p):
        line = line.rstrip()
        if line.startswith(">"):
            n = line[1:]
            d[n] = []
        elif n:
            d[n].append(line)
    return dict((k, "".join(v)) for k, v in d.items())


def mutate(s, r, indel=0.06):
    out = []
    for c in s:
        x = rng.random()
        if x < r * indel:
            continue
        if x < r * indel * 2:
            out.append(rng.choice(BASES)); out.append(c)
        elif x < r:
            out.append(rng.choice([b for b in BASES if b != c]))
        else:
            out.append(c)
    return "".join(out)


def main():
    src = read("%s/sets_c/POS__saq__s3_43seqs.fa" % OUT)
    cons = [v for k, v in src.items() if k.startswith("CONSENSUS_")][0]
    rnd = read("%s/sets_c/NEGRAND__saq__r00.fa" % OUT)
    bg = [v for k, v in rnd.items() if not k.startswith("CONSENSUS_")]
    L = len(cons)
    b = [round(i * L / NSEG) for i in range(NSEG + 1)]
    seg = [cons[b[i]:b[i + 1]] for i in range(NSEG)]

    def flank(n):
        d = rng.choice(bg)
        o = rng.randint(0, max(0, len(d) - n - 1))
        return d[o:o + n]

    def foreign(n):
        d = rng.choice(bg)
        o = rng.randint(0, max(0, len(d) - n - 1))
        return d[o:o + n]

    def build(kind):
        s = list(seg)
        if kind == "DUP":
            s[2] = s[1]                                   # c <- b
        elif kind == "SWAP":
            s[1], s[4] = s[4], s[1]                       # b <-> e
        elif kind == "DEL":
            s = s[:2] + s[3:]                             # drop c
        elif kind == "INS":
            s = s[:2] + [foreign(45)] + s[2:]             # foreign block after b
        elif kind == "TANDEM":
            s = s[:3] + [s[2]] + s[3:]                    # c c
        elif kind == "KALEID":
            s = list(seg)
            op = rng.choice(["dup", "swap", "del", "ins", "tandem"])
            i, j = rng.sample(range(len(s)), 2)
            if op == "dup":
                s[j] = s[i]
            elif op == "swap":
                s[i], s[j] = s[j], s[i]
            elif op == "del":
                s.pop(i)
            elif op == "ins":
                s.insert(i, foreign(rng.randint(30, 60)))
            else:
                s.insert(i, s[i])
        return "".join(s)

    for kind in ["DUP", "SWAP", "DEL", "INS", "TANDEM", "KALEID"]:
        recs = []
        for i in range(100):
            body = build(kind) if (kind == "KALEID" or rng.random() < FRAC) else cons
            t = "".join(rng.choice(BASES) for _ in range(12))
            recs.append(flank(58) + t + mutate(body, AGE) + "A" * 14 + t + flank(58))
        name = "MOSAIC%s__saq__s3" % kind
        with open("%s/sets_c/%s.fa" % (OUT, name), "w") as fh:
            fh.write(">CONSENSUS_s3\n%s\n" % cons)
            for i, r in enumerate(recs):
                fh.write(">m%03d\n%s\n" % (i, r))
        print("wrote", name)


if __name__ == "__main__":
    main()
