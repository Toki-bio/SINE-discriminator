#!/usr/bin/env python3
"""Factorial pairs: where two perturbations are mistaken for a third.

Single gradients showed each statistic in isolation. The dangerous cases are
interactions - most of all age x subfamily divergence, because sub_gap responds
to age at rho +0.99, so an OLD SINGLE family may be indistinguishable from TWO
YOUNG subfamilies. Each pair is a 10 x 10 grid.
"""
import os, sys, json
sys.path.insert(0, "/data/W/toki")
import gradients as G

PAIRS = [("AGE", "age", 0.02, 0.32, "SUBDIV", "subdiv", 0.0, 0.28),
         ("AGE", "age", 0.02, 0.32, "CONTAM", "contam", 0.0, 0.55),
         ("AGE", "age", 0.02, 0.32, "TRUNC", "trunc", 0.0, 0.60),
         ("AGE", "age", 0.02, 0.32, "POLYA", "polya", 0, 38),
         ("NCOPIES", "n", 25, 190, "SUBDIV", "subdiv", 0.0, 0.28),
         ("CONTAM", "contam", 0.0, 0.55, "SCRAM", "scram", 0.0, 1.0)]
N = 10


def main():
    out = G.OUT + "/combo_sets"
    os.makedirs(out, exist_ok=True)
    src = G.read("%s/sets_c/POS__saq__s3_43seqs.fa" % G.OUT)
    consA = [v for k, v in src.items() if k.startswith("CONSENSUS_")][0]
    src2 = G.read("%s/sets_c/POS__saq__s8_225seqs.fa" % G.OUT)
    consB = [v for k, v in src2.items() if k.startswith("CONSENSUS_")][0]
    rand = G.read("%s/sets_c/NEGRAND__saq__r00.fa" % G.OUT)
    bg = [v for k, v in rand.items() if not k.startswith("CONSENSUS_")]
    host = "".join(G.rng.choice(G.BASES) for _ in range(600))

    truth = {}
    for n1, p1, lo1, hi1, n2, p2, lo2, hi2 in PAIRS:
        for i in range(N):
            for j in range(N):
                v1 = lo1 + (hi1 - lo1) * i / (N - 1)
                v2 = lo2 + (hi2 - lo2) * j / (N - 1)
                if p1 in ("n", "tsd", "polya"):
                    v1 = int(round(v1))
                if p2 in ("n", "tsd", "polya"):
                    v2 = int(round(v2))
                p = dict(G.DEFAULTS)
                p[p1], p[p2] = v1, v2
                if "subdiv" in (p1, p2):
                    p["subfrac"] = 0.5
                recs = G.build(p, consA, consB, bg, host)
                sid = "COMBO_%s_%s_%02d_%02d" % (n1, n2, i, j)
                with open(os.path.join(out, sid + ".fa"), "w") as fh:
                    fh.write(">CONSENSUS_sim\n%s\n" % consA)
                    for kk, r in enumerate(recs):
                        fh.write(">c%03d\n%s\n" % (kk, r))
                truth[sid] = {"pair": n1 + "x" + n2, "p1": n1, "v1": v1,
                              "p2": n2, "v2": v2}
    json.dump(truth, open(os.path.join(out, "truth.json"), "w"), indent=0)
    print("wrote %d combo sets" % len(truth))


if __name__ == "__main__":
    main()
