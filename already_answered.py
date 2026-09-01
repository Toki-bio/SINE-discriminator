#!/usr/bin/env python3
"""Refuse to ask about an alignment whose answer he has already given.

He judged eight NEGCHIM sets - five of them "consensus too long but good sine" -
and I put a ninth NEGCHIM set in the next review batch. He answered it with the
same words and said, rightly, that it "perfectly matches what i have already
answered".

So: before asking about an alignment, check whether he has already answered
enough alignments of the same construction, consistently. If he has, predict the
answer from his previous ones and do not ask.

Prediction is only offered where his previous answers AGREE. Where they differ,
the alignment is genuinely worth asking about.
"""
import io
import json
import os
import sys
from collections import Counter, defaultdict

CALLS = "calls.tsv"
MIN_AGREE = 3        # this many previous answers of the same class...
AGREE_FRAC = 0.6     # ...and at least this share saying the same thing


def coarse(call):
    """The substance of an answer, ignoring which extra property it names.

    His nine NEGCHIM answers were SINE, GREY, SINE_CONSENSUS_TOO_LONG x4,
    SINE_CONSENSUS_TOO_LONG_NOT_MIXTURE, NOT_MIXTURE x2. Grouped by exact label
    that looks like disagreement (top answer only 44%), but eight of the nine
    say the same thing: it IS a SINE. Only the extra property differs.
    """
    c = call.upper()
    if c.startswith("NOT_SINE") or c in ("UNUSABLE", "UNUSABLE_UNALIGNED"):
        return "not a SINE"
    if c.startswith("SINE") or c in ("NOT_MIXTURE", "MOSAIC", "NO_FLANKS"):
        return "a SINE"
    if "MIXTURE" in c:
        return "a mixture"
    if "GREY" in c or "BADLY" in c or "ANCIENT" in c:
        return "unclear or badly presented"
    return c.lower()


def klass(name):
    """The construction a set comes from - the part before the first '__'."""
    return name.split("__")[0]


def load_calls(path=CALLS):
    rows = io.open(path, encoding="utf-8").read().rstrip("\n").split("\n")[1:]
    out = []
    for r in rows:
        f = r.split("\t")
        if len(f) >= 3:
            out.append((f[0], f[1], f[2]))
    return out


def build(calls):
    by = defaultdict(list)
    for name, corp, call in calls:
        by[klass(name)].append(coarse(call))
    known = {}
    for k, v in by.items():
        if len(v) < MIN_AGREE:
            continue
        c = Counter(v).most_common()
        top, n = c[0]
        if n / float(len(v)) >= AGREE_FRAC:
            known[k] = {"answer": top, "n_seen": len(v), "agreement": round(n / float(len(v)), 2),
                        "others": [a for a, _ in c[1:]]}
    return known


def check(names, calls=None):
    calls = calls or load_calls()
    judged = {n for n, _, _ in calls}
    known = build(calls)
    ask, skip = [], []
    for n in names:
        if n in judged:
            skip.append((n, "already judged directly"))
            continue
        k = known.get(klass(n))
        if k:
            skip.append((n, "%s: he answered %d of these, %d%% saying '%s'"
                         % (klass(n), k["n_seen"], int(k["agreement"] * 100), k["answer"])))
        else:
            ask.append(n)
    return ask, skip, known


def main():
    calls = load_calls()
    known = build(calls)
    print("Constructions he has answered enough times to predict (%d of %d judged):"
          % (len(known), len({klass(n) for n, _, _ in calls})))
    for k in sorted(known, key=lambda k: -known[k]["n_seen"]):
        v = known[k]
        print("  %-14s n=%-3d %3d%% agree -> %-34s %s"
              % (k, v["n_seen"], int(v["agreement"] * 100), v["answer"],
                 ("also: " + ", ".join(v["others"])) if v["others"] else ""))

    if len(sys.argv) > 1:
        names = json.load(open(sys.argv[1]))
        names = [x[1] if isinstance(x, list) else x for x in names]
        ask, skip, _ = check(names, calls)
        print()
        print("of %d proposed, DO NOT ASK %d:" % (len(names), len(skip)))
        for n, why in skip:
            print("   %-38s %s" % (n, why))
        print("still worth asking (%d):" % len(ask))
        for n in ask:
            print("   %s" % n)


if __name__ == "__main__":
    main()
