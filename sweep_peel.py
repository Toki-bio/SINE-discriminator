#!/usr/bin/env python3
"""Sweep the peel over a spread of settings — MSA-viewer's analyzeClusterability.

His reason, verbatim from script.js:

    A single run answers "did these settings cluster?", which cannot distinguish
    settings that are too strict from data with no structure to find. This runs
    the clusterer over a spread of settings and reports what changed the outcome
    and what did not, so a zero can be read for what it is.

Every subfamily I have called unseparable came from one run at one setting. t3
already turned out to be an artefact of the gap-skip, not a property of t3. This
exists so t1-4 and t6-1 are not reported the same way.

The grid mirrors his: vary each knob independently around the current setting,
then one "everything loosest" run at the panel's own minimums, because his own
comment records that the survey once passed minOccurrences 2 while the panel
accepted 1 and so "could conclude nothing clusters without ever trying the
lowest setting a user can pick".

Usage:
  sweep_peel.py <chunks.fa> <workdir> <truth.json> [focus_group ...]
"""
import json
import os
import re
import subprocess
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))


def grid():
    runs = [("current", {})]
    for v in (3, 4, 5, 8):
        runs.append(("MIN_SET %d" % v, {"MIN_SET": v}))
    for v in (2, 3, 5, 8):
        runs.append(("MIN_BLOCK %d" % v, {"MIN_BLOCK": v}))
    for v in (0.30, 0.45, 0.60, 0.75):
        runs.append(("FEAT_JACCARD %.2f" % v, {"FEAT_JACCARD": v}))
    for v in (0.40, 0.50, 0.60, 0.75):
        runs.append(("CARRY %.2f" % v, {"CARRY": v}))
    for v in (0.15, 0.25, 0.35, 0.50):
        runs.append(("EXCL_MIN %.2f" % v, {"EXCL_MIN": v}))
    for v in (0.50, 0.70, 0.90):
        runs.append(("MAX_FRAC %.2f" % v, {"MAX_FRAC": v}))
    runs.append(("no trim", {"TRIM": 0}))
    runs.append(("everything loosest", {
        "MIN_SET": 3, "MIN_BLOCK": 1, "FEAT_JACCARD": 0.20,
        "CARRY": 0.40, "EXCL_MIN": 0.05, "MAX_FRAC": 0.90, "MIN_GROUP": 3}))
    # dedupe runs whose override equals the default
    seen, out = set(), []
    for label, o in runs:
        key = tuple(sorted(o.items()))
        if key in seen:
            continue
        seen.add(key)
        out.append((label, o))
    return out


def run_one(fa, work, truth, over):
    env = dict(os.environ)
    for k, v in over.items():
        env["PEEL_" + k] = str(v)
    r = subprocess.run(
        [sys.executable, os.path.join(HERE, "peel_features.py"), fa, work, truth],
        capture_output=True, text=True, env=env, timeout=3600)
    return r.stdout


def parse(out, truth_path, focus):
    placed = purity = 0.0
    m = re.search(r"weighted mean purity: ([\d.]+) over (\d+) seqs", out)
    if m:
        purity, placed = float(m.group(1)), int(m.group(2))
    ngroups = 0
    m = re.search(r"peeled (\d+) groups", out)
    if m:
        ngroups = int(m.group(1))
    resid = 0
    m = re.search(r"residue: (\d+)", out)
    if m:
        resid = int(m.group(1))
    # best cluster per focus group, from the per-group summary block
    best = {}
    for line in out.split("\n"):
        f = line.split()
        if len(f) == 3 and f[0] in focus:
            try:
                size, pur = int(f[1]), float(f[2])
            except ValueError:
                continue
            if f[0] not in best or size > best[f[0]][0]:
                best[f[0]] = (size, pur)
    return ngroups, placed, purity, resid, best


def main():
    fa, work, truth = sys.argv[1], sys.argv[2], sys.argv[3]
    focus = sys.argv[4:] or []
    total = len(json.load(open(truth)))
    runs = grid()
    hdr = "%-22s %6s %7s %7s %7s" % ("setting", "groups", "placed", "purity", "resid")
    for g in focus:
        hdr += " %13s" % g
    print("%d chunks, %d runs\n" % (total, len(runs)))
    print(hdr)
    print("-" * len(hdr))
    for i, (label, over) in enumerate(runs):
        w = os.path.join(work, "s%02d" % i)
        try:
            out = run_one(fa, w, truth, over)
        except subprocess.TimeoutExpired:
            print("%-22s  timeout" % label)
            continue
        ng, pl, pu, rs, best = parse(out, truth, focus)
        row = "%-22s %6d %7d %7.3f %7d" % (label, ng, pl, pu, rs)
        for g in focus:
            row += " %13s" % ("%d@%.0f%%" % (best[g][0], 100 * best[g][1])
                              if g in best else "-")
        print(row, flush=True)


if __name__ == "__main__":
    main()
