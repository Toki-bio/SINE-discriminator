#!/usr/bin/env python3
"""Microsatellite content - the criterion he asked for by name.

On hyd_SINE_7: "looks like combination of microsatellites to me, not sine - not
caught by your filters! - needs new criteria on microsatellite content?"
On hyd_SINE_17: "looks like its not a sine cause surrounded by long 2-nt
microsatellites".

Two different places, so two different numbers:

  msat_elem    how much of the ELEMENT is inside a tandem repeat. A SINE built
               out of microsatellite is not a SINE (hyd_SINE_7).
  msat_flank   how much of the FLANK is. Copies landing in microsatellite tracts
               are a mapping artefact, not independent insertions (hyd_SINE_17).

Measured per copy on ungapped sequence, then the median over copies, so one
repeat-rich copy cannot carry the set.

Runs on the raw long-flank alignments, like the island scan: the trimmed
display alignment cuts away most of the flank.
"""
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np

GAP = "-"
MAX_PERIOD = 6
MIN_UNITS = 4          # a run must repeat the unit at least this many times
MIN_RUN_BP = 12        # ... and cover at least this many bases


def read_fa(p):
    names, seqs, cur, buf = [], [], None, []
    for line in open(p):
        line = line.rstrip("\n\r")
        if line.startswith(">"):
            if cur is not None:
                seqs.append("".join(buf))
            cur = line[1:]
            names.append(cur)
            buf = []
        else:
            buf.append(line.strip())
    if cur is not None:
        seqs.append("".join(buf))
    return names, seqs


def msat_fraction(seq):
    """Fraction of bases inside a tandem repeat of period 1..6.

    Straightforward scan rather than a library: for each period, mark every
    position whose base equals the one a period back, then keep marked stretches
    long enough to be a real tract. Marks from different periods are unioned, so
    a (CA)n inside an (A)n region is not counted twice.
    """
    s = seq.replace(GAP, "").upper()
    n = len(s)
    if n < MIN_RUN_BP:
        return 0.0, 0
    covered = np.zeros(n, dtype=bool)
    arr = np.frombuffer(s.encode(), dtype=np.uint8)
    for p in range(1, MAX_PERIOD + 1):
        if n <= p:
            break
        same = arr[p:] == arr[:-p]          # position i+p continues the period
        # find runs of True
        idx = np.flatnonzero(np.diff(np.concatenate(([0], same.view(np.int8), [0]))))
        for a, b in zip(idx[::2], idx[1::2]):
            length = b - a + p              # the tract includes the first unit
            if length >= max(MIN_RUN_BP, MIN_UNITS * p):
                covered[a:a + length] = True
    return float(covered.sum()) / n, n


def split_rows(path):
    """Consensus, then each copy as (left flank, element, right flank)."""
    names, seqs = read_fa(path)
    ci = [i for i, x in enumerate(names) if "CONSENSUS_" in x.upper()]
    if not ci or len(seqs) < 5:
        return None
    k = ci[0]
    cons = seqs[k]
    nz = [i for i, c in enumerate(cons) if c != GAP]
    if len(nz) < 40:
        return None
    lo, hi = nz[0], nz[-1]
    rows = []
    for i, s in enumerate(seqs):
        if i == k:
            continue
        rows.append((s[:lo], s[lo:hi + 1], s[hi + 1:]))
    return cons[lo:hi + 1], rows


def measure(path):
    sp = split_rows(path)
    if sp is None:
        return None
    cons_el, rows = sp
    el, fl = [], []
    for l, e, r in rows:
        f, n = msat_fraction(e)
        if n >= 40:
            el.append(f)
        for side in (l, r):
            f, n = msat_fraction(side)
            if n >= 40:
                fl.append(f)
    if not el:
        return None
    cf, _ = msat_fraction(cons_el)
    return {"msat_elem": round(float(np.median(el)), 4),
            "msat_elem_p90": round(float(np.percentile(el, 90)), 4),
            "msat_flank": round(float(np.median(fl)), 4) if fl else None,
            "msat_flank_p90": round(float(np.percentile(fl, 90)), 4) if fl else None,
            "msat_cons": round(cf, 4),
            "n_copies": len(el)}


def _one(f):
    try:
        return os.path.basename(f).replace(".aln.fa", ""), measure(f)
    except Exception as exc:
        return os.path.basename(f).replace(".aln.fa", ""), {"error": str(exc)}


def main():
    src = sys.argv[1]
    out = sys.argv[2]
    files = sorted(glob.glob(os.path.join(src, "*.aln.fa")))
    import multiprocessing as mp
    res = {}
    with mp.Pool(max(1, mp.cpu_count() - 4)) as pool:
        for name, r in pool.imap_unordered(_one, files, chunksize=4):
            if r and "error" not in r:
                res[name] = r
    json.dump(res, open(out, "w"), separators=(",", ":"))
    print("measured %d of %d" % (len(res), len(files)))

    by = defaultdict(list)
    for k, v in res.items():
        by[k.split("__")[0]].append(v)
    print("\n%-16s %5s %10s %10s %10s %10s"
          % ("group", "n", "elem med", "elem p90", "flank med", "flank p90"))
    print("-" * 66)
    for c in sorted(by, key=lambda c: -np.median([x["msat_elem"] for x in by[c]])):
        v = by[c]
        fk = [x["msat_flank"] for x in v if x["msat_flank"] is not None]
        print("%-16s %5d %10.3f %10.3f %10s %10s"
              % (c[:16], len(v),
                 np.median([x["msat_elem"] for x in v]),
                 np.median([x["msat_elem_p90"] for x in v]),
                 "-" if not fk else "%.3f" % np.median(fk),
                 "-" if not fk else "%.3f" % np.median(
                     [x["msat_flank_p90"] for x in v if x["msat_flank_p90"] is not None])))


if __name__ == "__main__":
    main()
