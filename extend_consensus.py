#!/usr/bin/env python3
"""Act on CONSENSUS_UNDEREXTENDED instead of only reporting it.

Sergei's three starting situations include "a library needing adjustment". This
is that case, and it is the one place where the tool can do the adjustment
rather than hand it back as advice.

hyd_SINE_0 is the worked example: RepBase SINE2-2B_HM is 208 bp, the AnnoSINE
seed covers 79 of them, and the flank-decay profile says similarity carries on
150 bp past the consensus before reaching background 0.30. So ~150 bp of the
element was being treated as context - which pushed the flank background to
0.614, above the element's own identity of 0.560, and made the scorer rule that
no element existed at all.

What this does: take the raw alignment, walk outward from the element edge on
each side as far as the decay profile says the similarity lasts, and take the
majority base at every column that enough copies reach. That gives a consensus
covering the whole element. The caller then re-searches the genome with it.

The flanks in these alignments are justified - de-gapped and pushed hard against
the element - so column j outside the edge is position j from the edge in every
copy, and a majority vote down that column is meaningful.
"""
import io
import json
import os
import sys
from collections import Counter

GAP = "-"
MIN_COV = 0.40        # a column needs this fraction of copies present to vote
MIN_MAJORITY = 0.34   # ... and the winning base needs at least this share
PAD = 25              # go this much past the decay distance, decay is a 25 bp grid


def read_fa(path):
    names, seqs, cur, buf = [], [], None, []
    for line in io.open(path, encoding="utf-8", errors="replace"):
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


def majority(col, n_rows):
    bases = [c for c in col if c in "ACGTacgt"]
    if len(bases) < MIN_COV * n_rows:
        return None
    cnt = Counter(b.upper() for b in bases)
    base, k = cnt.most_common(1)[0]
    if k < MIN_MAJORITY * len(bases):
        return None
    return base


def extend(aln_path, decay, out_fa, tag):
    names, seqs = read_fa(aln_path)
    ci = [i for i, x in enumerate(names) if "CONSENSUS_" in x.upper()]
    if not ci:
        return None, "no consensus row"
    k = ci[0]
    cons = seqs[k]
    nz = [i for i, c in enumerate(cons) if c != GAP]
    if len(nz) < 40:
        return None, "consensus too short to anchor"
    lo, hi = nz[0], nz[-1]
    rows = [s for i, s in enumerate(seqs) if i != k]
    n = len(rows)

    # element part: the consensus as it stands, gaps removed
    elem = "".join(c for c in cons[lo:hi + 1] if c != GAP)

    # left: de-gap each copy's flank and index from the element edge outward
    want_L = min(400, (decay.get("decay_L") or 0) + PAD)
    want_R = min(400, (decay.get("decay_R") or 0) + PAD)
    left_rows = [s[:lo].replace(GAP, "") for s in rows]
    right_rows = [s[hi + 1:].replace(GAP, "") for s in rows]

    addL = []
    for j in range(want_L):
        col = [r[-(j + 1)] for r in left_rows if len(r) > j]
        b = majority(col, n)
        if b is None:
            break
        addL.append(b)
    addL = "".join(reversed(addL))

    addR = []
    for j in range(want_R):
        col = [r[j] for r in right_rows if len(r) > j]
        b = majority(col, n)
        if b is None:
            break
        addR.append(b)
    addR = "".join(addR)

    new = addL + elem + addR
    with io.open(out_fa, "w", encoding="utf-8") as fh:
        fh.write(">%s\n%s\n" % (tag, new))
    return {"old_len": len(elem), "added_5": len(addL), "added_3": len(addR),
            "new_len": len(new)}, None


def main():
    aln_dir, decay_json, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
    tags = sys.argv[4:]
    decay = json.load(open(decay_json))
    os.makedirs(out_dir, exist_ok=True)
    for tag in tags:
        aln = os.path.join(aln_dir, "%s__top100.aln.fa" % tag)
        if not os.path.exists(aln):
            print("  %-16s no alignment" % tag)
            continue
        d = decay.get("%s__top100" % tag, {})
        out = os.path.join(out_dir, "%s.ext.fa" % tag)
        r, err = extend(aln, d, out, "%s_ext" % tag)
        if err:
            print("  %-16s %s" % (tag, err))
            continue
        print("  %-16s %d bp -> %d bp   (+%d at 5', +%d at 3')"
              % (tag, r["old_len"], r["new_len"], r["added_5"], r["added_3"]))


if __name__ == "__main__":
    main()
