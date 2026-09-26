#!/usr/bin/env python3
"""Correct a published alignment's orientation and element edge.

Two failures this repairs:

1. MAFFT --adjustdirection reverse-complemented every copy, and the consensus
   was then rebuilt from those copies. The seed consensus is the authority:
   if the alignment matches its reverse complement better than itself, the
   whole alignment is turned back.

2. The element edge stops where flanks become unique, and does not take a
   polyA run with it. A run already sitting at the 5' end of the consensus
   is moved back into the flank. The 3' A-tail that is already inside the
   seed is not stripped.

Flanks are the ungapped sequence outside the current consensus, position 0
touching the edge. Run this on an alignment whose flanks are already packed
against the element (the published file), so a distance in bases is a column.
"""
import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import extend_unique_flank as U
from fix_alignments import justify

COMP = str.maketrans("ACGTNacgtn", "TGCANTGCAN")


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def ungapped(seq):
    return "".join(c for c in seq if c not in "-.").upper()


def best_identity(query, target):
    """Identity of the shorter sequence placed on the longer one."""
    if not query or not target:
        return 0.0
    q, t = query, target
    if len(q) > len(t):
        q, t = t, q
    if len(t) < len(q):
        return 0.0
    best = 0.0
    n = len(q)
    for i in range(0, len(t) - n + 1):
        same = 0
        for k in range(n):
            if q[k] == t[i + k]:
                same += 1
        ident = same / float(n)
        if ident > best:
            best = ident
            if best == 1.0:
                return best
    return best


def should_flip(cons, ref):
    body = ungapped(cons)
    fwd = best_identity(ref, body)
    rev = best_identity(revcomp(ref), body)
    return rev >= 0.65 and rev > fwd + 0.12, fwd, rev


def uppercase_bounds(seq):
    upper = [i for i, c in enumerate(seq) if c.isupper()]
    if upper:
        return upper[0], upper[-1]
    nz = [i for i, c in enumerate(seq) if c not in "-."]
    return (nz[0], nz[-1]) if nz else (0, max(0, len(seq) - 1))


def trim_leading_polya(cons, lo, hi):
    """Move lo past a 5' polyA. A base or two in front of that run goes too.

    The 3' A-tail is not touched. A real 5' end that merely contains an A
    is left alone: the polyA has to sit at the start, after at most a stub
    shorter than the run itself.
    """
    cols = [k for k in range(lo, hi + 1) if cons[k] not in "-."]
    i = 0
    stub = 0
    while i < len(cols) and cons[cols[i]].upper() != "A":
        stub += 1
        i += 1
        if stub >= U.POLYA:
            return lo
    run = 0
    while i < len(cols) and cons[cols[i]].upper() == "A":
        run += 1
        i += 1
    if run >= U.POLYA and i < len(cols):
        return cols[i]
    return lo


def fill_consensus(cons, lo, hi, lefts, rights, left_n, right_n):
    row = list(cons)
    maj_l = U.majority(lefts, left_n)
    for k in range(left_n):
        col = lo - 1 - k
        if 0 <= col < len(row) and maj_l[k] in "ACGT":
            row[col] = maj_l[k]
    maj_r = U.majority(rights, right_n)
    for k in range(right_n):
        col = hi + 1 + k
        if 0 <= col < len(row) and maj_r[k] in "ACGT":
            row[col] = maj_r[k]
    return "".join(row)


def correct(names, seqs, ref):
    ci = U.consensus_index(names)
    flipped = False
    fwd = rev = None
    if ref:
        flipped, fwd, rev = should_flip(seqs[ci], ref)
        if flipped:
            seqs = [revcomp(s) for s in seqs]
            names = [n[3:] if n.startswith("_R_") else n for n in names]
    cons = seqs[ci]
    lo, hi = uppercase_bounds(cons)
    lo2 = trim_leading_polya(cons, lo, hi)
    if lo2 != lo:
        row = list(cons)
        for j in range(lo, lo2):
            if row[j] not in "-.":
                row[j] = row[j].lower()
        cons = "".join(row)
        seqs[ci] = cons
        lo = lo2
    # flanks against the (possibly trimmed) edge
    seqs[ci] = cons
    lefts, rights, lo, hi = U.flanks_of(seqs, ci)
    # Both sides stop at a polyA run. The A-tail already inside the seed
    # stays; a run out in the flank is not added to the consensus.
    left_n = U.extension_edge(lefts, polya="stop")
    right_n = U.extension_edge(rights, polya="stop")
    left_n = min(left_n, lo)
    right_n = min(right_n, len(cons) - hi - 1)
    cons = fill_consensus(cons, lo, hi, lefts, rights, left_n, right_n)
    seqs[ci] = cons
    new_lo = lo - left_n
    new_hi = hi + right_n
    width = len(cons)
    out = []
    for s in seqs:
        s = s.ljust(width, "-")[:width]
        out.append(justify(s, new_lo, new_hi))
    info = {
        "flipped": flipped,
        "fwd": fwd,
        "rev": rev,
        "left": left_n,
        "right": right_n,
        "consensus": ungapped(out[ci]),
    }
    return names, out, info


def write_fa(path, names, seqs):
    ci = U.consensus_index(names)
    order = [ci] + [i for i in range(len(names)) if i != ci]
    with io.open(path, "w", encoding="utf-8") as fh:
        for i in order:
            fh.write(">%s\n" % names[i])
            s = seqs[i]
            for k in range(0, len(s), 80):
                fh.write(s[k:k + 80] + "\n")


def subfam_of(path):
    base = os.path.basename(path)
    # rsi_r1_9seqs_top100.aln.fa -> r1_9seqs
    parts = base.split("_")
    # rsi, r1, 9seqs, top100...
    if len(parts) >= 3 and parts[1].startswith("r") and parts[2].endswith("seqs"):
        return parts[1] + "_" + parts[2]
    return None


def load_bank(path):
    names, seqs = U.read_fa(path)
    return {n.split()[0]: ungapped(s) for n, s in zip(names, seqs)}


def ref_for(bank, path, names):
    """The forward seed for this alignment.

    The consensus header is the name SINEderella put in the file. The rsi
    filename pattern is only a fallback for a header that no longer matches.
    """
    ci = U.consensus_index(names)
    key = names[ci].split()[0]
    if key in bank:
        return bank[key]
    sf = subfam_of(path)
    if sf and sf in bank:
        return bank[sf]
    return None


def main(argv):
    write = "--write" in argv
    args = [a for a in argv[1:] if not a.startswith("--")]
    if len(args) < 2:
        raise SystemExit("usage: correct_published_aln.py BANK.fa ALIGNMENT.fa... [--write]")
    bank = load_bank(args[0])
    for path in args[1:]:
        names, seqs = U.read_fa(path)
        ref = ref_for(bank, path, names)
        names2, seqs2, info = correct(names, seqs, ref)
        # info['consensus'] includes lowercase flank; show the element only
        ci = U.consensus_index(names2)
        body = "".join(c for c in seqs2[ci] if c.isupper())
        print("%s  flip=%s fwd=%s rev=%s  L+%d R+%d" % (
            os.path.basename(path), info["flipped"],
            None if info["fwd"] is None else round(info["fwd"], 3),
            None if info["rev"] is None else round(info["rev"], 3),
            info["left"], info["right"]))
        print("  5': %s" % body[:50])
        print("  3': %s" % body[-40:])
        if write:
            write_fa(path, names2, seqs2)


if __name__ == "__main__":
    main(sys.argv)
