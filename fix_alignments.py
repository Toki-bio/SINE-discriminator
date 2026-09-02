#!/usr/bin/env python3
"""Apply his display convention to every alignment: flanks degapped, not aligned.

Two defects he pointed at on NEGTRUNC5__saq__s5_5seqs:

  1. the consensus is over-extended on the left and not supported by the copies
  2. the left flank is aligned, with gaps through it, instead of degapped

They are one defect. His postprocess_flanks() in extract_alignments.sh defines
the element as the consensus row's first and last non-gap column, and degaps
everything outside it. That is correct only when the consensus is right. On
s5_5seqs the consensus's leftmost 12 bases sit in columns where only 34 of 100
copies have any base at all, so the body boundary lands 40+ columns too far
left, and the true left flank ends up INSIDE the body - where his degapping
never reaches. Measured across the corpus: the consensus spans 1.75x its own
length even in clean POS sets (43 % of the body is interior gap columns), so
this is systematic, not one bad file.

So: find the element edges by copy SUPPORT, then apply his rule.

  support(col) = fraction of copies with a non-gap base in that column

Walk in from each end over the consensus's own non-gap columns and stop at the
first one reaching MIN_SUPPORT. Everything outside becomes flank and is
degapped, lowercased and butted against the element exactly as his awk does:
internal gaps removed, bases pushed to the element edge, outer side padded so
the column count never changes.

The consensus bases that get trimmed are NOT deleted - they move into the flank
and are lowercased with everything else, so an over-extended consensus stays
visible as lowercase leader rather than being silently hidden.
"""
import glob
import io
import os
import sys

MIN_SUPPORT = 0.50      # a consensus column is real if half the copies reach it
GAPS = "-."


def read_fa(p):
    names, seqs, cur, buf = [], [], None, []
    for line in io.open(p, encoding="utf-8", errors="replace"):
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


def consensus_index(names):
    for i, h in enumerate(names):
        if "CONSENSUS" in h.upper():
            return i
    return 0


def element_bounds(cons, others):
    """His rule (first/last non-gap consensus base), corrected by copy support."""
    nz = [i for i, c in enumerate(cons) if c not in GAPS]
    if not nz or not others:
        return (nz[0], nz[-1]) if nz else (0, len(cons) - 1)
    n = float(len(others))
    sup = {}
    for j in nz:
        sup[j] = sum(1 for s in others if j < len(s) and s[j] not in GAPS) / n

    lo = next((j for j in nz if sup[j] >= MIN_SUPPORT), nz[0])
    hi = next((j for j in reversed(nz) if sup[j] >= MIN_SUPPORT), nz[-1])
    if hi <= lo:                       # nothing supported - keep his plain rule
        return nz[0], nz[-1]
    return lo, hi


def justify(seq, lo, hi):
    """His postprocess_flanks, verbatim in behaviour.

    left flank  : gaps dropped, lowercased, right-justified against the element
    right flank : gaps dropped, lowercased, left-justified against the element
    body        : forced UPPERCASE - see below
    width       : unchanged

    His awk lowercases the flanks and leaves the body alone, which is enough in
    his pipeline because his extracted copies arrive uppercase. This corpus does
    not: NEGLINEORF__teu__r00 is 98788 lowercase against 19674 uppercase because
    the genome's soft-masking was carried straight through, while
    NEGTRUNC5__saq__s5_5seqs is uppercase throughout. If the body is left as-is,
    lowercase stops meaning "flank" and the viewer's case cue is noise. So the
    body is uppercased explicitly and only the flank is lowered.
    """
    left_len = lo
    lf = "".join(c.lower() for c in seq[:lo] if c not in GAPS)
    lf = "-" * (left_len - len(lf)) + lf

    body = seq[lo:hi + 1].upper()

    right_len = len(seq) - hi - 1
    rf = "".join(c.lower() for c in seq[hi + 1:] if c not in GAPS)
    rf = rf + "-" * (right_len - len(rf))

    return lf + body + rf


def fix(path, out_path):
    names, seqs = read_fa(path)
    if len(seqs) < 2:
        return None
    ci = consensus_index(names)
    cons = seqs[ci]
    others = [s for i, s in enumerate(seqs) if i != ci]

    nz = [i for i, c in enumerate(cons) if c not in GAPS]
    if not nz:
        return None
    lo, hi = element_bounds(cons, others)
    trimmed_l = sum(1 for j in nz if j < lo)
    trimmed_r = sum(1 for j in nz if j > hi)

    out = []
    for i, s in enumerate(seqs):
        s = s.ljust(len(cons), "-")
        # the consensus is justified too: its trimmed ends become lowercase
        # flank, so an over-extended consensus stays visible instead of vanishing
        out.append(justify(s, lo, hi))

    with io.open(out_path, "w", encoding="utf-8") as fh:
        for h, s in zip(names, out):
            fh.write(u">%s\n" % h)
            for k in range(0, len(s), 80):
                fh.write(s[k:k + 80] + u"\n")
    return trimmed_l, trimmed_r, len(nz)


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "alignments"
    files = sorted(glob.glob(os.path.join(d, "*.aln.fa")))
    n = 0
    tl = tr = 0
    worst = []
    for p in files:
        r = fix(p, p)
        if not r:
            continue
        n += 1
        a, b, ln = r
        tl += a
        tr += b
        if a + b:
            worst.append((a + b, a, b, ln, os.path.basename(p)))
    worst.sort(reverse=True)
    print("rewrote %d alignments" % n)
    print("consensus bases moved into the flank: %d left, %d right" % (tl, tr))
    print("alignments with an over-extended consensus: %d" % len(worst))
    print("\nworst 15:")
    print("  %-46s %5s %5s %6s" % ("set", "left", "right", "conslen"))
    for t, a, b, ln, nm in worst[:15]:
        print("  %-46s %5d %5d %6d" % (nm[:46], a, b, ln))


if __name__ == "__main__":
    main()
