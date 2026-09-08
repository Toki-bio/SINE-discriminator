#!/usr/bin/env python3
"""Spot-check row-0 element on DRAGEN after rebuild."""
import sys
from pathlib import Path

sys.path.insert(0, "/staging/tmp/sinedisc")
from fix_alignments import read_fa

ALN = Path("/staging/tmp/scorpions/oma/run_oma/alignments")
for name in sys.argv[1:] or ["oma_SINE10_top100", "oma_SINE18_top100"]:
    p = ALN / (name + ".aln.fa")
    _, seqs = read_fa(str(p))
    e = "".join(c for c in seqs[0] if c.isupper())
    print(name, len(e), e[:12], "|", e[-12:])
