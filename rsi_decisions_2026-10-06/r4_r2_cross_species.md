# r4 and r2 on two other Rhinolophus genomes (run 2026-10-06/07)

The 21-consensus bank of `rsi_final` was run through SINEderella (same code and settings as the rsi run) on the other two Rhinolophus assemblies on therioserver:
`~/rhin/rda_cross/run_20261006_124835` (96 min) and `~/rhin/rre_cross/run_20261006_124835` (142 min). Tables: `rda_rre_assignment_and_audit.txt` (assigned copies per family, length-version table, consensus audit).

| family | rsi (rsi_final) | rda | rre |
|---|---|---|---|
| r4_32seqs, assigned copies | 12 | 13 | 14 |
| r2_3seqs, assigned copies | 16 | 14 | 16 |
| r3_58seqs, assigned copies | 49 | 46 | 52 |
| MEG-RS | 1,715 | 34 | 26 |
| MEG-RL | 12 | 12 | 12 |
| r9_15seqs | 6,844 | 7,129 | 6,810 |

1. **r4 and r2 are equally rare in all three genomes** (12-14 and 14-16 assigned copies). The length-version test therefore stays NOT_TESTED for r4/r3 and MEG-RS/MEG-RL in rda and rre as well (it needs 100 copies per family): another species with the same Rhinolophus history does not help; a species with more copies of these families is needed.
2. **Shared flanks (crude, `cross_flank.py`):** share of one species' plate copies that have a 16-mer in common with a flank of a copy in the other species (either strand), top-100 plates (which contain 25-35 copies for r4 and 56-73 for r2, more than the assigned counts):
   r4: rsi-rda 12 %, rsi-rre 32 %, rda-rre 11 %; r2: 60 %, 68 %, 71 %; r3 (100 copies each): 62 %, 57 %, 59 %.
   So r2 behaves like r3 (about 60 % of copies sit in flanks that are also found in the other species), r4 does not (11-32 %).
   Caveat: this is not an orthology test. Flanks that are themselves repeats would be shared as well, and no control family was run; treat the numbers as a pointer for the r4 and r2 question, not an answer.
3. **Consensus audit** in rda and rre: r4 and r2 SKIPPED (too few copies); MEG-RS SHORTER (rda rebuilt 89 bp, rre 120 bp against 135 bp), r7 LONGER (as in rsi), r5_r5_P48 DIVERGED in rda, MATCH in rre.
