# rsi decisions, 2026-10-06 - direct links to the data

Alignment links open in ViewAlign; table/file links open on GitHub. Whole run: [rsi_final README](https://github.com/Toki-bio/Tal/blob/main/rsi_final/README.md), [report.html (download, open locally)](https://github.com/Toki-bio/Tal/blob/main/rsi_final/report.html).

## 1. MEG-RS: tandem-array family, or dispersed elsewhere? (TESTED)

Share of MEG-RS copies lying in tandem arrays (`tools/array_flag.py`, same code on each existing run):

| genome | MEG-RS copies | in arrays | arrays | largest | median spacing |
|---|---|---|---|---|---|
| rsi (*Rhinolophus sinicus*) | 1,715 | **92.8 %** | 21 | 307 | 2,160 bp |
| rle (*Rousettus leschenaultii*, megabat) | 6,718 | **17.9 %** | 28 | 141 | 899 bp |
| rda (*Rhinolophus*) | 10 | 0 % | 0 | - | - |
| rre (*Rhinolophus*) | 21 | 0 % | 0 | - | - |

So in rle MEG-RS is a dispersed family (82 % of 6,718 copies are not in arrays); the tandem-array state is specific to rsi. The other two Rhinolophus genomes have almost no MEG-RS (earlier runs, 2026-09-27; they are being re-run with the 21-consensus bank, see 7).

- rsi top 100: [rsi MEG-RS top100](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rsi_final/alignments/rsi_MEG-RS_top100.aln.fa&title=rsi%20MEG-RS%20top100) (32 of 100 rows are marked `[array]`); random 100: [rsi MEG-RS rand100](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rsi_final/alignments/rsi_MEG-RS_rand100.aln.fa&title=rsi%20MEG-RS%20rand100)
- rle (dispersed): [rle MEG-RS top100](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rle/alignments/rle_MEG-RS_top100.aln.fa&title=rle%20MEG-RS%20top100), [rle MEG-RS rand100](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rle/alignments/rle_MEG-RS_rand100.aln.fa&title=rle%20MEG-RS%20rand100)
- array table of rsi: [array_flag.tsv](https://github.com/Toki-bio/Tal/blob/main/rsi_final/array_flag.tsv); MEG-RS vs MEG-RL consensus pair: [MEG-RS vs MEG-RL](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rsi_final/consensus_pairs/MEG-RS__MEG-RL.aln.fa&title=MEG-RS%20vs%20MEG-RL)
- length-version control in rle (Gogolevsky 2009): [report](https://github.com/Toki-bio/Tal/blob/main/chiroptera/length_variants/rle_MEG-RS_MEG-RL.report.txt)
- The README of rsi_final lists `MEG-RS_array_unit_2155bp.fa`; that file is NOT in the repository (commit 640af55 added only the README line). I could not find it on therioserver either.

## 2. r8 + r8 dimer (P42 + P43)

- Combined alignment: [r8 dimer P42+P43](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rsi_v7/r8_dimer_P42_P43.aln.fa&title=r8%20dimer%20P42%2BP43) (same file from the earlier run: [earlier run](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rsi_final/r8_dimer_P42_P43_earlier_run.aln.fa&title=earlier%20run))

## 3. The four new unit pairs that pass the 70 % rule

- r10 + P48: [r10 + P48 (P10)](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rsi_final/flankscan/cand/r10_19seqs__r5_r5_P48_P10.aln.fa&title=r10%20%2B%20P48%20%28P10%29)
- P48 + r8: [P48 + r8 (P15)](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rsi_final/flankscan/cand/r5_r5_P48__r8_83seqs_P15.aln.fa&title=P48%20%2B%20r8%20%28P15%29)
- C11 + P1: [C11 + P1 (P3)](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rsi_final/flankscan/cand/r1_r2_r3_r3_C11__r1_r3_P1_P3.aln.fa&title=C11%20%2B%20P1%20%28P3%29), [C11 + P1 (P4)](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rsi_final/flankscan/cand/r1_r2_r3_r3_C11__r1_r3_P1_P4.aln.fa&title=C11%20%2B%20P1%20%28P4%29)
- P26 + P34: [P26 + P34 (P8)](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rsi_final/flankscan/cand/r5_r6_P26__r5_r3_P34_P8.aln.fa&title=P26%20%2B%20P34%20%28P8%29)

## 4. Peel the subfamilies, then rerun with `--add`

- The 600 rsi SubFam chunk consensuses you peeled (10 calls): [rsi_subfam_input_30k](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rhin/alignments/rsi_subfam_input_30k.aln.fa&title=rsi_subfam_input_30k) - calls: https://github.com/Toki-bio/SINE-discriminator (`RSI_PEEL_LOG.md`)
- Subfamily plates of the current run, per family: [r6 subfam](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rsi_final/alignments/rsi_r6_210seqs_subfam.aln.fa&title=r6%20subfam), [r7 subfam](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rsi_final/alignments/rsi_r7_133seqs_subfam.aln.fa&title=r7%20subfam), [r5 subfam](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rsi_final/alignments/rsi_r5_27seqs_subfam.aln.fa&title=r5%20subfam)

## 5. G-rich 3' end of P18

- [P18 top100](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rsi_final/alignments/rsi_r10_r8_P18_top100.aln.fa&title=P18%20top100), [P18 subfam](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rsi_final/alignments/rsi_r10_r8_P18_subfam.aln.fa&title=P18%20subfam), [P18 rand100](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rsi_final/alignments/rsi_r10_r8_P18_rand100.aln.fa&title=P18%20rand100)

## 6. CpG-corrected divergence (TESTED on the final plates)

Same tool as the 2026-10-01 survey of 414 plates (`SINEderella/tools/cpg_div`, RepeatMasker rule: a CpG transition counts 1/10), run on the 21 top-100 plates of rsi_final. Table with every family: `cpg_rsi_final.tsv` (this folder).

| family | raw median divergence | CpG-adjusted | drop | copies changing divergence bin (of 5) |
|---|---|---|---|---|
| r9 | 0.0381 | 0.0210 | 45 % | 55 % |
| r10 | 0.0753 | 0.0441 | 41 % | 49 % |
| r7 | 0.0325 | 0.0201 | 38 % | 56 % |
| r6 | 0.0444 | 0.0289 | 35 % | 52 % |
| r5 | 0.0565 | 0.0373 | 34 % | 45 % |
| r8 | 0.0556 | 0.0403 | 27 % | 47 % |
| r1 | 0.0850 | 0.0667 | 22 % | 12 % |
| r3 | 0.1676 | 0.1347 | 20 % | 14 % |
| MEG-RS | 0.0435 | 0.0357 | 18 % | 2 % |
| MEG-T2 | 0.4000 | 0.3575 | 11 % | 40 % |
| MEG-RL | 0.2875 | 0.2650 | 8 % | 0 % |
| r4 | 0.3892 | 0.3623 | 7 % | 0 % |

Median over the 21 plates: divergence drops 27 %, and 42 % of copies move to a different divergence bin. The young, CpG-rich families (r9, r10, r7, r6, r5) are affected most (35-45 %), the old ones least (7-12 %); the order of the young families is mostly kept (r7 stays slightly younger than r9: raw 0.0325 vs 0.0381, adjusted 0.0201 vs 0.0210, nearly tied). Full survey and caveats: [docs/CPG_DIVERGENCE_TEST.md](https://github.com/Toki-bio/SINEderella/blob/main/docs/CPG_DIVERGENCE_TEST.md) (371 plates: median drop 10 %, mammals 27-35 %).

## 7. r4 / r2 on another species (RUNNING)

The 21-consensus bank of rsi_final is being run on the other two Rhinolophus genomes (therioserver `~/rhin/rda_cross`, `~/rhin/rre_cross`, started 2026-10-06 12:48, about 70 min each). Until they finish, the data for the question: the r4 / r2 pair alignment [r4 vs r2](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rsi_final/consensus_pairs/r4_32seqs__r2_3seqs.aln.fa&title=r4%20vs%20r2) and their rsi plates [r4 top100](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rsi_final/alignments/rsi_r4_32seqs_top100.aln.fa&title=r4%20top100), [r2 top100](https://toki-bio.github.io/MSA-viewer/?url=https://raw.githubusercontent.com/Toki-bio/Tal/main/rsi_final/alignments/rsi_r2_3seqs_top100.aln.fa&title=r2%20top100).

## 8. MEG-RS copies that share flanks

After the arrays are removed (current code, run of 2026-10-06): 96 MEG-RS copies assigned, 77 flank twins, 3 masked, 16 unique (candidates for independent insertions). Details: [OLD_VS_NEW_RUNS section 5c](https://github.com/Toki-bio/Tal/blob/main/OLD_VS_NEW_RUNS_2026-10-05.md).
