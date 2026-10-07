# CpG-corrected divergence: the problem, the proposed fixes, and what the rsi/bat data say (2026-10-07)

Prepared for decision D2 (should SINEderella's divergence profiles be CpG-adjusted?). Sources: PubMed records and abstracts (read, not full texts, except where
stated), the RepeatMasker source (`calcDivergenceFromAlign.pl`, `SearchResult.pm`, read in full), and measurements on the published plates (`cpg_nonCpG_ratio.tsv`,
script `cpg_ratio.py` in this folder). Every reference below was resolved in PubMed or the publisher record; where only an abstract was read, the text says "abstract".
This is a literature review plus one new measurement, not a proof that the correction is right or wrong for any one family.

## 1. The problem in one paragraph

Divergence of a repeat copy from its family consensus is the standard clock for TE age. That clock assumes substitutions accumulate at one rate. They do not: in
methylated genomes a cytosine in a CpG dinucleotide is often 5-methylated, and 5-methylcytosine deaminates to thymine (C>T on one strand, G>A on the other), a lesion
the repair machinery fixes less well than the uracil from ordinary cytosine deamination. CpG sites therefore decay roughly an order of magnitude faster than other sites,
and SINEs are CpG-rich, so a young CpG-rich SINE can look as old as a much older CpG-poor one. Two families of the same real age sit at different divergences, and the
divergence profile (the age axis of every SINE paper) is distorted in a family-specific way.

## 2. Mechanism and size of the effect

- **Mechanism.** Hotspots of base substitution sit at 5-methylcytosine and arise from its spontaneous deamination to thymine (Coulondre et al. 1978, E. coli, [10.1038/274775a0](https://doi.org/10.1038/274775a0));
  in animals 5-methylcytosine "tends to mutate abnormally frequently to T", the likely cause of CpG deficiency in heavily methylated genomes (Bird 1980, [10.1093/nar/8.7.1499](https://doi.org/10.1093/nar/8.7.1499)).
- **Size in humans.** Among collated disease mutations, 35 % lay in CpG dinucleotides and C>T / G>A there occurred 42-fold more often than random expectation (Cooper & Youssoufian 1988, [10.1007/BF00278187](https://doi.org/10.1007/BF00278187));
  from human-chimpanzee pseudogenes, transition and transversion rates at CpG are about one order of magnitude above those at other sites, mean rate about 2.5 x 10^-8 per site (Nachman & Crowell 2000, [10.1093/genetics/156.1.297](https://doi.org/10.1093/genetics/156.1.297)).
- **It is not one number.** Across 19 mammals, CpG transitions accumulated in a relatively clock-like way in contrast to other substitution types, and context-dependent replication errors, cytosine deamination and biased gene conversion contribute differently in different lineages (Hwang & Green 2004, [10.1073/pnas.0404142101](https://doi.org/10.1073/pnas.0404142101)).
  Rates also vary with the flanking bases and with methylation level, and differ between primate lineages (Chandra & Gao 2026, PLoS Genetics 22:e1011957, [10.1371/journal.pgen.1011957](https://doi.org/10.1371/journal.pgen.1011957)); general review: Hodgkinson & Eyre-Walker 2011, [10.1038/nrg3098](https://doi.org/10.1038/nrg3098).
  CpG transition mutability is the main driver of mutation-spectrum variation across 108 eukaryotes, yet genome-wide CpG methylation level does not predict the CpG transition rate across species (Ramos-Almodovar et al. 2025, bioRxiv preprint, not peer reviewed, title "Methylation-associated mutagenesis underlies variation in the mutation spectrum across eukaryotes"; DOI not resolved here).
  A mutational-signature model attributes the CpG signature mainly to DNA replication rounds in the germline (Spisak et al. 2024, PLoS Biology 22:e3002678, [10.1371/journal.pbio.3002678](https://doi.org/10.1371/journal.pbio.3002678)), which would tie it to cell divisions per generation.

## 3. How it hits TE age and subfamily work

- Ages of Alu subfamilies were calculated with Kimura's distance (Kapitonov & Jurka 1996, [10.1007/BF00163212](https://doi.org/10.1007/BF00163212); Kimura 1980, [10.1007/BF01731581](https://doi.org/10.1007/BF01731581)); the human Alu subfamily division itself rests on diagnostic positions (Jurka & Smith 1988, [10.1073/pnas.85.13.4775](https://doi.org/10.1073/pnas.85.13.4775)); overview in Batzer & Deininger 2002, [10.1038/nrg798](https://doi.org/10.1038/nrg798).
- Alu subfamily ages estimated from CpG or from non-CpG substitution density disagree (abstract of Xing et al. 2004, [10.1016/j.jmb.2004.09.058](https://doi.org/10.1016/j.jmb.2004.09.058)): 5,296 Alus in 20 subfamilies gave a roughly constant CpG : non-CpG substitution ratio of about 6 for young (AluY) and intermediate (AluS) subfamilies, but a non-linear relationship once old (AluJ) subfamilies were included. The authors suggest a slowdown of the neutral rate during primate evolution and/or a higher CpG rate after a burst of retrotransposition about 35 million years ago.
- Divergence from a derived consensus can be unreliable, especially for older elements (Giordano et al. 2007, [10.1371/journal.pcbi.0030137](https://doi.org/10.1371/journal.pcbi.0030137)); many "subfamilies" contain several source elements, so subfamily-consensus-based ages and mutation estimates need reassessment (Wacholder et al. 2014, PLoS Genetics 10:e1004482, [10.1371/journal.pgen.1004482](https://doi.org/10.1371/journal.pgen.1004482)).
- TE methylation is what lets TEs persist and is a cause of CpG loss in TEs; CpG O/E in TEs and host DNA falls as genome size and TE share rise across 53 organisms (Zhou et al. 2020, PNAS 117:19359, [10.1073/pnas.1921719117](https://doi.org/10.1073/pnas.1921719117)). In human spermatogenesis SINEs show differential methylation while LINEs appear protected (Siebert-Kuss et al. 2024, Am J Hum Genet 111:1125, [10.1016/j.ajhg.2024.04.017](https://doi.org/10.1016/j.ajhg.2024.04.017)): the methylation state of a SINE in the germline is therefore itself family- and stage-dependent.

## 4. The proposed solutions

| # | approach | what it does | cost / limit | source |
|---|---|---|---|---|
| A | none (raw p-distance or K2P) | leave CpG sites in | CpG-rich young families look old; family-specific distortion | default of many pipelines |
| B | drop CpG sites | count substitutions at non-CpG sites only | needs enough non-CpG sites; throws information away; right choice when CpG sites are saturated | used in Alu work (Xing et al. 2004 compare CpG and non-CpG density) |
| C | RepeatMasker CpG-adjusted Kimura | at consensus CpG sites two transitions count as one, one transition as 1/10 of a standard transition, transversions as usual; `-noCpGMod` switches it off | the 1/10 is a convention, one value for every family and species; "applicable to organisms for which CpG methylation is applicable" (RepeatMasker documentation, `util/calcDivergenceFromAlign.pl`) | RepeatMasker source |
| D | context-dependent substitution model | rates depend on the flanking nucleotides, branch and position; fitted by MCMC | needs a tree and a neutral alignment; not a per-copy distance | Hwang & Green 2004 |
| E | two clocks | CpG and non-CpG substitution density kept as two variables and their ratio studied | the ratio is not constant (old subfamilies) | Xing et al. 2004 |
| F | clocks that do not use divergence | relative age from TE-into-TE insertion order (defragmentation), or from explicit ancestry/source-element models | needs a whole genome / large copy numbers; no absolute scale | Giordano et al. 2007; Wacholder et al. 2014 |

## 5. When the correction should NOT be applied

The RepeatMasker rule presumes germline CpG methylation. Evidence that it does not carry over unchanged: in Drosophila and mosquito (very low methylation) CpG deficiency and TpG/CpA excess are similar to human (Jabbari & Bernardi 2004, Gene 333:143, [10.1016/j.gene.2004.02.043](https://doi.org/10.1016/j.gene.2004.02.043), abstract); CpG depletion in honeybee genes follows its methylation (Elango et al. 2009, PNAS 106:11206, [10.1073/pnas.0900301106](https://doi.org/10.1073/pnas.0900301106));
oyster CpG methylation targets young repeats (Wang et al. 2014, BMC Genomics 15:1119, [10.1186/1471-2164-15-1119](https://doi.org/10.1186/1471-2164-15-1119)). So "invertebrate" does not mean "no CpG effect", and a fixed 1/10 cannot be assumed: it has to be measured.

## 6. New measurement: what the CpG : non-CpG ratio is in these data

For every top-100 plate I counted substitutions per site at consensus CpG sites and at all other element sites (consensus defined as in RepeatMasker: C of a CG, and G of a CG), per copy, then averaged
(`cpg_nonCpG_ratio.tsv`, 74 plates; plates with fewer than 10 usable copies are left out). "Young" = non-CpG substitution density below 0.08.

| genome | plates | young plates | median CpG : non-CpG density, young plates | transitions among CpG-site changes |
|---|---|---|---|---|
| human (Alu and others) | 14 | 7 | **5.8** | 84-91 % |
| rsi (*Rhinolophus sinicus*) | 21 | 12 | **5.8** (r9 8.8, r7 8.2, r6 7.5, r8 5.8, r5 5.8, r10 5.8) | 78-92 % |
| rle (*Rousettus leschenaultii*) | 5 | 3 | **6.3** (MEG-TR 7.7, MEG-RL 6.3, MEG-RS 5.0) | 87-92 % |
| zebrafish | 13 | 7 | 2.7 | 72-91 % |
| hedgehog (eri) | 3 | 2 | 2.1 | n/a (few copies) |
| Timema (insect) | 14 | 3 | 1.5 | 66-98 % |

1. The human value (5.8) reproduces the published Alu figure of about 6 (Xing et al. 2004), which checks the measure.
2. Young bat SINE families behave like young Alus: CpG sites change 5-9 times faster than other sites, and almost all their changes are transitions. The RepeatMasker 1/10 is the right order, if slightly strong (1/6 would fit better).
3. Old families lose the signal. rsi MEG-T2, MEG-TR, MEG-RL, r2, r3, r4 and MEG-RS have ratios of 0.9-1.9 and CpG-site substitution density of 0.3-0.67: CpG sites are saturated (already decayed), so down-weighting them by 1/10 is the wrong operation; they carry no age information and should be dropped (approach B), not discounted.
4. Zebrafish and Timema show a much weaker CpG excess (2-3 and 1.5): a flat 1/10 over-corrects there.

Effect on the final rsi plates (SINEderella `tools/cpg_div`, same rule as RepeatMasker): median drop in divergence 27 %, 42 % of copies move to another of five divergence bins; r9 -45 %, r10 -41 %, r7 -38 %, r6 -35 %, r5 -34 %, against -7 % for r4, -8 % MEG-RL, -11 % MEG-T2 (`../rsi_decisions_2026-10-06/cpg_rsi_final.tsv`).
The 371-plate survey of 2026-10-01 is in SINEderella `docs/CPG_DIVERGENCE_TEST.md`.

## 7. Recommendation for D2 (not applied; the pipeline is unchanged)

1. Keep the raw divergence profile as the primary track.
2. Add a second track, and make it family-aware instead of using a universal 1/10: for each plate estimate the CpG : non-CpG density ratio r from the plate itself (as in section 6) and report it as a "CpG state" next to the divergence:
   r >= 4 young, the 1/10 rule (or 1/r) is appropriate; 2 <= r < 4 intermediate, use a milder weight 1/r; r < 2 saturated or not methylation-driven, report non-CpG-only divergence (drop CpG sites).
3. State in the manuscript that the correction presumes germline CpG methylation, and that the ratio is measured here (human 5.8, bats about 6, fish 2.7, insect 1.5), not assumed.
4. Never use consensus-defined CpG status for families whose consensus has itself lost its CpGs: the number of CpG sites in the consensus (14 of 135 positions in MEG-RS) must be reported with every adjusted value.
5. Treat the age ordering of families as the robust output: the order of the young rsi families stays almost unchanged (r7 0.0325 vs r9 0.0381 raw, 0.0201 vs 0.0210 adjusted, nearly tied), while absolute ages shift by up to 45 %.

## 8. Not verified / open

- Full texts were not read for the papers marked "abstract" or listed by link only; conclusions about them are limited to the abstract.
- I did not find, in this search, a paper that measures CpG hypermutability in bat germlines directly; the bat numbers above are from the plates only. (Searches of PubMed and a semantic index returned bat TE-landscape papers, e.g. Paulat et al. 2022 and Ricci et al. 2023, but nothing on bat CpG mutation rates.)
- The ratio uses the plate consensus as the ancestral state and ignores gaps, multiple hits and biased gene conversion; it is a descriptive measure, not a rate estimate. Plates carry the top 100 copies (the most similar to the consensus by construction), so young families are over-represented.
- Whether MEG-RS in rsi (array) and in rle (dispersed) differ in their CpG decay was not tested; they are close (1.6 vs 5.0), but the rsi array plate is dominated by old, saturated copies.
