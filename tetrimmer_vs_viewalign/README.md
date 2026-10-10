# TEtrimmer vs ViewAlign: clustering and cleaning, and what is worth borrowing (2026-10-11)

Sources. (1) TEtrimmer paper, Qian et al. 2025 Nat Commun, [10.1038/s41467-025-63889-y](https://doi.org/10.1038/s41467-025-63889-y), full text read via Europe PMC (PMC12462492), Methods v1.5.1.
(2) The source-level profile of v1.7.4 in the SubFam repo, kept here as `TETRIMMER_SubFam_profile.md` (read, not re-verified against the TEtrimmer code).
(3) ViewAlign at `C:\work\MSA-viewer-new` (HEAD 5d637fb): `cluster.js`, `peel.js`, `kmer-tree.js`, `block-mask.js`, `block-bicluster.js`, `script.js`, `index.html`.
(4) One small test, `crop_compare.js` (this folder), on four published plates. It is my re-implementation from the Methods text, not TEtrimmer's code.

## 1. What each tool is for

TEtrimmer curates a library: one seed consensus in, genome BLAST hits, extension, cleaning, boundaries, one graded consensus out. It needs a genome (`--genome_file` is required, profile section 1), so it cannot be pointed at an existing alignment of copies.
ViewAlign is an alignment workbench: the alignment is the input, and the user decides. They overlap in two places: splitting an alignment into variants, and cleaning an alignment before a consensus is made.

## 2. Version drift inside TEtrimmer (do not quote one number as "the" default)

| item | paper (v1.5.1) | SubFam profile (v1.7.4) |
|---|---|---|
| copies per MSA | 70 longest + 30 random | 100 longest, then random fill |
| DBSCAN min samples | 2 (text and Fig. 2) | 3 (`MSAcluster.py:366`) |
| Noise kept as a cluster | more than 15 sequences and over 60 % | at least 12 sequences and at least 60 % |
| clustering gate | not stated | more than max(50, 5 % of length) divergent columns |

I did not check which is right for the released code; the point is that the thresholds below move between versions.

## 3. Clustering, side by side

| step | TEtrimmer | ViewAlign |
|---|---|---|
| columns used | divergent columns (major allele below 0.8) plus "gap blocks" (indel regions flanked by over 80 % conserved bases, adjacent gap fractions within 10 %) | Peel: diagnostic columns by your definition (one character, gap included, fills at least 90 % of members and at most 2 % outside; N/IUPAC missing). SNP grouping: exclusive columns with 5 -> 1 relaxation. Indels get weight x2 in the refine pass |
| distance | ML tree on those columns (IQ-TREE K2P+I), patristic distance | aligned-column p-distance, or k-mer Jaccard (canonical option); UPGMA, no model |
| number of groups | DBSCAN, eps 0.1, so "how far apart is far" is a fixed relative branch length; at most 5 clusters; floors 10 / 18 / 20 sequences | tree cut by most persistent count, or Peel loop with no cap; minimum group size 3 |
| leftovers | noise dropped (or used if it is the majority) | left unassigned, or returned to the pool (misfit return), shown, not deleted |
| separate outliers | part of DBSCAN noise | explicit outlier step in Peel (no neighbour as close as a typical nearest neighbour), set aside first and after each peel |
| scale | needs a MSA and IQ-TREE; authors say it cannot handle thousands of sequences | thousands of rows in the browser; deterministic |
| validation | PCA of the tree distance (Fig. 2C) | scored against your calls: ccr ARI 0.55, oma 0.69-0.84; 13 calibration cases |

Reading. For the same job, ViewAlign is already further along: no fixed epsilon, no cluster cap, explicit outliers, deterministic, and your own diagnostic-column rule is a better feature definition than "major allele below 0.8" (it counts a group-specific character, not variability). TEtrimmer's clustering is built to protect a consensus from being a chimera, not to find subfamilies.
Two ideas are worth a test, neither a replacement:
- **Gap blocks as features.** TEtrimmer treats a contiguous run of gap columns as one indel feature, so a 12-bp indel counts once, not 12 times. Your rule already says "indels have higher diagnostic value"; the refine pass doubles their weight, but each gap column of a long indel still votes separately. Collapsing an indel block into one feature, and ranking it above a SNP, would match what you said. Cheap to try on the 13 calibration cases.
- **A noise-aware fallback.** When most sequences are old and fall in no group (their Noise cluster above 60 %), TEtrimmer analyses the noise as one group. Peel's "remainder as a group only if as tight as the peeled groups" is the same idea with a test, so no change is needed; it is a confirmation.

## 4. Cleaning, side by side (this is where TEtrimmer has something ViewAlign lacks)

| operation | TEtrimmer | ViewAlign |
|---|---|---|
| gappy columns | removed when over 80 % gap or fewer than 5 nucleotides; then 40-80 % gap with major allele below 70 % | "Empty cols" removes only columns that are a gap in every sequence (`script.js:24130`) |
| ragged ends | per row: a 40-column window slides in from each end until the mean column proportion of that row's residues reaches 0.7; everything before it is deleted; second pass window 4 at 1.0 | one cut for all rows: a sliding window over gap percentage (Left 50 %, Right 80 %, window 15), `cluster.js:975`; Hard deletes, Soft hides from grouping |
| gap-based end crop | `crop_end_by_gap`: per row, window 250, stops when gaps below 10 %; used for truncated LINE 5' ends | same as the global cut above, not per row |
| boundaries | terminal repeats by self-BLAST, then poly(A)/microsatellite tail for LINE/SINE, then ambiguity window | TSD finder with border slack and column-conservation scoring (`_findSineBoundaryColumns`), reference-row or manual |
| extension into flanks | iterative, genome needed | none (alignment only; the flank is whatever the alignment contains) |
| near-duplicate input | CD-HIT-EST at 95 % identity and coverage (`--dedup`) | I could not find a redundancy filter in `script.js` or `index.html` (only an exact same-name-and-source check at `script.js:1343`) |
| consensus | per column, A/C/G/T only, winner at 0.7 (0.65 divergent), else N; gap-free | threshold plus a separate coverage minimum, plurality or IUPAC, gap or keep-best fallback; group consensus; replace with consensus |
| grading, reports | Perfect / Good / Reco_check / Need_check, PDF plots, GUI | per-plate statistics and verdict live in SINEderella, not in ViewAlign |

The gap that matters: **ViewAlign cannot clean one row's ragged end without touching every other row.** TEtrimmer's per-row crop is what a curator does by eye when one copy has a non-homologous flank and the others do not.

## 5. Test: does TEtrimmer's per-row crop work on SINE plates? (`crop_compare.js`)

Four published top-100/rand100 plates, 102 rows each. "Agreement" is the mean, over residues kept, of the share of the column that carries the same base.

| plate | cols | ViewAlign global trim (kept residues, agreement) | TEtrimmer gap columns + pass 1 only (kept residues) |
|---|---|---|---|
| rsi MEG-RS rand100 | 3104 | 97.0 %, 0.784 | 79.3 % |
| rle MEG-RS rand100 | 438 | 97.4 %, 0.477 | 37.2 % |
| rsi r6 top100 | 516 | 97.5 %, 0.620 | 56.5 % |
| rsi P18 rand100 | 1367 | 98.6 %, 0.605 | 57.4 % |

- The ViewAlign trim removes 1-3 % of residues: it takes only empty ends, which is what it is meant to do.
- With its default thresholds (window 40, mean proportion 0.7), TEtrimmer's pass 1 deletes 21-63 % of the residues on these plates. The SINE plates have mean column agreement 0.5-0.8, so a 40-column window rarely reaches a mean of 0.7: the thresholds were set for young, well-conserved elements and delete real SINE sequence on old ones. My second pass (window 4, threshold 1.0, with "reached or exceeded" read as at least) deleted everything on three of four plates; I do not know whether the real code does the same, so I make no claim about pass 2.
- So the idea (per-row crop by column conservation) is sound, but the threshold has to be relative to the plate (for example a quantile of the plate's own agreement), not a constant 0.7.

Limits. One implementation written from the Methods text, four plates, no manual reference for "correct cut". Not a benchmark of TEtrimmer, and the numbers say nothing about it on the genomes it was designed for.

## 6. What I would borrow, ranked

1. **Per-row end crop by column conservation** (TEtrimmer `crop_end_by_divergence`), as an option next to Trim, with the threshold set relative to the plate and a preview in the same style (soft by default: dim the cropped residues, do not delete). Use: remove a non-homologous flank from single rows before a consensus or a peel. New code; about the size of `getTrimBoundaries`.
2. **Threshold-based gappy-column removal** (over 80 % gap or under 5 residues, optional second rule), next to "Empty cols". Trivial. Needed so that group consensuses and column-based distances are not driven by columns that only two rows fill.
3. **Indel blocks as single features** in Peel's diagnostic-column count (see section 3). Test on the 13 calibration cases first; keep only if ARI does not fall.
4. **Redundancy filter** (collapse sequences identical, or over a chosen % identity, in the aligned columns). Your earlier idea "double all sequences" for group-size tolerance is the opposite direction, so this should be an option, off by default. Not in TEtrimmer's clustering, only for its input.
5. **Grade-style summary per group** (size, number of diagnostic columns, fraction of members with an exception, mean agreement): ViewAlign already has most of these numbers in the peel explainer; a one-line grade would make them comparable between groups. Low priority.

Not worth borrowing: IQ-TREE/DBSCAN (needs a model, a tree and a fixed epsilon; your tree is already deterministic and fast), genome-dependent extension, Pfam/ORF classification (not relevant to SINEs).

## 7. Open points for you

- Do you want items 1-2 built (they touch the Trimming panel), or only noted?
- For item 1, the right test is your manual cleaning of a plate: if you can name two or three alignments where you removed a flank from single rows by hand, I can score the crop against those, as with the peel calls.
