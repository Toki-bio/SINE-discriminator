# SINE discriminator — project state, for whoever picks this up next

Written 2026-08-31, end of build 2. This is the document to read first. It is
written for a session that has none of the earlier context.

- `SINE_discriminator_spec.md` — Sergei's original design. Still the reference,
  with the amendments in §6 below.
- `FINDINGS.md` — the chronological record, including two retractions. Read §11
  before believing anything in §3–§5.
- `PIPELINE.md` — how the corpus is built (build 1 geometry; superseded flank size).
- `HANDOFF.md` — this file.

---

## 1. What the project is

Given a set of genomic loci proposed as one SINE family — the output of a
`sear`/`sear2k` search from a consensus, or of a tool like AnnoSINE — decide
whether it is a real family, refine its boundaries, and say which copies are
contaminants. The manual criterion being formalised is the three-tier eyeball
test in spec §1.1: alignability, a sharp three-part boundary structure, and
internal homogeneity that is not mosaic.

## 2. The two corrections that define the current design

Build 1 got both of these wrong and produced a spurious "blocker". They are the
most important things in this document.

**The element is delimited by the known subfamily consensus, added to the
alignment as a row.** Not a consensus recomputed from the alignment, and not a
threshold on any profile. Its first and last nucleotide are the boundary. No
free parameters, and nothing that moves when the copies are resampled. This is
already the Tal site's own convention — its published `rand100` alignments carry
the consensus as row 1.

Do not read the boundary as the consensus row's *span in alignment columns* —
copies' insertions stretch that from 1.2x its true length at n=25 to 3.9x at
n=200. The readout is each copy's ungapped base count between the consensus's
first and last nucleotide.

**Flanks must be short — 70 bp, not the spec's 300.** 300 bp forces MAFFT to
align 600 bp of unrelated DNA per copy and invent thousands of junk columns: a
250 bp element gets smeared over 1200–1475 alignment columns, and the smearing
grows with copy number. That was the sole cause of build 1's boundary
instability. At 70 bp the median alignment is 942 columns instead of 3025.

The background identity never needed to be in the alignment. It is measured
separately, on ungapped sequence, and equals ~0.25 for unrelated DNA.

Effect on reproducibility, same side quest both times:

| boundary rule | SD of element length over K=20 subsets |
|---|---|
| build 1, profile crossing | 74.2 |
| count of columns with occupancy > 0.9 | 13.3 |
| **consensus row, 70 bp flanks** | **10.7 bp** |

## 3. What works, measured on the build-2 corpus

AUC against real families (POS + side quest, 428 sets). 1.000 or 0.000 is
perfect separation; the direction only says which way the feature runs.

| failure mode | feature | AUC |
|---|---|---|
| random genomic loci | `cliff`, `cons_identity_med` | 1.000 |
| 10 % contamination | `frac_supported` | 0.996 |
| 30 % contamination | `frac_supported` | 0.998 |
| 5′ truncation | `elem_len_cv` | 0.000 |
| 5′ truncation (independent) | `res_asymmetry` | 0.001 |
| chimera | `tsd_len_med` | 0.956 |

Two of the spec's own ideas were dismissed in build 1 and work now, both because
they need exact edges, which only the consensus anchor provides:

- **`res_asymmetry`** — the ragged-5′-edge LINE signature. Separated nothing in
  build 1; AUC 0.001 against truncation now.
- **`tsd_frac`** — 0.74 of copies in real families carry a detectable TSD,
  median 9 bp, against 0.04 in truncated sets and 0.07 in random loci. Build 1
  reported this as noise.

`flank_id` reads 0.30 in every class, i.e. the true unrelated-DNA background.
Build 1's 0.50 was an artifact of comparing sequences after aligning them.

## 4. What is not to be trusted

- **`rank1_excess` is confounded.** It only computes after imputing missing
  windows, and the imputed fraction tracks class directly (0.26 real, 0.55
  random). Part of the separation is the imputation. The spec's central Tier-3
  statistic is still not honestly tested.
- **`NEGSPLICE` is mis-designed.** Every chimera splices at the same relative
  midpoint, producing a coherent new family rather than a mosaic. It is caught
  by TSD destruction, not by mosaicism detection. Rebuild with per-copy-varying
  breakpoints before drawing any Tier-3 conclusion.
- **`NEGJITTER` is not a negative.** ±20–100 bp edge displacement changes
  nothing biological; it is caught only by `flank_bp_med`, which is mechanical.
  Treat it as a robustness control. There is currently no real boundary-shift
  negative.
- **All negatives except `NEGRAND` are synthetic.** No LINEs, satellites or
  processed pseudogenes yet. Do not fit a model until they exist.

## 5. Prior art in the Tal repo that must not be re-derived

`eri/LOG.md` (34 KB) is the important one and was missed on the first pass.

- Flank-uniqueness background established empirically: 2 000 pairs from 300
  random unrelated 70 bp genomic windows → **25.0 % mean pairwise identity**.
  Uniqueness threshold ≤ 35 % on both flanks.
- A recorded failure: testing only for *exact-duplicate* flanks passed a cluster
  that actually had 66.5 % mean pairwise right-flank identity with a shared
  motif. Do not use exact-match flank tests.
- Contaminated clusters read 0–1.8 % singleton-flank; clean ones 60–82 %.
- Periodicity in SINE-hit spacing — essentially the spec's `coord_cluster` — was
  tested against direct TRF tandem annotation and found **unreliable**, wrong on
  four of ten tribes in both directions. Use TRF intersection instead.

`*/tribes/tribes_metadata.txt` documents the site's alignment convention:
50 bp left / 70 bp right flanks, attached as unaligned blocks to a body-only
MAFFT alignment. That architecture is the right one and predates this work.

## 6. Amendments to propose to the spec

1. §2 "widen flanks to 300 bp; costs nothing" — false. It degrades the
   alignment measurably. Use 70 bp, measure background separately.
2. §2 per-family permutation null by flank reassignment — degenerate as stated.
   Rows in the flank block are exchangeable, so any row-symmetric statistic is
   invariant under it.
3. §2/§4 `rank1_frac` — saturated (0.99 for everything). Needs the null
   subtracted, and the null needs to be unconfounded by imputation.
4. §4 boundary from consensus IC crossing — replace with the known-consensus
   row, which has no free parameter.
5. §4 `cliff_z` — sample-size dependent, `corr(n_copies, cliff_z) = −0.34`.
   Prefer the raw `cliff`.
6. §F `coord_cluster` for satellites — the hedgehog log shows it is unreliable;
   use TRF overlap.
7. §7 open item "should MIXED be a third verdict class" — yes. Contamination
   destroys the analysis before any model runs: at 10 % contamination build 1's
   pipeline failed outright in 57 % of sets.

## 7. Where everything is

**Local — `C:\work\SINE_discriminator\`**

| path | what |
|---|---|
| `measure_c.py` | MEASURE(), consensus-anchored. The current implementation. |
| `profiles.py` | positional profiles (per-nucleotide tracks) |
| `analyze.py` | feature separation tables; reads a features jsonl |
| `build_viewer_data.py` | picks the sets worth eyeballing, packs them for the page |
| `features_c.jsonl` | feature vectors, 580 sets, build 2 |
| `profiles.json`, `summary_c.json` | data behind the report page |
| `aln_c/` | 580 consensus-anchored alignments (current) |
| `aln/` | 580 build-1 alignments, 300 bp flanks — superseded, kept for comparison |
| `measure.py`, `features_v2.jsonl` | build 1 — superseded, do not use |
| `site/` | the GitHub Pages tree, ready to push |
| `report2.html` | the build-2 report; `site/index.html` is a copy |
| `data/Tal/` | clone of github.com/Toki-bio/Tal |

**KIT — `/data/W/toki/SINE_disc/`**

| path | what |
|---|---|
| `sets/` | unaligned sets, 300 bp flanks (build 1) |
| `sets_c/` | unaligned sets, consensus row + 70 bp flanks (build 2) |
| `aln_c/` | build-2 alignments |
| `manifest.json`, `prep.log`, `mafft_c.log` | provenance |

Scripts uploaded to `/data/W/toki/`: `prep_sets.py`, `add_consensus.py`,
`run_mafft.sh`. Note KIT's default `python3` is **3.5** — use
`/usr/local/bin/python3.12`. No scipy. Reach KIT through the SSH tunnel on port
2223, never per-command direct connections (see the `server-connections` skill).

## 8. How to reproduce from nothing

```bash
# on KIT, via the tunnel
python3.12 /data/W/toki/prep_sets.py        # extract pools, emit 580 sets
python3.12 /data/W/toki/add_consensus.py    # prepend consensus, trim to 70 bp
bash       /data/W/toki/run_mafft.sh        # FFT-NS-2, 24-way parallel, ~4 min
# locally
python measure_c.py aln_c features_c.jsonl  # ~6 min, 14 workers
python profiles.py                          # positional tracks
python build_viewer_data.py                 # selection for the report
```

Corpus: 28 curated subfamilies from saq (GCA_004024925.1), ccr
(GCF_000260355.1), teu (GCA_964194135.1), dmo (GCA_051107935.1), drawn from
`results/assignment_full.tsv` filtered to `Status == "assigned"`. Note
`regions.by_subfam.bed` in every run directory is **broken** — every line lists
every subfamily. Use `assignment_full.tsv` or `all_hits.labeled.bed`.

Classes: `POS` 28, `SQ` 400 (5 families × n ∈ {25,50,100,200} × K=20),
`NEGRAND` 20, `NEGJITTER` 28, `NEGTRUNC5` 28, `NEGSPLICE` 20, `MIXED10/30` 56.

## 9. Residual circularity, stated plainly

Element anchors come from `sear`'s own alignment to the consensus, so the
centring is seed-derived. The consensus row re-derives the boundary
independently of the copies, but not independently of the consensus. Removing
this fully needs the spec's independent nhmmer re-search at a permissive
threshold. `NEGJITTER` was the stand-in and does not bite.

## 10. Next steps, in order

1. **Compare the boundary calls against Sergei's manual calls** on a few
   families. This closes spec §6 step 1, and the calls are now stable enough
   (SD 10.7 bp) to be worth comparing.
2. **Rebuild `NEGSPLICE`** with per-copy-varying breakpoints, and test Tier 3
   with a statistic that does not depend on imputation.
3. **Add natural negatives** — LINE candidates already exist at
   `/data/W/toki/Genomes/Mammalia/Eulipotyphla/teu/line/` on KIT; plus
   satellites and processed pseudogenes.
4. **Port the nesting check** to TRF overlap, following `eri/LOG.md`.
5. **Then fit a model.** Not before 3.

## 11. One open question for Sergei

`POS__ccr__a_ccr` is a curated real family that scores like a negative:
`cliff` 0.285, `cons_identity_med` 0.572, `tsd_frac` 0.03,
`res_asymmetry` −1.79. Either it is a genuinely degenerate old family and the
thresholds must accommodate it, or it is mis-curated. This needs his judgement,
not more computation.

---

## 12. Exploratory edge refinement (added same day, Sergei's suggestion)

The consensus anchor is stable but inherits whatever is wrong with the
consensus. `edges.py` treats the consensus edge as a starting guess and searches
offsets around it.

**Objective.** At a candidate edge, boundary quality is the step in mean
pairwise identity between copies across it: `mean pair_id inside window` minus
`mean pair_id outside window`, W = 25, scanned over d = -60..+60.

**Avoiding circularity.** The optimised step value is deliberately NOT used as a
discriminating variable - optimising a statistic then scoring on it inflates
negatives too. What is reported is the SHAPE of the landscape: `prominence`
(peak minus median over the scan) and `d_best`. A real family should have a
sharp optimum; junk should have a flat one whatever its best step happens to be.
That prediction held.

**Ground-truth validation** (`test_edges.py`). The consensus edge is deliberately
moved inward by 10/20/40 bp and the search asked to recover it. Left edge: error
within 5 bp in 19/20 cases, usually within 1 bp, and it returns 0 when nothing
is wrong. This is the check to re-run after any change to the objective.

**Landscape shape discriminates, as designed:**

| class | prominence L | prominence R |
|---|---|---|
| POS | 0.440 | 0.455 |
| MIXED10 | 0.397 | 0.396 |
| MIXED30 | 0.305 | 0.325 |
| NEGTRUNC5 | **0.086** | 0.463 |
| NEGRAND | **0.082** | **0.080** |

Random loci have no preferred edge on either side. 5'-truncated sets are flat on
the truncated side and sharp on the intact side - the asymmetry in *prominence*
is a second, independent truncation detector alongside `elem_len_cv` and
`res_asymmetry`.

**Two substantive findings about the SINEderella consensuses themselves:**

1. **The 5' ends are right.** All 28 curated families put `d_best_L` within
   +-1 bp of the consensus edge. The search confirms the existing boundary
   rather than moving it. This also settles a reading of the positional
   profiles: the identity bump at -20..-10 in the 5' flank is the A-rich
   insertion signature the spec predicts, NOT the element extending further.
2. **The 3' ends are systematically ~7 bp short.** 20 of 28 families want to
   extend outward, median -7 bp, many at exactly -7. Consistent with poly-A
   being trimmed during consensus construction while the copies retain it -
   the spec's own "residual tail features bleeding into the right flank".
   The 8 families that instead want to trim inward are mostly the ones with the
   weakest right-edge prominence (`ccr__a_ccr` 0.239, `saq__s2_38seqs` 0.284,
   `dmo__d1_16seqs` 0.279), i.e. their right edge is poorly defined either way.

**Known limits.**
- The scan reaches +-60 bp, bounded by the 70 bp flank in the alignment. A
  consensus wrong by more than that cannot be detected. Fix cheaply by carrying
  a longer *analysis* flank that never enters the alignment - flank length only
  hurts when MAFFT sees it.
- `width` (offsets within 90 % of the peak) came out at 6 for every class. It is
  degenerate as computed; `prominence` carries the information. Do not use it.
- The right edge is intrinsically harder because poly-A keeps pairwise identity
  high past the element, so the step is ill-defined there. The recovery test
  reflects this: right-edge errors are larger and mostly equal the real ~7-9 bp
  extension rather than being noise.

**Not yet done:** the refined edges are reported but not fed back into
`measure_c.py`. The obvious next move is a second pass - refine edges, then
re-measure everything against the corrected boundary - and to check whether the
+7 bp 3' extension improves `tsd_frac` and `polyA_score`, which it should if the
interpretation is right.

---

## 13. Structural feature annotation (`annotate.py`)

Drawn as a track under the positional profile, colour-coded, hover for sequence
and score. The important distinction, which the track encodes:

**Consensus features** are properties of the query consensus, not of the copy
set. Two sets searched with the same consensus get identical boxes, so these
cannot discriminate between such sets — consistent with the spec treating them
as one-sided evidence. Confirmed empirically: `NEGRAND__dmo__r00` shows the same
A and B boxes as the real families, because it inherits the same query.

**Copy features** are properties of the loci and do vary between sets.

### Detected offline and working

| feature | kind | result |
|---|---|---|
| `abox` | consensus | Pol III type-2 A box `TRGCNNARYGG`. **Exact match `TGGCGCAGCGG` at position 12–22** in both saq/s5 and teu/t2 — independent genomes, canonical position. |
| `bbox` | consensus | B box `GWTCRANNC`. Exact `GTTCGACGC` at 64–72, `GTTCGATCC` at 67–75. |
| `trna_region` | consensus | A-box start to B-box end, i.e. the tRNA-derived head. **Operational definition from the boxes, not a database match to a real tRNA.** |
| `simple_repeat` | consensus | homopolymer / di- / tri-nucleotide runs. saq/s5 carries (A)13 at 233–245, i.e. the poly-A inside the consensus. |
| `tail_repeat` | copies | the same scan on the majority base per offset in the copies' 3′ flank. saq/s5: (A)9 and (A)14 beyond the consensus end. |
| `internal_dup` | consensus | near-identical segment pairs, seed-and-extend on 12-mers, drawn joined by a dashed line. |
| `conserved_core` | consensus | the most conserved 30 bp. **Empirical — NOT the CORE-SINE domain**, which needs its reference consensus to assert. |
| `tsd` | copies | 0.83 of copies, median 12 bp in `POS__saq__s5_5seqs`; 0.12 and 6 bp in `NEGRAND__dmo__r00`. The random-locus hits are minimum-length spurious matches. |

The poly-A result corroborates §12 independently: the consensus ends inside a
poly-A run that continues into the copies' flank, which is exactly why the 3′
edge search wants ~7 bp more.

### Requires reference data, deliberately not attempted

Listed in each record's `not_attempted` field rather than faked:

- **tRNA / 5S rRNA identity** — needs a tRNA database (tRNAscan-SE, or ssearch36
  against a tRNA set). The `trna_region` above is the box-derived proxy only.
- **LINE 3′ end** — needs LINE sequences. Candidates already exist on KIT at
  `/data/W/toki/Genomes/Mammalia/Eulipotyphla/teu/line/LINE_candidates.fa`, so
  this is the cheapest of the three to add: ssearch36 the consensus against them
  and mark the matching span.
- **CORE-SINE domain** — needs the CORE reference consensus.

### Documentation discipline

Sergei asked that important points be written into the documentation as part of
each exchange rather than batched at the end. Findings go here (state) and into
`FINDINGS.md` (chronological), then get committed and pushed in the same turn.

---

## 14. Flank handling — variants tested, and why step 2 turned out unnecessary

Sergei raised this on `POS__ccr__g5_7seqs`: the flanks are genuinely
non-alignable, so MAFFT scatters them and their properties cannot be read off
the alignment. Measured on that set: **70 bp of flank per copy spread over 219
columns, 32 % occupancy, 3.1 columns per base**, against 1.74 for the element.

He proposed (1) de-gap and justify the flanks, then (2) re-align them with very
high gap-open penalties so gaps appear only when certain. Six strategies were
built and scored (`make_variants.py`).

**Scoring criterion.** The honest reference is flank identity measured ungapped,
walking outward from each copy's own edge — no aligner involved. A strategy is
good when identity measured *on its alignment* matches that reference rather
than exceeding it, while the element stays compact and well aligned.

| variant | width | flank cols/base | flank id (aligned) | ungapped truth | element id |
|---|---|---|---|---|---|
| **POS__ccr__g5_7seqs** | | | | | |
| v1_current | 867 | 2.89 | 0.351 | 0.344 | 0.929 |
| **v2_justify** | **610** | **1.05** | **0.344** | 0.344 | **0.929** |
| v3_op3 | 643 | 1.29 | 0.359 | 0.344 | 0.929 |
| v4_op10 | 630 | 1.19 | 0.357 | 0.344 | 0.929 |
| v5_op20 | 620 | 1.12 | 0.347 | 0.344 | 0.929 |
| v6_whole_op10 | 604 | 1.11 | 0.333 | 0.336 | 0.892 |
| **POS__saq__s5_5seqs** | | | | | |
| v1_current | 836 | 2.92 | 0.334 | 0.289 | 0.864 |
| **v2_justify** | **578** | **1.09** | **0.287** | 0.289 | **0.864** |
| v6_whole_op10 | 476 | 1.07 | 0.291 | 0.293 | 0.879 |
| **NEGRAND__dmo__r00** | | | | | |
| v1_current | 1783 | 7.12 | 0.346 | 0.270 | 0.331 |
| **v2_justify** | **1379** | **4.18** | **0.263** | 0.270 | 0.331 |
| v6_whole_op10 | 617 | 2.55 | 0.264 | 0.264 | 0.265 |

**Conclusion: v2 (de-gap and justify, no re-alignment) wins on every criterion.**
It reproduces the ungapped truth almost exactly (0.344 vs 0.344; 0.287 vs 0.289;
0.263 vs 0.270), compacts the display by ~30 %, leaves the element alignment
byte-identical, and has no parameters to tune.

**Step 2 is unnecessary, and the data show why.** Re-aligning the de-gapped
flanks (v3–v5) is strictly worse than not re-aligning them: it puts columns back
(1.29 → 1.19 → 1.12 as the penalty rises) and slightly re-inflates identity
(0.359 → 0.357 → 0.347 against a truth of 0.344). The trend with increasing gap
penalty points straight at v2 — the limit of "very high gap penalty" *is* pure
justification. There is nothing homologous in the flanks to align, so any
alignment there manufactures similarity.

**How much the old geometry was lying.** The inflation is worst exactly where it
matters most — the null. `NEGRAND` read 0.346 aligned against 0.270 ungapped,
+0.076. Real families were inflated less (+0.007 and +0.045). So aligned-flank
statistics systematically flatter the negatives, narrowing the very gap the
discriminator depends on.

**Note.** `measure_c.py` already computes flank identity ungapped and
edge-anchored, so no published statistic is affected by this. It is a display
and re-alignment problem, and it matters for manual checking — which is the
point of the viewer.

**v6 (re-align everything at high gap penalty) is not recommended but is not
uninteresting:** it improved the element alignment for saq/s5 (identity
0.864 → 0.879, width 836 → 476) while degrading it for ccr/g5 (0.929 → 0.892).
Inconsistent across sets; would need a proper test on the whole corpus before
adopting. It is a question about the *element* alignment, separate from flanks.

Files: `variants/<set>__<variant>.aln.fa` in the repo, viewable in MSA-viewer.

### 14a. Follow-up — the 5′ insertion-site motif, and why it settles v2 vs v4

Sergei spotted that v4 appears to highlight the first two T residues of a
TTAAAA-like motif in the 5′ flank and asked whether that biological signal
should decide the choice. It does decide it — for v2.

**The motif is real and specific.** Majority base per offset going 5′ from the
element edge, four real families and one null:

```
POS ccr/g5   ATTTAAAAAAAAAAAAAATTTTTTTAATTAAAAAATAAAA
             3333443435667876544434333343333333333333   <- majority fraction x10
POS saq/s5   ATATAAAAATAAAAAAATATATAATAAAAATAAAAAATAT
             3234344343566764433434334334343333343333
POS dmo/d4   TTTTAAAAAAAAAAATTTTTTTAAAATTTAAAAAATAAAA
             3434343346677744433333333333333343333343
NEGRAND dmo  TTAATCTTTTTTTTTTTTTTCTATTTTTTTATTTTTCATT
             3333333434333333333323233433322343332333   <- flat, no peak
```

Every real family has a conservation peak at offsets ~10–17 reaching 70–83 %
majority; the null is flat at ~30 %. Read in genomic 5′→3′ order the ccr/g5
flank is `…TTTTTTT AAAAAAAAAAAAAA TTTA |element` — a T-tract followed by an
A-tract. This is almost certainly the L1 endonuclease target site (L1 EN nicks
at 5′-TT|AAAA-3′, and SINEs mobilised in trans by L1 inherit that preference),
and it is what spec §1.1 calls "an A-rich zone ~10 bp upstream of the core".
**It is now confirmed empirically and shown to be absent from random loci.**

**v4 does not resolve this motif better — it resolves it worse.** Column
composition for the 24 columns before the element:

| | v2_justify | v4_op10 |
|---|---|---|
| gaps in the element-adjacent columns | **0 %** | rising to **99 %** |
| A-tract peak purity | **83 %** | 80 % |
| T-tract peak purity | **47 %** | 33 % |
| proximal 20 columns, mean identity | **0.383** | 0.331 |
| distal flank, mean identity | **0.269** | 0.300 |

v4 inserts a gap block immediately against the element — the last eight columns
are 88–99 % gaps — which destroys the edge anchoring the motif's position is
defined against. It simultaneously raises identity in the *distal* flank, where
there is nothing real to find. So it weakens the true proximal signal and
manufactures a false distal one: exactly the wrong trade. What looks like a
cleanly highlighted TT column in the viewer is the aligner gathering T residues
from different distances into one column.

**A tested alternative that also failed.** If the motif sat at the outer end of
the TSD, its distance from the element edge would vary with TSD length (median
14 bp here, IQR 11–15, found in 95/100 copies) and an aligner would beat pure
justification. Anchoring each copy on its own outer TSD boundary instead of the
element edge was tried and is worse: peak purity 0.55 against 0.83. So the motif
is anchored to the **element edge**, not to the TSD — which is precisely the
anchor v2 preserves exactly and v4 breaks.

**Recommendation: v2.** It shows the biology more sharply, not less, and it is
the only variant that keeps the element-edge anchoring intact.

**Worth pursuing separately (biology, not display):** the T-tract/A-tract
insertion signature is a per-copy feature and a candidate discriminating
variable in its own right — real families have it, random loci do not. It should
be added to `annotate.py` as a scored motif and to `measure_c.py` as a one-sided
bonus, replacing the current crude `arich_score`.

---

## 15. v2 adopted corpus-wide; report restructured (Sergei confirmed)

`justify_all.py` applied the v2 strategy to all 580 alignments → `aln_v2/`,
which is now the current corpus. `aln_c/` is kept as the pre-justification
source (the element alignment is byte-identical between them; only the flanks
differ). Everything downstream was re-run on `aln_v2`:
`features_v2c.jsonl`, `profiles.json`, `annotations.json`, `commentary.json`.
The repo's `alignments/` directory now holds the v2 files, so all MSA-viewer
links serve the justified alignments.

The report is now a **vertical list** — one panel per alignment, each with the
positional profile, the structural feature track in distinct hues, the measured
parameters, a link into MSA-viewer for all copies, the alignment itself
collapsed, and per-set notes on what works and what does not.

`commentary.py` generates those notes **from the measured values**, so they
cannot drift away from the data. One bug found and fixed while checking the
output: the "3′ edge wants to move outward → poly-A" note fired on the
random-locus set, where there is no poly-A and no element. It is now gated on
`cliff > 0.3`. Generated prose still has to be read against the data it
describes.

## 16. SCOPE — everything here is ONE SINE family type

Sergei's caution, and it qualifies every biological claim in this document:
**the entire corpus is Tal, a tRNA-derived SINE, in four talpid genomes.**
Biology differs between SINE families, so the following are family-specific
findings and must not be stated as general properties of SINEs:

- the Pol III A box at 12–22 and B box at 64–75 — positions are family-specific,
  and a 5S-derived or non-tRNA-derived SINE differs entirely;
- the T-tract/A-tract L1 endonuclease insertion signature 10–17 bp upstream —
  families mobilised by a different LINE partner have a different target site,
  and the spec already notes that many SINEs lack TSDs altogether;
- TSD frequency ~0.74 and median length 9–14 bp;
- the ~7 bp 3′ consensus under-call, which is an artifact of how these
  particular consensuses were built;
- element length ~250–280 bp.

What should generalise, because it is geometry rather than biology: the
consensus-row boundary anchor, short flanks, ungapped edge-anchored flank
statistics, and de-gap-and-justify display.

**The obvious next test is a second family.** The Tal repo already carries
Erinaceidae (e1/e2, hedgehog), scorpion and Timema sets with their own
consensuses and run directories. Running the same pipeline on Erinaceidae would
show which of the numbers above are properties of SINEs and which are properties
of Tal. Until that is done, no threshold derived here should be applied to
another family.
