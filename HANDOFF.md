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

### 15a. Chart colour and layout fixes

Three problems Sergei reported on the vertical list, all fixed:

1. **Overlapping feature labels.** Lanes were packed by span width alone, but a
   9 bp A box carries a ~30 px label, so labels from adjacent narrow features
   collided. Lanes are now packed by the space the *label* occupies, and a label
   that would run past the right margin is hung to the left of its block instead.
2. **Track toggles were invisible.** They were an unlabelled inline row. There is
   now a sticky control bar at the top of the section with two labelled groups —
   profile tracks and structural features — with colour swatches, and unchecked
   items dimmed. Feature types are individually toggleable, which they were not
   before.
3. **Colours were not distinct enough.** Replaced ad-hoc hues with the validated
   8-slot categorical palette (blue, orange, aqua, yellow, magenta, green,
   violet, red), stepped separately for light and dark surfaces rather than
   flipped. Checked with the palette validator rather than by eye:

   - light: all pass. Worst adjacent CVD ΔE 9.1 (target ≥ 8), worst adjacent
     normal-vision ΔE 19.6 (floor ≥ 15). Contrast WARN on three hues, which
     obligates visible labels — every feature block carries one, so that relief
     is satisfied.
   - dark: all pass, including contrast ≥ 3:1 on all eight.

   Feature slots are fixed per type and never cycled, so a feature keeps its
   colour across panels.

---

## 17. MIXED is a verdict, not a negative — the refinement loop (spec §4.3)

Sergei, on `MIXED30__saq__s3_43seqs`: it should come out clearly positive,
because a real SINE is recoverable from that alignment. He is right, and the
error was in framing, not in the statistics.

Spec §1 asks for **two** failure modes with separate verdicts — "not a family at
all" and "a mixture: real copies plus contaminants" — and §4.3's required output
for the second is the family, refined bounds, and a *cleaned copy list*. Builds 1
and 2 instead scored MIXED as a negative class to be separated from POS, which
answers the wrong question. The statistics were already saying the right thing:
that set has `cliff` 0.559 (**above** the POS median of 0.563 by a hair, i.e.
squarely in the real range) and `frac_supported` 0.740 against a planted truth of
70 % real. The family was never in doubt; only the reporting was.

`prune.py` implements the loop: measure per-copy identity to the known consensus,
drop the copies that do not support it, re-measure, repeat to three passes.

### Which cut rule — and a dangerous one to avoid

Two rules were compared against ground truth (planted contaminants carry
`contam` in the name, so precision and recall are measurable, not assumed):

| rule | class | removed | precision | recall | identity of kept |
|---|---|---|---|---|---|
| **mad** (median − 3·MAD) | POS | **0.2** | — | — | 0.832 |
| | MIXED10 | 9.2 | 0.973 | 0.893 | 0.850 |
| | MIXED30 | 26.9 | **0.993** | 0.892 | 0.863 |
| | NEGRAND | **0.1** | — | — | **0.178** |
| gap (largest gap in sorted identity) | POS | 6.0 | — | — | 0.836 |
| | NEGRAND | **89.4** | — | — | **0.803** |

**Use the MAD rule.** The gap rule fails in the worst possible way: on random
loci it removes 89 of 100 copies and leaves a residue with median identity 0.803,
i.e. **it manufactures an apparent family out of noise**. Any pruning rule that
splits on the data's own largest gap will do this, because a large enough random
set always contains a few loci that resemble the query. The MAD rule leaves
random loci alone (0.1 copies removed, identity stays 0.178) and leaves clean
families alone (0.2 copies removed). Not manufacturing families from noise is the
property to protect; recall is secondary.

### Result — pruning recovers the family

| class | cliff | cons_identity | frac_supported | elem_len_cv | tsd_frac |
|---|---|---|---|---|---|
| POS (control) | 0.563 → 0.563 | 0.870 → 0.870 | 1.000 → 1.000 | 0.052 → 0.052 | 0.740 → 0.740 |
| MIXED10 | 0.548 → 0.551 | 0.855 → 0.865 | 0.910 → **1.000** | 0.107 → **0.059** | 0.680 → **0.741** |
| MIXED30 | 0.542 → 0.555 | 0.846 → 0.871 | 0.730 → **1.000** | 0.144 → **0.069** | 0.545 → **0.729** |
| NEGRAND | −0.089 → −0.089 | 0.184 → 0.184 | 0.110 → 0.110 | 0.180 → 0.180 | 0.070 → 0.070 |

After one pass both MIXED classes are statistically indistinguishable from clean
families, POS is untouched, and random loci are untouched — they stay junk.

The specific set: 100 copies → 74 kept; `frac_supported` 0.740 → 1.000,
`elem_len_cv` 0.137 → 0.065, `tsd_frac` 0.60 → 0.78, against 0.041 and 0.85 for
the never-contaminated original. It converges on the right answer.

### Consequences for the earlier reporting

- **The AUC tables in §3 treat MIXED10/MIXED30 as negatives. That framing is
  wrong** and should be read as "how well contamination is *detected*", not as
  "how well a bad set is rejected". The correct pipeline output for these sets is
  SINE + cleaned copy list, and `frac_supported` is the detector.
- §7's readiness table listed the pruning loop as "not built". It is built and
  validated now.

### Honest limits

- **Recall is 0.89**, so roughly one contaminant in ten survives — about 3 of 30
  in a MIXED30 set. Precision 0.99 means almost nothing real is thrown away,
  which is the right trade for a cleaned copy list.
- **Rows are dropped without re-aligning.** Spec §4.3 says prune, realign,
  rescore. The element alignment is anchored on the unchanged consensus row so
  re-measuring on the remaining rows is defensible, but the fuller version needs
  a MAFFT round on KIT per pass.
- Contaminants here are *random genomic loci*, the easy case. A contaminant that
  is a diverged member of a related family will sit much closer to the cut and
  has not been tested.

---

## 18. Weighted verdict with named sub-cases (`verdict.py`)

Three changes Sergei asked for, all built.

### A weighted score, not a binary call

Four evidence groups, each 0–1, combined **geometrically** with stated weights:
element 0.45, homogeneity 0.25, uniqueness 0.20, and insertion 0.10 applied as a
one-sided bonus (it can lift a borderline score, never reject — as spec §1.2
requires). Geometric rather than a weighted sum, because a sum lets one strong
group mask a fatal weakness: a set with no element is not a family however unique
its flanks are. Every component is shown in the report so a disagreement can be
traced to the variable that caused it.

**The weights are stated, not fitted.** Fitting now would fit the synthetic
negatives.

| class | score | element | homogeneity | uniqueness | insertion |
|---|---|---|---|---|---|
| POS | 96.7 | 0.93 | 0.99 | 0.92 | 0.78 |
| MIXED10 | 97.5 | 0.88 | 0.97 | 0.94 | 0.86 |
| MIXED30 | 87.7 | 0.68 | 0.95 | 0.96 | 0.89 |
| MIXSUBFAM | 94.5 | 0.93 | 0.91 | 0.94 | — |
| NEGTRUNC5 | 69.0 | 0.71 | 0.60 | 0.89 | 0.00 |
| NEGRAND | **0.1** | 0.00 | 0.62 | 1.00 | 0.18 |

**A circularity caught while building this.** The first version scored random
loci at 66.8, because the cliff was measured over the *core* copies — the ones
selected for matching the consensus. Selecting copies by identity and then
measuring their identity is the same trap that made the gap-based pruning rule
manufacture families out of noise (§17). The cliff is now measured over all
copies, and the mean rather than the median, since with 30 % contamination the
median is still a real copy.

### Sub-cases, because "negative" is not one thing

| flag | meaning | POS | MIXED30 | MIXSUBFAM | NEGRAND |
|---|---|---|---|---|---|
| `NO_ELEMENT` | nothing above background | 0/28 | 0/28 | 0/8 | **20/20** |
| `CONTAMINATED` | element present, some copies do not support it | 2/28 | **28/28** | 0/8 | 16/20 |
| `SUBFAMILY_STRUCTURE` | supported copies split into structurally different groups | **0/28** | 0/28 | **6/8** | 0/8 |
| `NESTED_COPIES` | copies share flanking sequence — inside another repeat, a duplication, or a satellite | 5/28 | 3/28 | 2/8 | 0/20 |
| `RECOVERABLE_CORE` | poor set, but ≥ 20 copies look genuine | — | — | — | — |

`RECOVERABLE_CORE` is Sergei's point that 20–30 real-looking copies inside a
heavily polluted candidate are worth reporting, because that is where the next
clean family comes from. **`n_core` — copies that both support the consensus and
sit in unique genomic sequence — is reported for every set regardless of
verdict.**

### The missing negative, built and used to calibrate

`SUBFAMILY_STRUCTURE` had nothing to fire on, so `kit/mixsubfam.py` builds it:
8 sets of 50 copies from subfamily A plus 50 from subfamily B of the same
species, searched with A's consensus only. Detector: copy × copy identity inside
the element, split on the leading eigenvector, scored as within-group minus
between-group identity.

Calibrated against the control rather than guessed — clean families reach a gap
of at most 0.018 (median 0.007); mixtures run 0.011–0.126 (median 0.071). The
threshold is 0.03: **6/8 mixtures caught, 0/28 clean families falsely flagged.**
The two misses (`ccr g1+g6` at 0.011, `dmo d3+d5` at 0.022) are pairs whose
subfamilies are genuinely close in structure — this statistic cannot separate
subfamilies that differ only by age.

### Motif profiles instead of solid boxes

Sergei's third point: a solid box says only where the best match is. The A box,
B box and tRNA-head similarity are now scored **at every position** and drawn as
lines, so a duplicated box, a decayed second copy of the head, or a box barely
better than background are all visible. On `POS__saq__s5_5seqs` the A box scores
1.00 at position 12 with the next-best peak at 0.64, and the tRNA head has no
second peak above 0.37 — i.e. no internal duplication of the head in that family.

The boxes are kept as well: the profile shows the landscape, the box marks the
call. Both are individually toggleable.

---

## 19. NEGSPLICE was never a negative; the real chimera and the real mosaic

Sergei, on `NEGSPLICE__saq_s5_5seqs__saq_s8_225seqs__02`: "very clear SINE,
almost perfect — reconsider your criteria." The criteria were right; the class
label was wrong.

**Every one of the 20 NEGSPLICE sets joins two Tal SUBFAMILIES.** All consensuses
in this corpus are Tal, so splicing s5 to s8 produces a within-family
recombinant — a real biological entity that *should* score as a SINE. Scoring it
~100 is the correct answer. The class has been reporting a success as a failure
since build 1, and every earlier statement that "chimeras are not caught" was
measuring the wrong thing.

`kit/negchim.py` builds the two cases it should have been:

- **NEGCHIM** — the 3′ half of each copy replaced with random genomic sequence,
  breakpoint varying per copy. Half element, half non-element.
- **NEGMOSAIC** — two subfamilies with the breakpoint varying **per copy**, so
  different copies genuinely carry different blocks. This is the spec's Tier-3
  kaleidoscope, which the old fixed-midpoint splice never produced.

| class | score | n_core | flags |
|---|---|---|---|
| POS | 96.7 | 89 | — |
| NEGSPLICE (recombinant, not a negative) | 96.7 | 47 | — |
| **NEGMOSAIC** | **97.6** | 91 | none |
| **NEGCHIM** | **19.1** | 14 | NO_ELEMENT 12/16 |
| NEGRAND | 0.1 | 10 | NO_ELEMENT 20/20 |

**Good news: a genuine chimera is rejected.** NEGCHIM scores 19.1 with
NO_ELEMENT on 12 of 16. The criteria do catch half-element/half-junk; they were
simply never shown one.

**Bad news, and this is the clean negative result: a real mosaic is NOT
detected.** NEGMOSAIC scores 97.6, indistinguishable from a clean family, and
neither statistic sees it:

| class | rank1_excess | rank2_frac | elem_len_cv |
|---|---|---|---|
| POS | 0.0017 | 0.0022 | 0.052 |
| **NEGMOSAIC** | **0.0017** | 0.0032 | 0.053 |
| NEGCHIM | 0.0071 | 0.0135 | 0.094 |
| NEGRAND | 0.1352 | 0.0246 | 0.180 |

`rank1_excess` on a real mosaic is **identical to a clean family to four decimal
places**. So spec §2's central Tier-3 idea — SVD rank as the mosaicism detector —
now has a properly constructed test and **fails it**. Earlier builds recorded
Tier 3 as "untested because the negative was mis-designed"; it is now tested, and
the statistic does not work as specified.

A likely reason, worth checking before discarding the idea: every NEGMOSAIC copy
is still a valid Tal element, so the copy × window identity matrix stays nearly
rank 1 — a copy that is 40 % subfamily A and 60 % B matches the A consensus at a
level that still looks like "an older copy". Separability may need the matrix
built against *two* candidate consensuses rather than one, or a phylogenetic
incongruence test (the PHI/GARD route the spec also lists) rather than SVD.

**Also note NEGMOSAIC is arguably not a negative either** — subfamily
recombinants occur. The right output for it is probably `SUBFAMILY_STRUCTURE`
with a note that the breakpoint varies, not rejection. That distinction has not
been built.

## 20. Motif tracks were plotting the chance baseline

Sergei: "in all cases what I see in small tRNA/boxes tracks they look as pure
noise." Correct, and measurably so. Raw fraction-of-positions-matched for a
degenerate motif containing `N` and `R/Y` sits at **median 0.45, IQR 0.33–0.56**
across the whole element — the eye sees that band, and the one real hit is a
single point at 1.0.

Now scored as a **z against the same motif scanned over shuffled versions of the
same consensus**, which puts background at 0:

| track | before (raw fraction) | after (z vs shuffled null) |
|---|---|---|
| A box | median 0.46, IQR 0.36–0.55, max 1.00 | median −0.01, IQR −0.74–0.73, **max 4.4** |
| B box | median 0.44, IQR 0.33–0.56, max 1.00 | median −0.13, IQR −0.96–0.70, **max 4.0** |
| tRNA-head self-similarity | median 0.29, max 1.00 (trivial self-match) | median −0.24, **max 2.4**, self-match region masked |

Only 0.4 % of positions now exceed |z| > 3. Two things this makes readable:

- the boxes are real but modest — z ≈ 4 for an 11 bp degenerate motif is a
  genuine signal, not an overwhelming one, and it should not be oversold;
- **no internal duplication of the tRNA head** in these families: the
  self-similarity track peaks at z ≈ 2.4 off-target, which is not significant.
  That is now a stateable negative result rather than an unreadable band.

---

## 21. Chart layout bug, and gradient rows instead of lines

**The bug.** §20 changed the motif tracks from 0–1 fractions to z-scores, but the
chart still mapped them with `(1 - v) * height` as if they were fractions. A z of
4.4 therefore drew *above* the strip and a z of −2.7 drew ~126 px *below* it,
straight across the feature track — which is what Sergei saw. Fixed by scaling
the strip to a stated z range (−3 … 6) with clamping, plus z=0 and z=3 reference
lines. A regression check now runs on every build: every `<path>` in every chart
is measured against its own viewBox, and any that escapes is reported.

**Gradient rows.** On Sergei's suggestion, the tracks whose *shape* matters more
than their exact value are now heat rows rather than lines — three overlapping
noisy z-score lines were much harder to read than three stacked bands:

- **A box / B box / tRNA-repeat**: sequential ramp, one hue each, darker = higher
  z (1.5 → 5). A hit is a dark tick; the background is blank rather than a noisy
  line at zero.
- **A+T composition**: diverging ramp — blue GC-rich, orange AT-rich, neutral
  grey at the genomic background of 0.62 — smoothed over 9 bp, because at 1 bp
  resolution it is speckle and the shape is the point.

Lines are kept for pairwise identity, identity to consensus and coverage, where
the actual value is read off the axis.

## 22. Why A+T behaves differently inside the element — it is a Pol III fingerprint

Sergei asked. Measured across the real families:

| region | A+T |
|---|---|
| 5′ flank | 0.642 |
| **element 0–15 (the A box)** | **0.223** |
| element 15–80 (tRNA head) | 0.396 |
| element 80 … L−30 | 0.389 |
| **element last 30** | **0.755** |
| 3′ flank | 0.629 |
| random loci: flank vs middle | 0.595 vs 0.545 — flat |

The element is **GC-rich because it is tRNA-derived**. tRNA genes need stable
stem pairing, so they are GC-rich, and the A box `TGGCGCAGCGG` is the most
GC-rich part of the whole element at 78 % GC. The 3′ end swings the opposite way
to 0.755 A+T — that is the poly-A tail. The flanks sit at mammalian genomic
background, ~0.63, and random loci show no structure at all (0.595 → 0.545).

So the A+T track is not noise: the *shape* — background, sharp GC dip at the 5′
end, moderate through the body, A-rich spike at the 3′ end, background again — is
a structural signature of a tRNA-derived Pol III element. It is also a candidate
discriminating variable that has not been formalised, and unlike the current
`arich_score` it would use the whole profile rather than one window.

Scope caveat from §16 applies with force here: these numbers are the fingerprint
of a **tRNA-derived** SINE. A 5S-derived or non-tRNA family will have a different
composition profile, and the A-box position that anchors it will not be at 12–22.

---

## 23. How the boxes are searched and drawn, and what the grey in A+T means

Both were fair questions about things the display was getting wrong.

### Searching

The A box pattern is `TRGCNNARYGG` — 11 positions of IUPAC codes, so `R` = A or
G, `N` = anything. At **every** position in the consensus, count how many of the
11 positions match their allowed set. That raw count is then converted to a
z-score against a null built by shuffling the consensus 120 times and rescanning
the same motif, so composition is controlled for. The track is therefore
continuous: one value per position, not a set of hits.

### Why they looked like discrete rectangles — a display bug, not the data

The heat row was blanking every cell below z = 1.5. With background at z ≈ 0 ±
0.8, about **7 % of positions clear 1.5 by chance**, so the row drew the real hit
plus ~15 chance hits at comparable weight, reading as scattered ticks. The floor
was hiding the continuity and giving noise the same visual weight as signal.

Now a continuous ramp from z = 0 to z = 5 with a gamma of 2.2, so chance-level
positions are faint, the real box is solid, and the row reads as a landscape with
one peak. Nothing about the underlying scan changed.

### The grey in A+T

A+T uses a **diverging** ramp, which has a meaningful midpoint: blue for GC-rich,
orange for AT-rich, and **grey where the composition sits at the genomic
background of 0.62** — i.e. grey means "ordinary DNA, no compositional bias".
That is why the flanks are largely grey and the element is strongly blue: see
§22. Cells within 8 % of the midpoint are painted neutral, and the row is
smoothed over 9 bp because at 1 bp resolution it is speckle.

## 24. Pol III terminator (Sergei's bonus feature)

Added to `annotate.py`. RNA Pol III terminates at a run of T on the non-template
strand:

- **strong** — `TTTT` or longer, the canonical signal
- **moderate** — `TCTTT` / `TGTTT` / `TATTT`, a T-run interrupted by one non-T,
  which still terminates but less efficiently, so read-through is likelier and
  the element can acquire a 3′ extension

Searched both in the consensus and, more usefully, **per copy** in each copy's
own 40 bp of downstream sequence, since the functional terminator is the first
one past the element end.

| set | strong | moderate | none | median distance |
|---|---|---|---|---|
| POS saq/s5 | 0.23 | 0.11 | 0.66 | 13 bp |
| POS ccr/g5 | 0.25 | 0.08 | 0.67 | 24 bp |
| POS teu/t2 | 0.23 | 0.11 | 0.66 | 23 bp |
| POS dmo/d4 | 0.19 | 0.04 | 0.77 | 24 bp |
| NEGRAND | 0.13 | 0.10 | 0.77 | 14 bp |
| NEGTRUNC5 | 0.17 | 0.10 | 0.73 | 20 bp |

**Report this honestly: it discriminates only weakly.** Real families carry a
strong terminator downstream in 19–25 % of copies against 13 % for random loci —
a real difference but a small one, nothing like `cliff` or `frac_supported`. Two
things follow:

1. Two thirds of copies have **no** terminator within 40 bp. That is a
   biological observation worth checking rather than a detector failure —
   read-through is expected to be common, and these are old degraded copies.
2. The search window starts at the consensus 3′ end, which §12 showed is
   systematically ~7 bp short. The window origin is therefore slightly wrong,
   and re-running it against the refined 3′ edge is the obvious next test.

### A palette consequence worth recording

The validated categorical palette is capped at eight hues, and the feature track
was already using all eight. Rather than invent a ninth, `trna_region` — which
is *defined* by the A and B boxes rather than being independent evidence — is now
drawn as a faint bracket in the A-box hue, freeing aqua for the terminator. Any
further feature type (LINE 3′ end, CORE domain, 5S similarity) must either reuse
a hue within a category or the track must split into rows by category:
promoter elements / repeats / insertion evidence / termination.

---

## 25. What mosaicism actually means here, and how it is estimated

Sergei corrected a misunderstanding that had been running through the whole
project. Record it before anything else in this section.

**What the code assumed:** two-parent recombination — copies built from
subfamily A and subfamily B with a breakpoint.

**What Sergei means:**

```
consensus   a b c d
copy 1      a b b d      slot 3 carries segment b
copy 2      a d c d      slot 2 carries segment d
```

The element's **own** segments reshuffled, duplicated, or **missing** between
slots, with the affected slot varying from copy to copy. Every copy is still
made entirely of this element's parts. He notes this is common in non-SINE
repeats, and that deletions — a long gap where other copies have sequence — are
as frequent as rearrangements.

### The three estimators, and what each can and cannot see

**1. rank-1 residual (spec §2, `rank1_excess`).** Copy × window identity should
factor as age × positional conservation. Cannot see any of this: a scrambled
copy still matches the consensus at roughly the right overall level.
Measured identical to a clean family (0.0017 both).

**2. partition congruence (tests Sergei's k-mer idea, `mosaic_kmer.py`).** In
each window, split the copies in two; ask whether the split is the same from
window to window, by adjusted Rand index. Congruent = one family or clean
subfamilies; incongruent = block-swapped.

| class | congruence, adjacent windows | all pairs |
|---|---|---|
| POS | 0.155 | 0.034 |
| **MIXSUBFAM** | **0.498** | **0.368** |
| NEGMOSAIC (2-parent) | 0.194 | 0.042 |
| NEGRAND | 0.478 | 0.350 |

It is a **good subfamily-structure detector** — 10× separation from clean
families, independent of the eigenvector-gap detector already in `verdict.py` —
but it does not detect the mosaic, and it has an artifact on random loci
(congruence driven by coverage, not sequence) that needs a guard.

**The k-mer SPECTRUM route fails, as expected.** Per-copy k-mer composition
correlations for k = 4/6/8/10: POS 0.681/0.405/0.272/0.200 against NEGMOSAIC
0.689/0.412/0.276/0.200 — identical. The spectrum discards position, and
mosaicism is entirely positional. k-mers are useful here as a way to ask *which
consensus segment does this slot resemble*, not as a composition summary.

**3. segment mapping (`segmap.py`) — the one that matches Sergei's definition.**
Split the consensus into segments; for each copy and each slot, ask which segment
that slot best matches. Clean copies map slot j → segment j.

**And a methodological finding that matters:** the aligner hides this.

| set | diag_frac on aligned columns | on unaligned copies |
|---|---|---|
| SIMCLEAN | 0.997 | 0.902 |
| SIMSCRAM swap f = 0.5 | 0.995 | 0.784 |
| SIMSCRAM swap f = 1.0 | 0.992 | **0.505** |
| POS (real family) | 0.999 | 0.808 |

MAFFT aligns a duplicated or swapped segment **back to the segment it came
from**, converting a slot substitution into an indel pair and destroying the
positional evidence. Any detector for this must work on unaligned sequence.
Note the real-family baseline is 0.808 rather than ~0.9, because indels shift
proportional slot boundaries — the threshold has to be set against real data,
not against the simulation.

### Deletions, which the alignment does keep

Unlike a swap, a missing segment survives as a shared gap. Two statistics,
calibrated against real families:

| class | copies with a long internal gap (≥20 bp) | position concentration |
|---|---|---|
| SIMCLEAN | 0.013 | 0.10 |
| **POS** | **0.046** | 0.20 |
| SIMDEL, same block in every affected copy | 0.26–0.54 | **0.65–0.94** |
| SIMDEL, a different block per copy | 0.15–0.46 | 0.26 |
| NEGCHIM | 0.596 | 0.44 |
| NEGRAND | 0.734 | 0.31 |

Clean families essentially do not carry long internal gaps (4.6 % of copies).
**Concentration separates the two modes**: one shared deletion concentrates at a
position (0.65–0.94), while per-copy deletions scatter (0.26). On
`SIMDEL__one__f050` the block was recovered exactly — position 110, length ~50,
in 51 of 100 copies, against a planted deletion of 0.2 L at 0.45 L.

This is spec §4.1 G's `big_indel_frac`, which had never been implemented.

## 26. Synthetic corpus (`kit/simulate.py`)

Every negative until now was carved out of real data, and NEGSPLICE showed how
badly that can go — it was never a negative at all. The generator gives exact
control and exact labels. Real genomic background from the same species, so
composition and repeat content are not idealised away.

| grid | parameter swept | sets |
|---|---|---|
| `SIMCLEAN` | age 0.05 → 0.30 | 4 |
| `SIMMOSAIC` | 2-parent recombinant fraction × breakpoint spread | 6 |
| `SIMSUBFAM` | subfamily divergence 0.05 / 0.12 / 0.25 | 3 |
| `SIMNEST` | fraction sharing a host flank | 2 |
| `SIMTRUNC` | 5′ truncation rate | 2 |
| `SIMSCRAM` | segment swap / duplication, fraction 0.2–1.0 | 6 |
| `SIMDEL` | segment deletion, same block or per-copy | 4 |

Already earning its keep: `SIMSUBFAM` shows the subfamily detector needs
divergence ≥ 0.12 to fire (d = 0.05 scores 100 and is missed), and `SIMTRUNC`
reproduces the truncation signature at a known rate.

## 27. Report UI

- Feature blocks reduced from 13 px to 8 px — they were far heavier than the
  information they carry.
- **Tooltips on everything.** Transparent hover columns across the plot report
  the position and every visible track's value at that point; each heat cell
  reports its own value and units; each feature block reports its sequence,
  score and note.
- **Double-click any panel** for a detail dialog: the consensus as nucleotides
  in 60-base lines with coordinates, coloured where a structural feature covers
  it, plus the full feature table and the scalar summary.

---

## 28. Gradient sweep: 1400 sets, one parameter each, and what it exposed

Sergei asked for a synthetic test of every feature on a gradient. Built
(`kit/gradients.py`): **14 grids × 100 sets**, each sweeping ONE parameter
continuously with everything else at a realistic default, on real saq genomic
background. Aligned on KIT in about two minutes.

| grid | parameter | range |
|---|---|---|
| AGE | substitution rate | 0.02 → 0.35 |
| NCOPIES | copies per set | 20 → 200 |
| TRUNC | 5′ truncation rate | 0 → 0.70 |
| CONTAM | random-locus fraction | 0 → 0.60 |
| SUBDIV | subfamily divergence | 0 → 0.30 |
| SUBFRAC | fraction from the second subfamily | 0 → 0.50 |
| SCRAM | segment swap fraction | 0 → 1.0 |
| DUP | segment duplication fraction | 0 → 1.0 |
| DELSHARED | fraction with the same block deleted | 0 → 0.80 |
| DELCOPY | fraction with a per-copy block deleted | 0 → 0.80 |
| NEST | fraction sharing a host flank | 0 → 0.80 |
| TSD | TSD length | 0 → 22 bp |
| POLYA | poly-A length | 0 → 40 bp |
| RECOMB | 2-parent recombinant fraction | 0 → 1.0 |

`gradient_analysis.py` computes Spearman rho between each swept parameter and
every statistic. This asks a sharper question than "does it separate": is the
statistic **sensitive** (monotone in the thing it claims to measure) and
**specific** (flat for everything else)?

### The specific diagnostics — safe to build verdicts on

| statistic | responds to | strongest |
|---|---|---|
| `res_asymmetry`, `resL_iqr` | **1** | TRUNC +0.79 |
| `frac_supported` | 3, dominated by one | **CONTAM −0.96** |
| `gap_concentration` | **1** | DELSHARED +0.74 |
| `long_gap_frac` | 6 | DELCOPY +0.98, DELSHARED +0.98 |
| `seg_diag` | 5 | **SCRAM −0.88, DUP −0.85** |
| `tsd_len_med` | 2 | TSD +0.84 |

### Four confounds this exposed, all of which invalidate current claims

**1. `rank1_frac` is a measure of AGE, not mosaicism.** rho = **−1.00 with
AGE** — a perfect rank correlation. `rank1_excess` +0.99, `rank1_null` −1.00.
Earlier sections recorded that the spec's Tier-3 SVD statistic "fails to detect
a mosaic". The sweep says why: it is an age meter. Any apparent mosaicism signal
from it is age leaking through. §19's conclusion stands but the reason is now
established, and the statistic should be dropped rather than repaired.

**2. `sub_gap` — my own subfamily detector — is confounded by age and sample
size.** AGE +0.99, NCOPIES **−0.97**, against SUBFRAC +0.94. So
`SUBFAMILY_STRUCTURE` will fire on an old family, and will quietly stop firing
as copy number grows. The 6/8 detection and 0/28 false-positive rate in §18 were
measured at fixed n = 100 and similar ages, which hid both. **This flag must not
be trusted until it is normalised for age and n.**

**3. `cliff` is confounded by poly-A length.** POLYA **−0.99**, as strong as its
AGE −0.99. A longer poly-A tail extends past the consensus 3′ end into the flank
window, raising measured flank identity and shrinking the cliff. Since §12 showed
the consensus 3′ end is systematically ~7 bp short, this is not hypothetical.

**4. `arich_score` measures the TSD, not the A-rich zone.** Its only strong
response is TSD −0.85. It should be dropped in favour of the insertion-site
motif described in §14a.

### Which statistics are diagnostics and which are quality scores

`g_element` responds to 9 of 14 gradients, `score` and `rank1_frac` to 7, and
`g_homogeneity` to 6. These are general quality measures — useful for ranking
candidates, useless for saying *what is wrong*. The named sub-case flags in
`verdict.py` must be driven only by the specific statistics in the table above.

### Next

Single-feature gradients are done. The combination sweep is the obvious
follow-up: pairs of parameters on a factorial grid, to find where two
perturbations are mistaken for a third — the case that matters most is
age × subfamily divergence, since confound 2 above means an old single family
may be indistinguishable from two young subfamilies.

Data: `grad_response.json` (the full rho matrix), `grad_stats.json` (every
statistic on every set), `grad_truth.json` (labels), `grad_report.txt`.
