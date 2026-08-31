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
