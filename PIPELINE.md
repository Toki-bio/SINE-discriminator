# SINE discriminator — build 1 (2026-08-31)

Executes step 3 and step 1 of the spec's §6 next-steps list: implement `MEASURE()`
alone, no model, and run the §3 subsampling side quest against known families.

## Where things are

| | |
|---|---|
| Spec | `SINE_discriminator_spec.md` |
| Tal data (cloned) | `data/Tal/` — github.com/Toki-bio/Tal |
| KIT scripts | `kit/prep_sets.py`, `kit/run_mafft.sh` → uploaded to `/data/W/toki/` |
| KIT working dir | `/data/W/toki/SINE_disc/` (`sets/`, `aln/`, `manifest.json`, logs) |
| MEASURE | `measure.py` (local, numpy/scipy) |
| Inspection | `analyze.py` |

## Why the Tal repo alone was not enough

The published alignments (`*/alignments/*_rand100.aln.fa`) are **element-only**.
`WORKFLOW.md` step 1 extracts with a 50 bp flank, but the flank is trimmed before
these files are written, and 50 bp is in any case below what the spec calls for.
Tier 2 — the whole boundary/cliff criterion — is unmeasurable without flanks.

What the repo *does* give is the key: headers such as
`PVIJ01006492.1:46528-46784(+)|s1_30seqs|1788` are genomic coordinates, and the
same coordinates live in `results/assignment_full.tsv` in each SINEderella run
directory on KIT, alongside the genomes. So the corpus was rebuilt from source
with 300 bp flanks.

Note the skill's warning was load-bearing here: `regions.by_subfam.bed` is broken
(every line lists every subfamily). `assignment_full.tsv` filtered to
`Status == "assigned"` is the correct source.

## Corpus (580 sets, 4 species, 28 subfamilies)

Extraction: `assignment_full.tsv` → BED → `bedtools slop -b 300` →
`getfasta -s` (strand-corrected, so left flank is really 5′) → pool of 2000
copies per subfamily, stratified over bitscore deciles, contig-edge and
N-rich extractions dropped. Plus 6000 random genomic loci per species at the
median copy length, extracted identically.

| class | n sets | what it is |
|---|---|---|
| `POS` | 28 | one 100-copy sample per curated subfamily (saq 9, ccr 8, teu 6, dmo 5) |
| `SQ` | 400 | side quest: 5 families × n∈{25,50,100,200} × K=20 independent subsets |
| `NEGRAND` | 20 | 100 random genomic loci — the pure null |
| `NEGJITTER` | 28 | real copies, each edge independently displaced 20–100 bp |
| `NEGTRUNC5` | 28 | geometric 5′ truncation, foreign genomic sequence pasted on as the 5′ flank — LINE-fragment mimic |
| `NEGSPLICE` | 20 | chimeras: left half of family A + right half of family B |
| `MIXED10/30` | 56 | real family diluted with 10 % / 30 % random loci |

Alignment: `mafft --auto --adjustdirection`, 24-way parallel on KIT.

Sampling is drawn from the **raw assigned hit set**, not from the curated
published copy lists — the spec's §3 caution about subsets of a boundary-curated
family agreeing with the boundary rule by construction.

## Residual circularity, stated plainly

The anchors still come from `sear`'s own alignment to the consensus, so element
coordinates are seed-derived. The extraction widens by 300 bp on each side and
`MEASURE()` re-derives the edge from consensus information content, so the
*measured* boundary is independent — but the *centring* is not. Removing this
fully needs the spec's independent nhmmer re-search; `NEGJITTER` is the stand-in
that shows what a seed-independent, badly-centred set looks like.
