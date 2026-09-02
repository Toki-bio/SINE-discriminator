---
name: sinederella
description: What SINEderella is and does - the orchestrator, every step, the standalone steps, and which of its tools NOT to reimplement. Read before any SINE search, consensus, boundary or alignment work.
---

# SINEderella

Toki's pipeline. **It already does search, consensus, boundary refinement,
subfamily assignment and alignment extraction.** Read this before writing any
code for those things — every time I have written my own version of one of its
tools it has been worse, and I have now done it with `SubFam`, `conse`,
`sine_consensus.sh`, `step7_boundary_refine.sh` and `step8a_extract_alignments.sh`.

Source of truth: `/staging/tmp/SINEderella/` on DRAGEN (a git checkout, so
`README.md` and the 70 KB `MANUAL.md` travel with the code). Also
`/data/V/toki/bin/SINEderella/` on KIT.

## The one command

```bash
SINEderella <GENOME_FASTA> <CONSENSUS_FASTA>       # full run
SINEderella --add|-a <newSINE.fas> [--run <dir>]   # add consensuses, search only the new ones
SINEderella --exclude|-e <SINE_name> [--run <dir>] # drop a family, reassign the rest
SINEderella --resume|-r [--run <dir>]              # resume an interrupted run
```

Env: `THREADS` (default nproc), `CHUNK_BP` (30000), `FLANK` (50).

Creates `run_<YYYYMMDD_HHMMSS>/` in the working directory, writes `manifest.txt`
(run, pwd, script dir, resolved genome and consensus paths, threads, chunk,
flank, mode, date), and **copies every step script into the run dir as a frozen
copy** so a run is reproducible from its own contents.

## What the orchestrator runs, and what it does not

`mode_full` → sanitize genome (`genome.clean.fa`) and consensuses
(`consensuses.clean.fa`) → assembly QC (diagnostic, never blocks; verdict
WARN/CAUTION/OK in `genome.clean_assembly_qc.tsv`) → **step1 → step2 → step3 →
step4 → `results/` symlink farm → step6 report**.

**Steps 5, 7, 8a, 8b are standalone and are NOT run by the orchestrator.**
Check for their outputs before calling a run complete.

| step | script | what it does |
|---|---|---|
| 1 | `step1_search_extract.sh <cons> <genome> [chunk=30000] [flank=50]` | `sear -k <query> <genome> 0.8 65 50` per consensus (hit must cover 80 % of query at 65 % identity), merge overlaps, extract with flank, sample to `sample_size` (30000), **SubFam** on the sample, then `mafft --localpair --maxiterate 1000` |
| 2 | `step2_asSINEment.sh <cons> <extracted> <threads> <outdir>` | assign each copy to a consensus: 10 independent ssearch36 voting cycles, unanimity + bitscore threshold |
| 3 | `step3_postprocess.sh <rundir> <threads>` | per-copy similarity, LEAK detection, CONFLICT flags, soft assignment, summary tables |
| 4 | `step4_plots.sh` | divergence histograms, nucleotide composition (optional, warns) |
| 5 | `step5_align_subfamilies.sh <bank> <multi_dir>` | subfamily alignments **without flanks** — see the warning below |
| 6 | `step6_report.sh <rundir>` | `results/report.html`. Must run AFTER the results symlink farm, or the farm's `rm -rf results` destroys it — a bug hit on three separate runs before being fixed |
| 7 | `step7_boundary_refine.sh <RUN_ROOT> [flank=50] [step=50] [max=1000]` | **boundary refinement — see below** |
| 8a | `step8a_extract_alignments.sh` | core/best50/subfam alignments, flanks sized per subfamily from step7's table |
| 8b | `step8b_publish_report.sh` | wires 8a's alignments into `report.html` as MSA-viewer links |

## step7 — boundary refinement. Do not reimplement this.

This is the tool for "where does the element actually end", and it is better
than anything I built.

- **Step A measures background empirically**: samples random genomic regions
  ≥200 bp, computes their pairwise identity and `elevated_frac`, and sets the
  pass threshold at **background + 5 %**. It does not assume 0.25.
- Then per **subfamily × population × side** it walks outward in `STEP_BP`
  windows up to `MAX_EXT_BP` (1000), reverse-complementing minus-strand copies
  so every window is read in element orientation, and stops when
  `elevated_frac` falls to background.
- Output `step2/step2_output/boundary_refinement.tsv`: `subfamily | population |
  side | boundary_bp | status | final_identity_pct | final_elevated_frac_pct |
  background_identity_pct | background_elevated_frac_pct`.
- Status is `confirmed` / `undetermined` / `insufficient`.

**Populations are tested separately** — top100 (bitscore-ranked) and
rand100/general — because "a boundary confirmed on a random sample is not
automatically valid for the bitscore-ranked set".

Results seen: Timema 32/32 sides confirmed, 0 undetermined; a later run 24/24.
Also run on hedgehog and scorpion.

**Why this matters:** using one symmetric flank for both ends, one population,
and an assumed 0.25 background is exactly how I produced a consensus with 21 bp
of junk on the 5′ end and its poly-A tail cut off.

## Alignments: `extract_alignments.sh`, not step5

`ALIGNMENT_DEPLOYMENT.md` is explicit: **"CRITICAL: Use extract_alignments.sh,
NOT step5_align_subfamilies.sh"** — step5's output has no flanking sequence, and
flanks are mandatory for boundary inspection.

```bash
extract_alignments.sh <consensus_bank.fa> <multi_run_dir>/ [<out_dir>]
```

| file | content |
|---|---|
| `{sf}.core.fa` | 50 random copies, 30 bp left + 70 bp right flanks, MAFFT |
| `{sf}.best50.fa` | top 50 by bitscore, same flanks |
| `{sf}.subfam.fa` | SubFam chunk consensuses — **only at ≥400 copies**, so optional, never a required link |
| `{sf}.evidence.txt` | zero-hit note |

Layout discovered as `<multi_dir>/<Species>/run_*/`, then that run's
`step2/step2_output/subfamilies/` or `results/subfamilies/`. For a single run,
symlink it into that shape.

## Consensus building — `conse` and SINE_consensus

**`conse <msa.fasta> [pct]`** (`/usr/bin/conse` on DRAGEN,
`/data/V/toki/bin/` on KIT): EMBOSS `cons` at a plurality of **pct % of the
number of sequences, default 35 %**, then **N → `-`**, one-line FASTA out.
Because uncalled columns become gaps, the boundary emerges by itself.
**A consensus threshold is a fraction of the sequences, never a count** — reusing
SubFam's `-plurality 18` (a third of a 50-copy chunk) on a 20-copy family
demands near-unanimity and produces chimeras.

**`sine_consensus.sh`** (repo `Toki-bio/SINE_consensus`): bootstrap — subsample,
align, mini-consensus, accumulate, re-align, converge on Hamming < 0.01. Gaps
**excluded** from the denominator, which is what makes it full-length. `-c 70 -g
30` is the documented way to trim flanks. Falls back to plurality ≥35 % with
gaps in denominator if it does not converge.

**`sine_pairwise_consensus.sh`** for high-copy families: bootstrap a reference,
align each copy to it with `--op 5`, project onto reference coordinates, count
with gaps **in** the denominator, track insertion hotspots. The reference is a
coordinate frame only; its bases do not vote.

His rule on which: **approximate consensus → the bootstrap/pairwise method;
robust and data-supported → `conse` over the top 100 hits by bitscore**, which is
less exposed to winner's curse. Run that final pass **once** — iterating creeps
outward (203 → 252 → 269 bp on hyd_SINE_0).

## SubFam — what it actually is

`SubFam <bank.fa> [chunk=50]`: rough mafft reorder → `seqkit split2 -s 50` →
align each chunk → **one consensus per chunk** (`cons -plurality 18`) → align
those to each other (`.clw`, then `.msf`).

**A "subfamily alignment" is an alignment of chunk-consensus rows, not of
copies.** The consensuses fed to a run are *families*, not subfamilies.

**`input.clw` is NOT a converged alignment** (MANUAL §6.1.1) — degap and realign
before using it:
`seqkit seq -w 0 input.clw | seqkit seq -g > d.fa && mafft --auto d.fa > r.fa`

**Do not orient subfam chunk consensuses.** His words: *"subfam after sinederella
usually doesnt need flip because consensuses used are expected to be properly
oriantated"*. Orienting them fused two of his Timema subfamilies and split
another in two. Orientation is for AnnoSINE seed candidates, where about half are
reverse-complemented — not for post-SINEderella chunks.

## Reading a run

`results/`: `assignment_full.tsv` (Sequence, Subfamily, Bitscore, Votes, Status,
Threshold, RunnerInfo), `summary.by_subfam.tsv`, `subfamilies/<sf>.fasta`
(element-only, headers `>chrom:start-end(strand)|subfam|bitscore`, `_` escaped
as `@U@`), `assigned.fasta`, `unassigned.fasta`, `consensuses.fa`.

- **`regions.by_subfam.bed` is broken** — every line lists all subfamily names.
  Use `assignment_full.tsv` or `all_hits.labeled.bed`.
- **No run has `manifest.txt` at the root in the shape a `run_2026*/manifest.txt`
  glob expects** — assignment output is under `step2/step2_output/`.
- **leak %** is a real-family signal: Timema families judged real ran 0.00–0.18 %,
  noisy ones 65–98 %.
- Cluster IDs `C<N>R` are **not stable** across runs — never key a lookup on them.

## Environment

Needs `ssearch36`, `mafft`, `bedtools`, `samtools`, `seqkit`, `cons` (EMBOSS),
Python 3 + matplotlib. Not all on the default PATH on DRAGEN:

```bash
export PATH=/staging/conda/envs/bioinfo/bin:/staging/miniconda3/bin:/usr/local/bin:$PATH
```

`SubFam` has a `/bin/sh` shebang but bash syntax — on KIT (dash) run it as
`bash /data/V/toki/bin/SubFam`. Never write run data to `/home` on DRAGEN.

Run locations are in `DATA_LOCATIONS.md` in the SINE_discriminator repo. See also
the `sine-orth-loc-sinederella` skill for `SINE_orth_loc.bash`, PM/MP/SINE
orthology, and the three non-interchangeable subfamily labelling systems.

## SINE-de-novo-genome-scan — the stage before SINEderella

`/data/V/toki/bin/SINE-de-novo-genome-scan/` on KIT (a copy of `sine_scan.sh`
also sits in `/data/W/toki/scorpio/denovo_scan/`). Finds SINEs in a genome with
**no library**, on the premise that most SINEs share at least faint similarity
with some already-known SINE.

Two stages:

**1. Shred a SINE database into fragments** (`01_sanitize_shred.sh`, then
`02_filter_lc_cluster_nr85.sh`). Input deduplicated at 95 % with vsearch. Each
SINE of length L is cut region-aware, because the head and tail are the
diagnostic parts:

| region | span | step |
|---|---|---|
| 5′ | bases 1-150 | 10 bp |
| middle | 151 to L-100, only if L > 250 | 25 bp |
| 3′ | last 99 bp, only if L > 150 | 25 bp |

Headers `>SINE_ID|REGION:START-END`. IDs sanitised (`|` and `:` → `_`, uniquified)
because EMBOSS and downstream parsers need it. Then low-complexity filtered and
clustered to `nr85`.

**2. Scan the genome** (`sine_scan.sh`): per-fragment Smith-Waterman with
`ssearch36`, filter on identity and **query** coverage, merge hits within 500 bp
into candidate loci, extract with flanks. Outputs `query_summary.tsv`,
`all_hits.filtered.m8`, `candidate_loci.bed`, `candidates.fa` → feeds SubFam.

Two versions exist: the newer takes `-q -g -o` and parallelises **per query**;
the scorpio copy takes positional args and chunks the **database** into 100 Mb
blocks under GNU parallel. Same constraint either way - `ssearch36` has no index,
so the work has to be split somewhere.

**"minus-bank"** is simply a search DB whose headers are coordinate-wrapped,
`Scaffold_8:0-77793632()` - what `bedtools getfasta` writes when the BED has no
strand column. It comes from the optional genome subsampling step (by fraction
or fixed bp, random scaffolds). The script detects it from the first header and
then *requires* `TARGET_GENOME` to be set, because hits must be mapped back to
original scaffold coordinates:
`shift = interval_start0 + chunk_off1 - 1`.

Consequence for the discriminator: these loci are **anchors expanded to
candidates**, not verified full elements. That is why a de novo family scores
FRAGMENT_OF_LONGER so often - it is describing the input honestly.

## Corrections to earlier notes in this file

**SubFam chunks by SIMILARITY, not input order.** `mafft --retree 0 --reorder`
emits in guide-tree order, so each 50-sequence block is a similarity
neighbourhood. Measured on Timema: adjacent chunks share a subfamily **92 %** of
the time, decaying with distance (0.89 at 2 apart, 0.85 at 10, 0.74 at 50, 0.53
at 200) against **0.31** for a random pair. So chunk index is real information
and a clustering that ignores it is discarding a strong prior.

**`-plurality 18` in SubFam is hardcoded, not scaled.** It is 18/50 = 36 %, which
matches `conse`'s 35 % only because chunks are 50. Change the chunk size and that
constant silently becomes wrong.

**SubFam's `.clw` is a plain `cat` of the per-chunk `.cons` files** - not an
alignment despite the extension. The aligned product is `.msf`, built with
`mafft --localpair --maxiterate 1000 --ep 0.123` (note the non-default gap
extension penalty).
