---
name: sinederella
description: What SINEderella is and does - the orchestrator, every step, the standalone steps, what a subfamily is defined as, and which of its tools NOT to reimplement. Read before any SINE search, consensus, boundary, clustering or alignment work.
---

# SINEderella

Toki's pipeline. **It already does search, consensus, boundary refinement,
subfamily assignment, subfamily clustering and alignment extraction.** Read this
before writing code for any of those — every time I wrote my own version it was
worse, and I have now done it with `SubFam`, `conse`, `sine_consensus.sh`,
`step7_boundary_refine.sh`, `step8a_extract_alignments.sh` and the peel loop.

Everything below was confirmed by reading the scripts line by line and, where a
number appears, by measuring it on a real run. Nothing here is inferred from
documentation alone; where his docs and his code disagree, both are recorded.

**Source of truth:** `/staging/tmp/SINEderella/` on DRAGEN (a git checkout, so
`README.md` and the 70 KB `MANUAL.md` travel with the code). Also
`/data/V/toki/bin/SINEderella/` on KIT. The de novo stage is a separate repo at
`/data/V/toki/bin/SINE-de-novo-genome-scan/` on KIT.

---

## What a subfamily *is* — the definition everything hangs on

MANUAL §6.1.6, his words:

> A *family* is a set of SINE copies that share the same sequence "parts" — the
> same structural composition, still meaningfully alignable/collinear across
> their full mutual length. A *subfamily* is a lineage within a family, defined
> by a **specific, shared, diagnostic pattern** of small indels/SNPs common to
> that lineage's copies (a synapomorphy, in effect) — **not** by the generic
> accumulation of private per-copy mutations from ordinary post-insertion decay,
> which is noise on top of the subfamily signal rather than the signal itself.

**Consequence: average pairwise identity is the wrong axis for subfamily
clustering.** It measures exactly the private decay the definition calls noise.
A clustering built on an identity matrix is answering a different question. His
`SINEClusterer` clusters on shared (position, character) features instead.

Distinguishing genuine diagnostic patterns from recurrent/convergent mutation
hotspots (homoplasy) is explicitly parked as future work, not solved.

---

## The one command

```bash
SINEderella <GENOME_FASTA> <CONSENSUS_FASTA>       # full run
SINEderella --add|-a <newSINE.fas> [--run <dir>]   # add consensuses, search only the new ones
SINEderella --exclude|-e <SINE_name> [--run <dir>] # drop a family, reassign the rest
SINEderella --resume|-r [--run <dir>]              # resume an interrupted run
```

Env: `THREADS` (default nproc), `CHUNK_BP` (30000), `FLANK` (50).

Creates `run_<YYYYMMDD_HHMMSS>/`, writes `manifest.txt`, and **copies every step
script into the run dir** so a run is reproducible from its own contents.

`mode_full` → sanitize genome (`genome.clean.fa`) and consensuses
(`consensuses.clean.fa`) → assembly QC (diagnostic, never blocks) →
**step1 → step2 → step3 → step4 → `results/` symlink farm → step6 report**.

**Steps 5, 7, 8a, 8b are standalone and are NOT run by the orchestrator.**
Check for their outputs before calling a run complete.

| step | script | what it does |
|---|---|---|
| 1 | `step1_search_extract.sh <cons> <genome> [sample=30000] [bin=50]` | `sear -k` per consensus, merge, extract, sample, SubFam, align chunks + consensus bank |
| 2 | `step2_asSINEment.sh <cons> <extracted> <threads> <outdir>` | 10-cycle unanimous voting assignment |
| 3 | `step3_postprocess.sh <rundir> <threads>` | per-copy similarity, LEAK, CONFLICT, soft assignment, summaries |
| 4 | `step4_plots.sh` | divergence histograms, nucleotide composition (optional) |
| 5 | `step5_align_subfamilies.sh <bank> <multi_dir>` | subfamily alignments — **has a defect, see below** |
| 6 | `step6_report.sh <rundir>` | `results/report.html`. Must run AFTER the symlink farm or the farm's `rm -rf results` destroys it |
| 7 | `step7_boundary_refine.sh <RUN_ROOT> [flank=50] [step=50] [max=1000]` | boundary refinement |
| 8a | `step8a_extract_alignments.sh <RUN_ROOT> <SPECIES_CODE>` | top100/rand100/subfam alignments |
| 8b | `step8b_publish_report.sh` | wires 8a's alignments into `report.html` as MSA-viewer links |

---

## `sear` — the search engine under step1

`sear [options] <query.fa> <bank.fa> [Slf=0.8] [Homology=65] [Flank=50] [KeepSplits=0]`

**Two engines by query length.** Queries ≥1000 bp → minimap2 (identity =
`matches / min(qaln, taln)`, so gaps are penalised). Shorter — i.e. all SINEs →
ssearch36. `--force-ssearch` overrides.

**Chunking, and why.** `seqkit sliding -g -s 2000 -W 2500` → 2500 bp windows at
2000 bp step, **500 bp overlap**, so a ~300 bp element is never cut across a
boundary without an intact copy in the neighbouring window. Then
`seqkit split2 -p 20`. Offsets are parsed back out of the `_sliding` headers and
added to hit coordinates at the end. `ssearch36` is Smith-Waterman with no index
and its statistics degrade on long sequences — that is the reason for chunking,
not sensitivity.

ssearch36 flags: `-g -3 -T <threads> -Q -n -z 11 -E 2 -w 95 -W 70 -m 8C`, ktup 3.

**`Slf` is a minimum hit LENGTH, not query coverage.** `flen = int(Slf × query_len)`;
a merged hit is kept if its span in the bank is ≥ `flen` **and** identity ≥
`Homology`. At 0.8 a hit must be ≥80 % of the consensus's length. Truncated
elements are dropped, not shortened.

**`-m/--minus` — a minus-bank is hit depletion, not a subsample.** After a round,
`bedtools complement` of the hits gives the **un-hit remainder**, `getfasta`
writes it, the next round searches that, repeating until a round returns zero.
Copies shadowed by stronger neighbours surface in later rounds. The
`chrom:start-end()` header form (empty parens = no strand column) comes from that
`getfasta`, which is why downstream scripts special-case it.

**`-b/--best N` builds the top-N alignment.** Top N unique by bitscore → ±50 bp
flanks → `getfasta -s` → **query prepended as row 1** → `mafft --localpair
--maxiterate 1000 --ep 0.123 --nuc --reorder` → `<tax>-<query>.best100.mafft.fa`.
Correctly oriented by construction. Do not rebuild this with blastn + bedtools +
`mafft --adjustdirection`; that is how orientation gets lost.

Other outputs: `<tax>-<query>.bed` (all merged hits), `<tax>-<query>.bnk` (hits +
flanks, stranded) — the `.bnk` that feeds SubFam. `@U@` is a header placeholder
for `_`, restored on output. `-s/--stop N` caps total hits, default 1000. `-k`
keeps splits so step1 reuses them across queries.

*Doc/code mismatch:* the usage text says the minimap2 cascade is
`asm5 -> asm10 -> asm20`; the loop runs only `asm10 asm20`. Harmless for SINEs.

---

## step1 — and the two things about it that matter

`sear -k "$query" "$GENOME" 0.8 65 50` per consensus — **no `-b`**, so step1 does
not produce top100 alignments.

**It merges hits across all consensuses without `-s`.** A locus hit by two
families on opposite strands therefore gets strand `+,-`.

**It ends by building the subfam alignment.** `input.clw` → degapped →
`input_reps.fasta` → concatenated with the **original consensus bank** → `mafft
--localpair --maxiterate 1000 --ep 0.123 --preservecase` →
`genome.clean_step1/subfam_input/input.clw.al`.

Timema `run_20260821_132226`: **605 rows = 597 chunk consensuses + the 8 family
consensuses**. That is the subfamily alignment with anchor rows, already built —
do not re-derive it from `*.bnk.cons`.

---

## `SubFam` — what it actually is

`SubFam <bank.fa> [chunk=50]`: `mafft --retree 0 --reorder` → `seqkit split2 -s 50`
→ align each chunk → **one consensus per chunk** (`cons -plurality 18`) → `.clw`
→ `mafft --localpair --maxiterate 1000 --ep 0.123` → `.msf`.

**A "subfamily alignment" is an alignment of chunk-consensus rows, not of copies.**
The consensuses fed to a run are *families*, not subfamilies.

**Chunking is by SIMILARITY, not input order.** `--reorder` emits in guide-tree
order, so each 50-block is a similarity neighbourhood. Measured on Timema's 550
labelled chunks:

| chunks apart | same subfamily |
|---|---|
| 1 | **0.916** |
| 2 | 0.89 |
| 10 | 0.846 |
| 50 | 0.74 |
| 200 | 0.530 |
| random pair | **0.313** |

Chunk index is real information; a clustering that ignores it discards a strong
prior.

**`-plurality 18` is hardcoded** — 18/50 = 36 %, correct only because chunks are
50. Change the chunk size and it silently becomes wrong.

**`.clw` is a plain `cat` of the per-chunk `.cons` files**, not an alignment
despite the extension. MANUAL §6.1.1 is explicit that it is not converged —
degap and realign before measuring anything off it:

```bash
seqkit seq -w 0 input.clw | seqkit seq -g > input.degapped.fasta
mafft --auto --quiet input.degapped.fasta > input.realigned.fasta
```

**Redundancy before SubFam is load-bearing, not an oversight** (§6.1.5). Guide-tree
order puts a family's near-identical copies in the same chunk, and `cons
-plurality 18` needs that volume to produce a well-evidenced chunk consensus.
An upstream CD-HIT-style dedup would starve it. Two distinct subfamilies that are
>80 % identical are precisely what step2's voting exists to separate *after*
extraction.

---

## `asSINEment` / step2 — how a copy gets its subfamily

`asSINEment <consensus.fa> <sines.fa> [threads] [outdir]`

1. Copies split into blocks of 20,000.
2. **10 cycles** of `ssearch36 ... -m 8`, consensuses as query; per cycle each
   copy's best hit by bitscore is one vote.
3. A copy is a candidate only on **10/10 unanimity** for one consensus.
4. Bar is **relative and per-subfamily**: sort that subfamily's unanimous
   bitscores descending, `N = min(10, count)`, threshold = **0.45 × the N-th
   best**. Below → `rejected_low_bitscore`.

Outputs `assigned.fasta` (`>seqID|subfamily|bitscore`), `unassigned.fasta`,
`assignment_stats.tsv`, `assignment_full.tsv`.

**The 10 cycles do real work — tested, not assumed.** Same command, same inputs,
three runs on 20,000 Timema copies: **43,396 / 43,484 / 43,657 rows, three
different md5s**, and **17 copies changed their winning consensus** between runs.
`ssearch36` is non-deterministic under `-T` with `-z 11`, so unanimity is a
genuine *stability* filter discarding copies whose assignment flips on
thread-scheduling noise.

---

## step3 — per-copy divergence, already computed

- **`sim_ratio`** = copy-vs-its-consensus ssearch36 bitscore ÷ that consensus's
  self-alignment bitscore. Timema: n=182,565, min 0.290, q1 0.599, **median
  0.748**, q3 0.897, max 1.028.
- **LEAK** = `runner_ratio >= 0.90` — the second-best consensus scores within 10 %
  of the best. Timema 9,340 copies. Leak % is a real-family signal: Timema
  families judged real ran 0.00–0.18 %, noisy ones 65–98 %.
- **CONFLICT** = the locus's region carries >1 subfamily label. Timema 82,365.
- Uses `-E 100` here, not `-E 2`.

---

## The `(+,-)` strand — measured, and what it really means

Timema `assigned.fasta`, 337,898 loci: `+` 128,115, `-` 128,908, **`+,-` 80,875
= 23.9 %**. step8a does `if (strand != "+" && strand != "-") strand = "+"`, so all
of those are extracted plus-strand.

Cost, measured as identity to the consensus row: `top100` 0 flipped rows in 7 of
8 subfamilies (t4: 13 low, 4 genuinely reverse); `rand100` 0 in 5 of 8 (t3: 25
low but only 4 reverse, t4: 4, t8: 1). **A real defect worth 0–4 rows per 100,
not a pipeline-wrecker.**

**`+,-` is the CONFLICT set.** 82,365 of 339,566 regions (24.3 %) carry more than
one subfamily label — the same population as the 23.9 % ambiguous strand. The
**800 highest-bitscore loci are all `+,-`**: the most complete copies are the ones
several related consensuses hit at once.

`regions.by_subfam.bed` is **not** broken in this run (contrary to the
`sine-orth-loc-sinederella` skill): 170,186 lines list `t1` alone, 33,393 `t8`,
21,199 `t7`. The multi-label sets are informative — the dominant one is
**`t3,t4,t5` at 69,585 loci**, which is his pipeline independently reporting the
t3/t4/t5 overlap that identity-based clustering also stumbles on.

---

## step7 — boundary refinement. Do not reimplement this.

- **Step A measures background empirically**: samples random genomic regions
  ≥200 bp, computes pairwise identity and `elevated_frac`, sets the pass
  threshold at **background + 5 %**. It does not assume 0.25.
- Then per **subfamily × population × side** it walks outward in `STEP_BP` windows
  to `MAX_EXT_BP` (1000), reverse-complementing minus-strand copies so every
  window is read in element orientation, stopping when `elevated_frac` reaches
  background.
- Output `step2/step2_output/boundary_refinement.tsv`: `subfamily | population |
  side | boundary_bp | status | final_identity_pct | final_elevated_frac_pct |
  background_identity_pct | background_elevated_frac_pct`. Status is
  `confirmed` / `undetermined` / `insufficient`.
- **Populations are tested separately** — top100 (bitscore-ranked) and
  rand100/general — because a boundary confirmed on a random sample is not
  automatically valid for the bitscore-ranked set.

Results: Timema 32/32 sides confirmed; a later run 24/24. Also run on hedgehog
and scorpion.

Using one symmetric flank, one population, and an assumed 0.25 background is
exactly how I produced a consensus with 21 bp of junk on the 5' end and its
poly-A tail cut off.

---

## step8a — the alignments to use

Base flanks **50 up, 70 down**, extended per subfamily from step7 but **only on
`status=="confirmed"` rows**; `undetermined` means the walk hit the 1000 bp cap
without confirming, and applying it once gave a 1070 bp flank that made
`--localpair --maxiterate 1000` crawl.

`top100` uses step7's **top100-population** rows; `rand100` and `subfam` use the
**general** rows.

- `top100`: 100 highest-bitscore members
- `rand100`: 100 shuffled members
- `subfam`: **only at ≥400 members** — seeded shuffle (42), up to 10,000, SubFam,
  degap `.clw`, align with the subfamily's consensus. Optional; never a required link.

`--reorder` can move the consensus off the top, so a python helper forces it back
to row 1 afterwards. `@U@` restored to `_` on write.

---

## `extract_alignments.sh` — and the flank trap

`extract_alignments.sh <consensus_bank.fa> <multi_run_dir>/ [<out_dir>]`

| file | content |
|---|---|
| `{sf}.core.fa` | 50 random copies, 30 bp left + 70 bp right flanks, MAFFT |
| `{sf}.best50.fa` | top 50 by bitscore, same flanks |
| `{sf}.subfam.fa` | SubFam chunk consensuses — **only at ≥400 copies** |
| `{sf}.evidence.txt` | zero-hit note |

**`postprocess_flanks()` destroys flank homology, deliberately.** The consensus
row's first and last non-gap column define the body. Outside it, per copy:
internal gaps deleted, bases lowercased, gaps padded on the **outer** side — left
flank right-justified against the body, right flank left-justified. Column count
preserved, column correspondence gone.

**This is the "no ragged tails" convention** — flanks butted against the element
edge for the eye. `ALIGNMENT_DEPLOYMENT.md` confirms the intent: the viewer "must
show sequences WITH lowercase flanking DNA before/after the SINE".

**Therefore any column-wise flank statistic — island z-scores, flank identity
decay, `decay_L`/`decay_R` — is meaningless on `.core.fa` / `.best50.fa`.** It is
valid on step8a's `*_top100.aln.fa` / `*_rand100.aln.fa`, which do not get this
treatment.

---

## step5 — two defects, use step8a instead

`ALIGNMENT_DEPLOYMENT.md` says to prefer `extract_alignments.sh` because step5
produces no flanks. The code shows two more concrete problems:

- Large subfamilies (≥400) get `mafft --add <bank> --keeplength` — **`--keeplength`
  forces the added consensus into the existing column frame, deleting consensus
  residues that would need new columns.**
- Small subfamilies get `--auto`, not `--localpair --maxiterate 1000 --ep 0.123`.

---

## Consensus building

**`conse <msa.fasta> [pct]`** — EMBOSS `cons` at a plurality of **pct % of the
number of sequences, default 35 %**, then **N → `-`**, one-line FASTA. Because
uncalled columns become gaps, the boundary emerges by itself. **A consensus
threshold is a fraction of the sequences, never a count** — reusing SubFam's
`-plurality 18` on a 20-copy family demands near-unanimity and yields chimeras.

**`sine_consensus.sh input.fasta [subsample=100] [max_iters=50] [thresh=0.01]`** —
**positional only, no flags.** Bootstrap: subsample → `mafft --auto` →
mini-consensus → append to master → realign → converge on Hamming/len < 0.01,
minimum 5 iterations. Column rule: **>50 % gaps → `-`**, otherwise plurality among
non-gap bases (gaps excluded from that denominator), ties → `N`. Says explicitly:
*"No automatic trimming. Trim manually after visual inspection."*

**`sine_consensus_smart.sh [options] input.fasta`** — adds variance detection,
early stopping and a `QUALITY: LOW_CONFIDENCE` flag. `-n` subsample, `-m` max
iters, `-s` stability, **`-t`** call threshold (50), **`-c`** min coverage (30),
`-k` keep. Slow convergence (>30 iters) or `NOT_CONVERGED` means mixed/FP data;
`analyze_convergence.sh` batch-summarises this.

**`sine_pairwise_consensus.sh [options] input.fasta`** — bootstrap a reference,
align each copy to it with `--op 5`, project onto reference coordinates, count
with gaps **in** the denominator, track insertion hotspots. The reference is a
coordinate frame only; its bases do not vote. `-r` reference, `-n` max copies,
`-t` call threshold, `-c` min coverage, `-j` jobs, `-N/-M/-S` bootstrap params.

*Doc/code mismatch:* its help says `-t` defaults to 50 and `-c` to 30; the code
sets `CONS_THRESH=30`, `MIN_COVERAGE=10`. Pass them explicitly.

**So "trim the flanks" is `-t 70 -c 30`.** There is no `-c 70 -g 30` — I invented
that flag pair; it does not exist in any of these scripts.

His rule on which to use: **approximate consensus → the bootstrap/pairwise
method; robust and data-supported → `conse` over the top 100 hits by bitscore**,
which is less exposed to winner's curse. Run that final pass **once** — iterating
creeps outward (203 → 252 → 269 bp on hyd_SINE_0).

**Do not bolt on an automated consensus-refinement loop** (§6.1.5, explicit):
RepeatModeler2 `Refiner`-style `refineUntil` or TEstrainer BEAT loops are
unreliable alone — they still need human judgement to catch a consensus that
"converged" confidently onto a chimeric or boundary-contaminated alignment.

---

## The manual subfamily step (§6.1) and `SINEClusterer`

This is the step being automated, and it is **his manual workflow**:

1. Degap + realign `input.clw` (§6.1.1 — never cluster off the raw `.clw`).
2. Open `input.realigned.fasta` in MSA-viewer, visually identify groups sharing
   diagnostic columns/motifs and a consistent indel pattern (§6.1.2).
3. Extract each candidate group to its own FASTA.
4. Build a consensus per group with `sine_consensus.sh` / `_smart` (§6.1.3).
5. Combine them into a new bank and re-run SINEderella, or run
   `step2_asSINEment.sh` directly against the existing `extracted.fasta` (§6.1.4).

**`step1b_cluster_subfamilies_assist.sh input.clw out_dir [--recurse]`** runs
`cluster_assist.js` → `subfam_cluster_lib.js`'s `SINEClusterer`, vendored from
MSA-viewer. It is **an assist, not a replacement**, and its limits are measured:

> Validated 2026-08-21 against a real 9-subfamily ground truth (Tal/saq,
> s1–s9 over 600 chunk consensuses): recovered a true 2-subfamily split
> correctly, but on the 9-subfamily case collapsed everything into **2 giant
> blobs (293 + 307) plus one 6-seq cluster**, and `--recurse` recovered nothing
> further.

**How `SINEClusterer` works** — this is the peel loop, already written:

- A candidate group is the **exact set of sequences sharing one (position,
  character)**. Gaps are never diagnostic. Columns where one base holds >80 % of
  the *whole* alignment are skipped as non-diagnostic (checked globally, not
  against the remaining pool — checking the pool made a clean remainder look
  non-diagnostic once earlier clusters were peeled).
- Candidate sets are fuzzy-merged at Jaccard ≥0.90 with size difference ≤5.
- Group size is capped at **50 % of available** to prevent degenerate
  "everything" clusters, with a relaxed-upper-bound retry when nothing is found.
- Features are scored: `outside == 0` → **3** if every member matches, **2** at
  ≥80 %, **1.5** otherwise; else `qual = inside% − outside%` above a size-dependent
  threshold → **1**. A group needs `minPerfect` good features.
- Outliers are pruned: members matching fewer than ~30 % of the group's features
  are dropped, then the group is re-validated.
- **Peel loop**: take the best group, remove it, repeat, up to 20 iterations,
  with thresholds **relaxing** each round (`minPerfect` 4 → 1, quality thresholds
  down by up to 50 points) and an ultra-relaxed rescue at ≤10 remaining.

Ground truth available for benchmarking a replacement: Timema t1–t8 (chunk
purity **0.953** via copy-majority — reliable) and Tal/saq s1–s9. **The saq
chunk→group labels are not on disk** in any form I have found: copy-majority
gives purity 0.757 with sizes that contradict the group names (s6 → 90 chunks
against the 7 its name states), and best-hit to the 9 consensuses lands s2 and
s6 exactly but leaves s4 and s9 empty. Ask before treating either as truth.

---

## SINE-de-novo-genome-scan — the stage before SINEderella

Finds SINEs with **no library**, on the premise that most SINEs share at least
faint similarity with some already-known SINE.

**1. Shred a SINE database** (`01_sanitize_shred.sh`, `02_filter_lc_cluster_nr85.sh`).
Input deduplicated at 95 % with vsearch. Each SINE of length L is cut
region-aware, because head and tail are the diagnostic parts:

| region | span | step |
|---|---|---|
| 5' | bases 1–150 | 10 bp |
| middle | 151 to L−100, only if L > 250 | 25 bp |
| 3' | last 99 bp, only if L > 150 | 25 bp |

Headers `>SINE_ID|REGION:START-END`; IDs sanitised (`|`, `:` → `_`) because
EMBOSS and downstream parsers need it. Then low-complexity filtered, clustered
to `nr85`.

**2. Scan** (`sine_scan.sh`): per-fragment Smith-Waterman with `ssearch36`, filter
on identity and **query** coverage, merge hits within 500 bp into candidate loci,
extract with flanks. Outputs `query_summary.tsv`, `all_hits.filtered.m8`,
`candidate_loci.bed`, `candidates.fa` → feeds SubFam.

Two versions exist: the newer takes `-q -g -o` and parallelises **per query**;
the scorpio copy takes positional args and chunks the **database** into 100 Mb
blocks under GNU parallel.

**Consequence for the discriminator:** these loci are **anchors expanded to
candidates**, not verified full elements. A de novo family scoring
FRAGMENT_OF_LONGER is describing its input honestly, not disagreeing with the scan.

---

## Two competing models of divergence (RESEARCH_DIRECTIONS.md)

Recorded so it is not lost, and because it frames what the discriminator is
allowed to assume.

Everything in the pipeline today assumes a **gradual, substitution-driven** model:
relatedness is a matter of degree, measured as identity. The competing model is
**modular / "repeat pangenome" (panconsensus)**: some repeats are composites
assembled from a library of reusable, independently-mobile modules, so
relatedness is a **presence/absence and arrangement matrix over a module
library**, not one number.

Motivating evidence already in hand: the `eri` e2-3 downstream flank is
**bimodal**, not gradually diverging — roughly half of sampled copies still share
an exact motif 100–120 bp past the called SINE end while the rest look like
background at the same offset. The gradual framing called that boundary
under-calling and needed a mean → fraction-of-pairs fix to detect it at all.

The two predict different signatures: gradual under-calling correlates boundary
distance with subfamily age; a present/absent module gives a clean two-cluster
split with no gradient. Proposed as an **additive standalone step** (segment →
cluster blocks not copies → build presence matrix → compare model fit), to be
tried first on `eri` e2-3. Status: idea only as of 2026-08-21, no code.

---

## Reading a run

`results/`: `assignment_full.tsv` (Sequence, Subfamily, Bitscore, Votes, Status,
Threshold, RunnerInfo), `summary.by_subfam.tsv`, `subfamilies/<sf>.fasta`
(element-only, headers `>chrom:start-end(strand)|subfam|bitscore`, `_` escaped as
`@U@`), `assigned.fasta`, `unassigned.fasta`, `consensuses.fa`.

- **No run has `manifest.txt` at the root in the shape a `run_2026*/manifest.txt`
  glob expects** — assignment output is under `step2/step2_output/`.
- Cluster IDs `C<N>R` are **not stable** across runs — never key a lookup on them.

## Environment

Needs `ssearch36`, `mafft`, `bedtools`, `samtools`, `seqkit`, `cons` (EMBOSS),
Python 3 + matplotlib. Not all on the default PATH on DRAGEN:

```bash
export PATH=/staging/conda/envs/bioinfo/bin:/staging/miniconda3/bin:/usr/local/bin:$PATH
```

`node` is **not** installed on either DRAGEN or KIT, so `cluster_assist.js` cannot
be run there as-is.

`SubFam` has a `/bin/sh` shebang but bash syntax — on KIT (dash) run it as
`bash /data/V/toki/bin/SubFam`. Never write run data to `/home` on DRAGEN.

Run locations are in `DATA_LOCATIONS.md`. See also the
`sine-orth-loc-sinederella` skill for `SINE_orth_loc.bash`, PM/MP/SINE orthology,
and the three non-interchangeable subfamily labelling systems.

---

## Traps, in one place

1. **Do not cluster subfamilies on pairwise identity** — §6.1.6 says that is the
   noise, not the signal. Cluster on shared diagnostic columns.
2. **Do not measure column-wise flank statistics on `.core.fa` / `.best50.fa`** —
   `postprocess_flanks` has already destroyed column correspondence there.
3. **Do not orient subfam chunk consensuses.** They inherit the input
   consensuses' orientation. Flipping them fused Timema t3 with t4 and split t1
   in two. Orientation *is* needed for AnnoSINE seed candidates, where about half
   are reverse-complemented — different situation.
4. **Do not use single linkage.** It chained Timema t3/t4/t5 into one 143-chunk
   cluster whose median within-cluster identity was 0.983 — the highest of any
   cluster — while being 40 % pure. No within-cluster statistic can see a join
   made one link at a time.
5. **Do not treat `input.clw` as an alignment.** Degap and realign first.
6. **Do not dedup before SubFam.** The redundancy is what `-plurality 18` needs.
7. **Do not iterate the final `conse` pass.** Once. It creeps outward otherwise.
8. **Do not apply an `undetermined` step7 boundary.** Fall back to the base flank.
9. **Check step7 population** before reusing a boundary: top100 ≠ general.
