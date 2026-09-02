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

## `sear` — the search engine under step1. Read it before touching hits.

`/staging/tmp/SINEderella/sear`, 653 lines of bash.
`sear [options] <query.fa> <bank.fa> [Slf=0.8] [Homology=65] [Flank=50] [KeepSplits=0]`

**Two engines, chosen by query length.** Queries ≥1000 bp go to minimap2
(cascade over presets, identity = `matches / min(qaln,taln)` so gaps are
penalised); shorter queries — i.e. all SINEs — go to the ssearch36 path.
`--force-ssearch` overrides.

**How the genome is chunked, and why.**
`seqkit sliding -g -s 2000 -W 2500` → 2500 bp windows at 2000 bp step, so
**500 bp of overlap**: a ~300 bp element can never be cut across a boundary
without an intact copy in the neighbouring window. Then `seqkit split2 -p 20`.
Window offsets are parsed back out of the `_sliding` headers and added to the
hit coordinates at the end. That is the answer to "why not an indexed search" —
`ssearch36` is Smith-Waterman with no index, and its statistics degrade on huge
sequences, so the bank is cut into pieces small enough to score properly and the
coordinates are re-assembled afterwards.

ssearch36 flags: `-g -3 -T <threads> -Q -n -z 11 -E 2 -w 95 -W 70 -m 8C`, ktup 3.

**`Slf` is a minimum hit LENGTH, not query coverage.** `flen = int(Slf × query_len)`
and a merged hit is kept if its span in the bank is ≥ `flen` **and** identity ≥
`Homology`. At the default 0.8 a hit must be at least 80 % of the consensus's
length. Matters for truncated elements: they are dropped, not shortened.

**`-m/--minus` — what a minus-bank actually is.** Not a subsample. After a round
of hits, `bedtools complement` of those hits against the bank gives the
**un-hit remainder**, `getfasta` writes it out, and the next round searches that.
Repeats until a round returns zero. It is **hit depletion**, so copies shadowed
by stronger neighbours surface in later rounds. The `chrom:start-end()` header
form comes from that `getfasta` (empty parens = no strand column), which is why
downstream scripts special-case it. *(My earlier note guessed genome subsampling.
Wrong.)*

**`-b/--best N` already builds the top-100 alignment.** Top N unique hits by
bitscore → ±50 bp flanks → `getfasta -s` → **the query prepended as row 1** →
`mafft --localpair --maxiterate 1000 --ep 0.123 --nuc --reorder`. Output
`<tax>-<query>.best100.mafft.fa`.

**So `sear -b100` IS the top100 consensus-anchored alignment.** Correctly
oriented by construction, because every copy is extracted stranded and the
consensus is row 1. Do not rebuild this with blastn + bedtools + `mafft
--adjustdirection` — that is exactly how the orientation was lost.

Other outputs: `<tax>-<query>.bed` (all merged hits), `<tax>-<query>.bnk`
(all hits + flanks, stranded) — the `.bnk` that feeds SubFam.

`@U@` is a header placeholder for `_`, restored on output; `-s/--stop N` caps
total hits at 1000 by default; `-k` keeps splits so step1 reuses them across
queries.

**Doc/code discrepancy:** the usage text says the minimap2 cascade is
`asm5 → asm10 → asm20`; the loop runs only `asm10 asm20`. Harmless for SINEs
(they never reach minimap2 mode) but the help is wrong.

## `asSINEment` / step2 — how a copy gets its subfamily

`asSINEment <consensus.fa> <sines.fa> [threads] [outdir]`, 440 lines.

1. Copies split into blocks of 20 000.
2. **10 cycles** of `ssearch36 ... -m 8`, consensuses as query, copies as bank;
   per cycle, each copy's best hit by bitscore is one vote.
3. A copy is a candidate only on **10/10 unanimity** for one consensus.
4. Threshold is **relative and per-subfamily**: take that subfamily's unanimous
   copies, sort bitscores descending, `N = min(10, count)`, and the bar is
   **0.45 × the N-th best bitscore**. Below it → `rejected_low_bitscore`.

Outputs `assigned.fasta` (`>seqID|subfamily|bitscore`), `unassigned.fasta`,
`assignment_stats.tsv`, `assignment_full.tsv`.

**Open question, not yet tested:** the 10 cycles run an identical deterministic
command, so unanimity can only fail if `ssearch36` is non-deterministic under
`-T <threads>` with `-z 11` (thread-order affecting the fitted statistics and
thus the `-E 2` cutoff for borderline hits). If that is the mechanism, the vote
is a *stability* filter and does real work; if ssearch36 is fully deterministic
here, the 10 cycles cost 10× runtime for nothing. Testable by running one cycle
twice and diffing. Ask before assuming either way.

## Read the rest (2026-09-02). Corrections and measurements.

### step1 already builds the subfam alignment you keep asking for

`step1_search_extract.sh` calls `sear -k <query> <genome> 0.8 65 50` — **no `-b`**,
so step1 does not make top100 alignments. What it *does* make, at the end, is:

`genome.clean_step1/subfam_input/input.clw.al` = SubFam's chunk consensuses
(degapped from `input.clw`) **concatenated with the family consensus bank** and
aligned `--localpair --maxiterate 1000 --ep 0.123 --preservecase`.

Timema `run_20260821_132226`: **605 rows = 597 chunk consensuses + the 8 family
consensuses** (`>t2` etc. are in the file). That is the subfam alignment with
anchor rows, already built. I was re-deriving it from `*.bnk.cons`.

### The `(+,-)` strand, measured

step1 merges hits across all consensuses **without `-s`**, so a locus hit by two
families on opposite strands gets strand `+,-`. step8a then does
`if (strand != "+" && strand != "-") strand = "+"`.

Timema, `assigned.fasta` — 337,898 loci: `+` 128,115, `-` 128,908, **`+,-` 80,875
= 23.9 %**, all extracted plus-strand.

Actual cost, measured by identity to the consensus row: `top100` 0 flipped rows
in 7 of 8 subfamilies (t4: 13 low, 4 genuinely reverse); `rand100` 0 in 5 of 8
(t3: 25 low but only 4 reverse, t4: 4, t8: 1). **So it is a real defect worth
0–4 rows per 100, not a pipeline-wrecker.** It was not the cause of the hydra
orientation problem — that was my own code.

**And `+,-` is not noise: it is the same set as CONFLICT.** 82,365 of 339,566
regions (24.3 %) carry more than one subfamily label — near-identical to the
23.9 % ambiguous strand. The 800 highest-bitscore loci are **all** `+,-`: the
strongest, most complete copies are the ones several related consensuses hit.

### `regions.by_subfam.bed` is NOT broken in this run

The `sine-orth-loc-sinederella` skill says every line lists all subfamilies.
Not here: 170,186 lines list `t1` alone, 33,393 `t8`, 21,199 `t7`. The multi-label
sets are informative — the dominant one is **`t3,t4,t5` at 69,585 loci**.

That is his pipeline independently reporting the same thing my peel loop found:
t3/t4/t5 overlap. t3 measured 0.685 internal identity and was closer to t5
(0.727) than to itself. Two independent methods agree, so it is a property of
the families, not a clustering failure.

### step3 already computes per-copy divergence — do not reinvent it

- `sim_ratio` = copy-vs-its-consensus ssearch36 bitscore ÷ that consensus's
  self-alignment bitscore. Timema: n=182,565, min 0.290, q1 0.599, **median
  0.748**, q3 0.897, max 1.028.
- **LEAK** = `runner_ratio >= 0.90` — the second-best consensus scores within
  10 % of the best. Timema 9,340 copies.
- **CONFLICT** = the locus's region carries >1 subfamily label. Timema 82,365.
- Uses `-E 100` here, not `-E 2`.

### The 10 voting cycles do real work — tested, not assumed

`ssearch36 -g -3 -T 32 -Q -n -z 11 -E 2 -w 95 -W 70 -m 8`, same inputs, three runs:
**43,396 / 43,484 / 43,657 rows, three different md5s**, and **17 of 20,000 copies
changed their winning consensus** between run 1 and run 2. So it is
non-deterministic under threading, and 10/10 unanimity is a genuine *stability*
filter that discards copies whose assignment flips on thread-scheduling noise.
Not wasted runtime.

### `extract_alignments.sh` post-processes flanks — its columns are NOT homologous

`postprocess_flanks()`: the consensus row's first and last non-gap column define
the body. Outside it, **per copy, internal gaps are deleted, bases lowercased,
and gaps padded on the outer side** — left flank right-justified against the
body, right flank left-justified. Column count preserved, homology destroyed.

**This is the "no ragged tails" convention.** Flanks are butted against the
element edge for the eye, not aligned.

Consequence: any column-wise flank statistic — island z-scores, flank identity
decay, my `decay_L`/`decay_R` — is **meaningless on `.core.fa` / `.best50.fa`**.
It is valid on step8a's `*_top100.aln.fa` / `*_rand100.aln.fa`, which do **not**
get this treatment. No such `.core.fa` files exist on DRAGEN yet, so nothing I
measured is affected — but this is a trap the moment extract_alignments is run.

### step5 truncates the consensus — a second reason to prefer step8a

Large subfamilies (>=400): `mafft --add "$CONS_BANK" --keeplength`. `--keeplength`
forces the added consensus into the existing column frame, **deleting consensus
residues that would need new columns**. Small subfamilies get `--auto`, not
`--localpair --maxiterate 1000 --ep 0.123`. step8a cats and fully re-aligns instead.

### Consensus-builder flags — my earlier note had them wrong

`sine_consensus.sh input.fasta [subsample=100] [max_iters=50] [thresh=0.01]`.
**Positional only, no flags.** Column rule: **>50 % gaps -> `-`**, otherwise
plurality of non-gap bases, ties -> `N`. Gaps *are* in the denominator for the
gap test. Converges on Hamming/len < 0.01, minimum 5 iterations. Says explicitly:
*"No automatic trimming. Trim manually after visual inspection."*

There is **no `-c 70 -g 30`** — I invented that. The real flags:

| script | flags |
|---|---|
| `sine_consensus_smart.sh` | `-n` subsample, `-m` max iters, `-s` stability, **`-t`** call threshold (50), **`-c`** min coverage (30), `-k` keep |
| `sine_pairwise_consensus.sh` | `-r` reference, `-n` max copies, **`-t`** call threshold, **`-c`** min coverage, `-j` jobs, `-N/-M/-S` bootstrap |

So "trim the flanks" is **`-t 70 -c 30`**, not `-c 70 -g 30`.

**Doc/code mismatch in `sine_pairwise_consensus.sh`:** its help says `-t` defaults
to 50 and `-c` to 30; the code sets `CONS_THRESH=30`, `MIN_COVERAGE=10`. Running
it with defaults gives a much looser consensus than the help claims. Pass them
explicitly.
