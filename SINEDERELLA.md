---
name: sinederella
description: What SINEderella is and does - the orchestrator, every step, the standalone tools, what a subfamily is defined as, and which of its tools NOT to reimplement. Read before any SINE search, consensus, boundary, clustering, alignment or report work.
---

# SINEderella

Toki's pipeline. **It already does search, consensus, boundary refinement,
subfamily assignment, subfamily clustering, diagnostic-position analysis,
alignment extraction and report generation.** Read this before writing code for
any of those.

Every statement here comes from reading the script, and every number from
measuring it on a real run. Where his docs and his code disagree, both are
recorded. Facts are stated once — this file was consolidated from four appended
reading passes on 2026-09-03, so there are no "correction to an earlier section"
notes left in it.

**Source of truth:** `/staging/tmp/SINEderella/` on DRAGEN (a git checkout, so
`README.md` and the 70 KB `MANUAL.md` travel with the code); also
`/data/V/toki/bin/SINEderella/` on KIT. The de novo stage is a separate repo,
`/data/V/toki/bin/SINE-de-novo-genome-scan/`. His wider toolbox is
`/data/V/toki/bin/` — about 200 executables.

**Before writing any new script, read `REINVENTED.md`** — 21 of his tools that
I rewrote instead of calling.

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

Distinguishing genuine diagnostic patterns from recurrent/convergent mutation
hotspots (homoplasy) is explicitly parked as future work.

---

## The orchestrator

```bash
SINEderella <GENOME_FASTA> <CONSENSUS_FASTA>       # full run
SINEderella --add|-a <newSINE.fas> [--run <dir>]   # add consensuses
SINEderella --exclude|-e <SINE_name> [--run <dir>] # drop a family, reassign
SINEderella --resume|-r [--run <dir>]              # resume an interrupted run
```

Env: `THREADS` (default nproc), `CHUNK_BP` (30000), `FLANK` (50).

`mode_full` → sanitize genome and consensuses → assembly QC (diagnostic, never
blocks) → **step1 → step2 → step3 → step4 → `results/` symlink farm → step6**.

**Steps 5, 7, 8a, 8b are standalone and are NOT run by the orchestrator.**

| step | script | what it does |
|---|---|---|
| 1 | `step1_search_extract.sh <cons> <genome> [sample=30000] [bin=50]` | `sear -k` per consensus, merge, extract, sample, SubFam, align chunks + consensus bank |
| 2 | `step2_asSINEment.sh <cons> <extracted> <threads> <outdir>` | 10-cycle unanimous voting assignment |
| 3 | `step3_postprocess.sh <rundir> <threads>` | per-copy similarity, LEAK, CONFLICT, soft assignment, summaries |
| 4 | `step4_plots.sh <rundir> [threads] [max_seqs]` | divergence histogram + nucleotide-frequency plots |
| 4d | `step4_diagnostic.sh` → `step4_diagnostic.py` | diagnostic positions, PCA, chimera flags |
| 5 | `step5_align_subfamilies.sh <bank> <multi_dir>` | subfamily alignments — **defective, see below** |
| 6 | `step6_report.sh <rundir> [args]` | `results/report.html` |
| 7 | `step7_boundary_refine.sh <RUN_ROOT> [flank=50] [step=50] [max=1000]` | boundary refinement |
| 8a | `step8a_extract_alignments.sh <RUN_ROOT> <SPECIES_CODE>` | top100/rand100/subfam alignments |
| 8b | `step8b_publish_report.sh <RUN_ROOT> <SPECIES_CODE>` | wires 8a's alignments into `report.html` |

### Things about the orchestrator that bite

- **`@U@` originates here.** `sanitize_fasta` in **genome** mode does
  `gsub(/_/,"@U@",hdr)` (`SINEderella:214`); **consensus** mode does not — it
  only strips after the first whitespace. Any downstream tool that writes
  headers must restore `_`.
- Every mode copies the step scripts into the run dir as a frozen copy, **but
  `--resume` overwrites them** from the current install. Reproducibility does not
  survive a resume.
- **`--resume` never runs step6** when step2+3 are already complete — that branch
  rebuilds the `results/` symlink farm and stops. Only `run_step2_step3` calls
  step6.
- `results/` is a symlink farm rebuilt with `rm -rf results`, which is why step6
  must run *after* it. A bug hit on three runs before being fixed.
- `rebuild_from_searches` (used by `add`/`exclude`) re-merges with
  `bedtools merge -c 4,5,6 -o max,max,distinct` — **no `-s`**, so `+,-` intervals
  are regenerated on every add.
- `add` and `exclude` do **not** copy `assembly_qc.sh` into the new run.
- `find_latest_run` requires a completed run: `consensuses.clean.fa`,
  `genome.clean.fa`, `extracted.fasta`, and an `assignment_full.tsv`.

---

## `sear` — the search engine

`sear [options] <query.fa> <bank.fa> [Slf=0.8] [Homology=65] [Flank=50] [KeepSplits=0]`

**Two engines by query length.** ≥1000 bp → minimap2 (identity =
`matches / min(qaln,taln)`, so gaps are penalised). Shorter — i.e. all SINEs →
ssearch36. `--force-ssearch` overrides.

**Chunking, and why.** `seqkit sliding -g -s 2000 -W 2500` → 2500 bp windows at
2000 bp step, **500 bp overlap**, so a ~300 bp element is never cut without an
intact copy in the neighbouring window; then `seqkit split2 -p 20`. Offsets are
parsed back out of the `_sliding` headers at the end. `ssearch36` has no index
and its statistics degrade on long sequences — that is the reason for chunking.

ssearch36 flags: `-g -3 -T <n> -Q -n -z 11 -E 2 -w 95 -W 70 -m 8C`, ktup 3.

**`Slf` is a minimum hit LENGTH, not query coverage.** `flen = int(Slf ×
query_len)`; a merged hit is kept if its span in the bank is ≥ `flen` **and**
identity ≥ `Homology`. Truncated elements are dropped, not shortened.

**`-m/--minus` is hit depletion, not subsampling.** After a round,
`bedtools complement` of the hits gives the un-hit remainder; the next round
searches that, until a round returns zero. Copies shadowed by stronger
neighbours surface later. The `chrom:start-end()` header form comes from that
`getfasta`.

**`-b/--best N`** → top N unique by bitscore, ±50 bp flanks, `getfasta -s`,
**query prepended as row 1**, `mafft --localpair --maxiterate 1000 --ep 0.123
--nuc --reorder` → `<tax>-<query>.best100.mafft.fa`. Correctly oriented by
construction.

Other outputs: `<tax>-<query>.bed`, `<tax>-<query>.bnk` (hits + flanks,
stranded — the `.bnk` that feeds SubFam). `-s/--stop N` caps hits (default 1000).
`-k` keeps splits for reuse. `--slop L,R` overrides the flanks.

**Usage trap:** `sear` builds its output name from the query **path**
(`qname="${query%.*}"`), so a query in a subdirectory writes to `gen-q/….bed`
and dies with "No such file or directory". Keep the query in the working
directory.

*Doc/code mismatch:* the usage text says the minimap2 cascade is
`asm5 → asm10 → asm20`; the loop runs only `asm10 asm20`.

## `sear_multi`

`sear_multi [-k] [--slop L,R] <multi_query.fa> <bank.fa> [Slf] [Hom] [Flank]`

Same splits and bank sanitisation. **One `ssearch36` call per split with the
whole multi-query FASTA**, output split by query name (column 1), each filtered
by *its own* `flen`. Outputs `<tax>-<queryname>.bed`/`.bnk` per query with `@U@`
restored. **Names come from the FASTA header, not the path**, so it has no
subdirectory trap.

**No `-b/--best`, no `-m/--minus`** — it cannot build top100 alignments.

---

## step1

`sear -k <query> <genome> 0.8 65 50` per consensus — **no `-b`**, so step1 does
not produce top100 alignments.

**It merges hits across all consensuses without `-s`**, so a locus hit by two
families on opposite strands gets strand `+,-`.

**It ends by building the subfam alignment**: `input.clw` → degapped →
`input_reps.fasta` → concatenated with the **original consensus bank** →
`mafft --localpair --maxiterate 1000 --ep 0.123 --preservecase` →
`genome.clean_step1/subfam_input/input.clw.al`.

Timema `run_20260821_132226`: **605 rows = 597 chunk consensuses + 8 family
consensuses**. That is the subfamily alignment with anchor rows, already built.

---

## `SubFam`

`SubFam <bank.fa> [chunk=50]`: `mafft --retree 0 --reorder` → `seqkit split2 -s
50` → align each chunk → **one consensus per chunk** (`cons -plurality 18`) →
`.clw` → `mafft --localpair --maxiterate 1000 --ep 0.123` → `.msf`.

**A "subfamily alignment" is an alignment of chunk-consensus rows, not copies.**

**Chunking is by SIMILARITY, not input order** — `--reorder` emits in guide-tree
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

**`-plurality 18` is hardcoded** — 18/50 = 36 %, correct only because chunks are
50.

**`.clw` is a plain `cat` of the per-chunk `.cons` files**, not an alignment.
MANUAL §6.1.1: degap and realign before measuring anything off it:

```bash
seqkit seq -w 0 input.clw | seqkit seq -g > input.degapped.fasta
mafft --auto --quiet input.degapped.fasta > input.realigned.fasta
```

**Redundancy before SubFam is load-bearing** (§6.1.5): guide-tree order puts a
family's near-identical copies in the same chunk, and `cons -plurality 18` needs
that volume. An upstream CD-HIT-style dedup would starve it.

---

## step2 / `asSINEment` — how a copy gets its subfamily

1. Copies split into blocks of 20 000.
2. **10 cycles** of `ssearch36 … -m 8`, consensuses as query; per cycle each
   copy's best hit by bitscore is one vote.
3. A copy is a candidate only on **10/10 unanimity**.
4. Bar is **relative and per-subfamily**: sort that subfamily's unanimous
   bitscores descending, `N = min(10, count)`, threshold = **0.45 × the N-th
   best**. Below → `rejected_low_bitscore`.

Outputs: `assigned.fasta` (`>seqID|subfamily|bitscore`), `unassigned.fasta`,
`subfamilies/{sf}.fasta`, `assignment_stats.tsv`
(`Subfamily, Assigned, TopN_Bitscore, Threshold`), `assignment_full.tsv`,
`summary.txt`.

**The 10 cycles do real work — tested.** Same command, same inputs, three runs on
20 000 Timema copies: **43 396 / 43 484 / 43 657 rows, three different md5s**,
and **17 copies changed their winning consensus**. `ssearch36` is
non-deterministic under `-T` with `-z 11`, so unanimity is a genuine *stability*
filter.

### Incremental mode

`--incremental --old-assignment <tsv> --new-beds <bed>`:

- **Tier 1** — firmly assigned (10/10) in the old run, still present, **and not
  overlapping the new consensus's BED hits** → carried forward untouched.
- **Tier 2** — everything else → the full 10-cycle vote.

So `--add` does not globally re-open old assignments: a locus re-competes only if
the new family actually hit it. Incremental additionally writes
**`all_votes.tsv`** and appends **`no_unanimous`** rows to `assignment_full.tsv`;
full mode does neither.

---

## step3 — postprocessing

- **`sim_ratio`** per copy = best `ssearch36` bitscore against its own consensus
  ÷ that consensus's self-alignment bitscore, both at `-E 100`. Timema:
  n=182 565, min 0.290, q1 0.599, **median 0.748**, q3 0.897, max 1.028.
- **LEAK** = `runner_ratio ≥ 0.90` — the second-best consensus scores within
  10 % of the best. Timema 9 340 copies. Leak % is a real-family signal: Timema
  families judged real ran 0.00–0.18 %, noisy ones 65–98 %.
- **CONFLICT** = the locus's region carries >1 subfamily label; `conf_alt=yes`
  means the list contains a subfamily **other than** the assigned one.
- **Soft assignment (stage 3.5)** — loci that fail assignment are intersected
  with `all_hits.labeled.bed` to recover which query found them; the
  highest-scoring query becomes `Soft_Subfamily`, all candidates are listed in
  `All_Queries`, with a `Reason`. Written to `unassigned.tsv`. Nothing is
  discarded.

`summary.by_subfam.tsv` columns: `subfam, firm_assigned, soft_assigned,
total_assigned, leak_n, conf_alt_n, firm_pct, total_pct, leak_pct,
conf_alt_pct, sim_mean, sim_median`. `firm_pct`/`total_pct` are fractions of
**total extracted**; `leak_pct`/`conf_alt_pct` of **firm assigned**.

Sanity counters printed at the end: `runner_ratio_tags`, `LEAK_rows`,
`CONFLICT_rows`, `conf_alt_yes_any`, `sim_ratio_scored`,
`unassigned_soft_assigned`.

---

## The `(+,-)` strand — measured

Timema `assigned.fasta`, 337 898 loci: `+` 128 115, `-` 128 908, **`+,-` 80 875
= 23.9 %**. step8a does `if (strand != "+" && strand != "-") strand = "+"`.

Cost, measured as identity to the consensus row: `top100` 0 flipped rows in 7 of
8 subfamilies (t4: 13 low, 4 genuinely reverse); `rand100` 0 in 5 of 8. **A real
defect worth 0–4 rows per 100.**

**`+,-` is the CONFLICT set.** 82 365 of 339 566 regions (24.3 %) carry more than
one subfamily label — the same population. The **800 highest-bitscore loci are
all `+,-`**: the most complete copies are the ones several consensuses hit.

The dominant multi-label set is **`t3,t4,t5` at 69 585 loci**.

`regions.by_subfam.bed` is **not** broken in this run (contrary to the
`sine-orth-loc-sinederella` skill): 170 186 lines list `t1` alone.

**He has already found and fixed this bug elsewhere.** From
`extract_top100_rand100_subfam.sh`'s header: `(+,-)` left un-stripped "contains a
literal `-` that corrupted the coordinate split below (silently dropping the BED
line, which meant these hits ended up unflanked in the alignments with no error
— **found by direct verification, not by trusting the success-looking logs of
the pipeline itself**)".

---

## step4 — plots and diagnostics

### `step4_plots.sh`

Per subfamily, subsampled to 10 000 with `seqkit sample -s 42`:

1. `ssearch36 -E 100 -m 8` consensus vs copies, best hit per copy, **column 3
   (% identity)** → divergence histogram.
2. `mafft --addfragments <copies> --adjustdirection` against the consensus alone
   — pairwise against a reference, **not** a full MSA → per-position nucleotide
   frequency.

`plot_subfamily.py` writes `{sf}_divergence.png/pdf`, `{sf}_nucfreq.png/pdf`.

### `step4_diagnostic.py` — the diagnostic-position machinery

Nine stages:

1. **Bitscore matrix** — one `ssearch36 -E 100 -m 8`, all consensuses vs all
   assigned copies → copies × subfamilies. Falls back to `assignment_full.tsv`.
2. **PCA** on that matrix, with loading arrows; writes `pca_coordinates.tsv`.
3. **Position weights**:

| method | what it captures |
|---|---|
| **Mutual information** `I(Position; Subfamily)` | how much a column tells you about the label |
| **Random Forest importance** (`n_estimators=200, max_depth=5`) | **interactions between positions** |
| **KL divergence** | pairwise discrimination between subfamilies |

`combined = 0.4·MI + 0.3·KL + 0.3·RF`, top-K (default 30) →
`diagnostic_positions.tsv`, `position_weights.tsv`, plus `synergy_pairs.tsv`.

4. **Per-copy diagnostic state** — copies aligned to the longest consensus with
   `mafft --add`; nucleotide read at each diagnostic position.
5. **Scoring and discordance flags:**

| flag | condition |
|---|---|
| `fp_risk` | assigned, diagnostic score < 0.5 with a higher-scoring runner-up; or < 0.3 outright |
| `fn_risk` | **not** assigned, but diagnostic score > 0.7 |
| `chimera_suspected` | 5′ half of the diagnostic positions best-matches one subfamily, 3′ half another, **each at ≥60 % of that half** |

6. copies × positions heatmap. 7. diagnostic-score vs `sim_ratio` scatter with
discordance quadrants. 8. volcano plots per subfamily pair. 9.
`sineplot_input.txt` — all-vs-all `ssearch36 -m 8` for `SINEplot.py`.

**`chimera_suspected` is the working mosaic detector.** It splits the diagnostic
pattern rather than averaging a similarity, which is why it works where
`verdict.py`'s `homogeneity` does not — that one reads 1.000 for 72 % of corpus
sets, including 8/8 NEGMOSAIC, 6/6 SIMMOSAIC and 6/6 SIMSCRAM, identical to
28/28 POS.

**This is MANUAL §6.1.6's definition, implemented.** It is **supervised** — MI is
computed against known labels — so it composes with an unsupervised peel that
proposes a partition: peel to propose, `step4_diagnostic.py` to weight, validate
and flag chimeras. Do not reimplement position weighting.

---

## step5 — defective, prefer step8a

- Large subfamilies (≥400) get `mafft --add <bank> --keeplength` — **`--keeplength`
  forces the added consensus into the existing column frame, deleting consensus
  residues that would need new columns.**
- Small subfamilies get `--auto`, not `--localpair --maxiterate 1000 --ep 0.123`.
- `ALIGNMENT_DEPLOYMENT.md` also notes it produces no flanks.

`step5_direct.sh <run_dir>` is the single-run version (same two branches, output
to `<run>/alignments/`). `run_step5_wrapper.sh` fakes the multi-species layout
step5 expects with a symlink in a temp dir, then moves results back.

---

## step6 — the report, and the page template

`step6_report.sh <RUN_ROOT> [args]` → `step6_report.py`. Picks the newest
python ≥3.8 on PATH. Key flags: `--out`, `--inline-plotly/--no-inline-plotly`,
`--no-embed-images`, `--max-rows`, `--sineplot/--no-sineplot`, `--sineplot-max`,
`--threads`, **`--tal-species-code <code>`**.

Sections: `overview`, `composition`, `divergence`, `pca`, then collapsible
`run`, `step1`, `assignment`, `thresholds`, `flags`, plus the step4 gallery.

Design system is in the `CSS` constant: `--accent:#4C72B0`, card sections with
8 px radius, `.metric` tiles, `.tbl` tables with hover-help headers, `.aln-link`
pill buttons (blue / orange / green), a lightbox for plot images, and
`details.card` collapsibles.

**Two different PCAs exist.** `build_pca_fig` here is a **mutation-space** PCA:
each copy is a binary vector — 1 where it differs from its assigned subfamily's
consensus — via `mafft --add --keeplength`, then SVD. Stratified sample up to 200
per subfamily, seed 42. `run_pca` in `step4_diagnostic.py` is PCA on the
**bitscore** matrix.

**`--tal-species-code <code>` builds its own alignments section** and links, per
subfamily: `{code}_{sf}_top100.aln.fa`, `{code}_{sf}_rand100.aln.fa`, and
**`{sf}.al` from `<base>/<code>/subfam/`** — the manual subfam alignments — plus
section-level links to `{code}_consensuses.fa` and `{code}_subfam_input.aln.fa`.
The section text states the flanks as **50 bp upstream + 70 bp downstream**.

Two cautions: the section links `{sf}.al` for **every** subfamily in
`assignment_stats.tsv`, but `.al` files only exist for subfamilies with **≥400
copies** and only for whatever partition the generator last ran on — regenerating
a report against a newer partition links files that do not exist. And the
cross-species nav bar is **hardcoded to `saq ↔ ccr`**; every other species gets
an empty link.

Charts are emitted as Plotly **dicts built by hand** — the plotly Python library
is never imported. `kde_curve` is a hand-rolled Gaussian KDE with Scott's
bandwidth `h = std · n^(-1/5)`, no scipy. `conservation_curve` plots the
per-position frequency of the consensus base. `stratified_sample_sim` samples up
to 3000 copies per subfamily for the density plots. `fig_funnel` shows
hits → unanimous 10/10 → assigned → unassigned.

### Two defects in the report, verified against data

**The LEAK column is always zero.** `count_flags_per_subfam` reads `parts[8]`
for both flags, but step3's ALL table puts `leak_flag` in column 8 (`parts[7]`)
and `conflict_flag` in column 9 (`parts[8]`), so the `elif "LEAK" in flag` branch
tests the conflict field and never fires. Measured on Timema
`run_20260821_132226`: the file holds **9 340 LEAK rows**, and the parser returns
`CONFLICT=82 365, LEAK=0, OK=257 201`. `flags_table` then renders a LEAK column
of zeros. That hides a signal his own notes call diagnostic — real Timema
families run 0.00–0.18 % leak, noisy ones 65–98 %.

**The "Threshold/Self" column is not on a usable scale.** `thresholds_table`
computes `(threshold / 100.0) / real_self_bits`, but `asSINEment` writes the
threshold as a plain bitscore (`0.45 × the N-th best`), not ×100. On Timema:
t1 `threshold=782`, `self_bits_real=179.1` → the column shows **0.0437**, while
`threshold/self` is **4.366**. The two are on different scales whichever way it
is divided — thresholds come out several times larger than the recorded
consensus self-bitscore. Ask which of the two numbers is the unexpected one
before changing anything.

## `plot_subfamily.py`, `analyze_convergence.sh`, `benchmark_sear.sh`

**`plot_subfamily.py`** — the step4 plotting helper. Divergence histogram of
`100 − %identity` with mean and median lines; per-position stacked nucleotide
frequency using the **consensus row's non-gap columns as the position axis**, so
the x-axis is consensus coordinates and each tick is the consensus base itself.
Everything non-ACGT counts as `gap`. It strips mafft's `_R_` prefix from names
(added by `--adjustdirection`).

**`analyze_convergence.sh <logdir>`** — batch-summarises `*_consensus.log` into
a CSV (`NAME,SEQS,ITERS,BP,STATUS,AVG_CHANGE`) and flags problems:
`ITERS > 30` = slow convergence, likely mixed or FP-contaminated data;
`NOT_CONVERGED` = hit the 50-iteration cap; `AVG_CHANGE > 0.01` = noisy.

**`benchmark_sear.sh <genomes.tsv> <consensus.fa> [max_mb] [n_queries]`** —
he has already benchmarked `sear` per query against `sear_multi`: truncates
genomes to a size cap, picks N queries, runs both approaches, and compares wall
time, hit counts and BED identity. Read this before re-deriving which is faster.

## SINE-KB — the downstream knowledge base

`SINEderella_multi` calls `generate_sinekb_report.py`, and
`import_squamata_run.py` loads the resulting `sinekb_report.json` into a system
under `sine-kb/` via `models.import_sinekb_report()`: per-genome copy counts, a
taxonomy tree extended per import, and a rebuilt consensus bank. The squamata
import covers **79 genomes and 9 subfamilies**. Neither `sine-kb/models.py` nor
`generate_sinekb_report.py` has been read.

## step8b — the other way alignments reach a report

Reads `results/alignments/manifest.tsv` (**not** a directory listing), builds
`<section id="alignments">` with `top100`/`rand100`/`subfam` columns, and splices
it at the placeholder `"alignments not yet available"` — **or, if that is absent,
appends before `</body>`, so a second run duplicates the section.**

It has **no slot for `{sf}.al`**. `--raw-base-url` turns links into MSA-viewer
URLs; without it they stay local relative paths.

---

## step7 — boundary refinement. Do not reimplement.

**Background is measured, not assumed**: 300 random 70 bp genomic regions, 2000
random pairs, ungapped positional identity; `elevated_frac` = fraction of pairs
above `ELEVATED_CUTOFF=45 %`; pass threshold = `background_frac + 5`.

**The walk** tests **one `STEP_BP`-wide window sitting `ext` bp out** from the
element edge — not a growing region — strand-aware (upstream/downstream swap for
minus copies), reverse-complemented so every window reads in element
orientation. `ext` starts at `FLANK_BASE`, grows by `STEP_BP` to `MAX_EXT_BP`.

Stop conditions: `elevated_frac ≤ threshold` → `confirmed`; running off the
contig → also `confirmed`; exhausting `MAX_EXT_BP` → `undetermined` with
`boundary_bp = MAX_EXT_BP`; fewer than 20 members → `insufficient_data`.

**Populations are separate**: **general** = random sample capped at 2000 (seed
42); **top100** = `sort -k6,6nr | head -100`, matching step8a's actual sample
rather than approximating it. His note: a random-100 sample confirmed background
at ext=400 on scorpion g01 while the real top100 stayed ~90 % identical in that
same window, because top100 members are the least-diverged copies.

Output `boundary_refinement.tsv`: `subfamily | population | side | boundary_bp |
status | final_identity_pct | final_elevated_frac_pct | background_identity_pct |
background_elevated_frac_pct`.

**Known defect: `background_identity_pct` is deflated by soft-masking.** The
comparison is case-sensitive (`a[i] != b[i]`). Timema's genome is 48.4 %
soft-masked, giving **13.640 %** — reproducing his output exactly — where
case-insensitive gives **27.176 %**. The decision variable is `elevated_frac`
(~0 either way against a 5 % threshold), so **the boundary calls are unaffected**;
only the reported number is wrong. One `.upper()` fixes it.

---

## The alignment products

### step8a

Base flanks **50 up / 70 down**, extended per subfamily from step7 but **only on
`status=="confirmed"` rows** — `undetermined` means the walk hit the cap, and
applying it once produced a 1070 bp flank that made mafft crawl. `top100` uses
step7's top100-population rows; `rand100` and `subfam` use the general rows.
`subfam` only at ≥400 members (seeded shuffle 42, up to 10 000). `--reorder` can
move the consensus off the top, so a helper forces it back to row 1; `@U@`
restored on write.

### `extract_top100_rand100_subfam.sh` — the species-page generator

```
extract_top100_rand100_subfam.sh <RUN_ROOT> <SPECIES_CODE> <OUT_DIR>
```

`FLANK_L=50 FLANK_R=70 TOPN=100 RANDN=100 SEED=42 BINSIZE=50 N_SAMPLE=10000
MIN_COPIES_SUBFAM=400`.

| output | how |
|---|---|
| `{sp}_{sf}_top100.aln.fa` | bitscore from the `>ctg:start-end(strand)\|sf\|bits` header, top 100, **re-extracted from the genome with 50 L + 70 R strand-aware flanks**, consensus appended, mafft, `postprocess_flanks`, `@U@`→`_` |
| `{sp}_{sf}_rand100.aln.fa` | same, sampled with his `/data/V/toki/bin/sample` (unseeded shuf; inline fallback elsewhere) |
| `{sp}_{sf}_subfam.aln.fa` | **≥400 copies only**: seed-42 sample up to 10 000, `SubFam input.fasta 50`, `.clw` degapped, consensus appended, mafft. **Not flanked** |

`MAFFT_OPTS = --thread N --localpair --maxiterate 1000 --ep 0.123 --nuc
--reorder --preservecase --quiet`. The consensus record is renamed
`{sf}_CONSENSUS`, which is what `postprocess_flanks` matches on.

Its header also records a SIGPIPE trap: `shuf | head` under `set -o pipefail`
kills the script — write intermediates to files.

### `run_subfam_per_sf.sh` — where `subfam/*.al` comes from

Per subfamily: sample up to 10 000 from `assigned.fasta`, **skip below 400
copies**, `SubFam input.fasta 50`, degap `.clw`, `cat input_reps.fasta {sf}.cons`
→ mafft → `{sf}.al` in `<run>/results/subfam_alignments`. These files are
pipeline output, not manual work.

### `extract_alignments.sh` — four tiers

| tier | file | content |
|---|---|---|
| A | `{sf}.core.fa` | random 50, **30 L + 70 R** |
| B | `{sf}.best50.fa` | top 50 by bitscore, same flanks |
| C | `{sf}.subfam.fa` | **>400 copies**, subsample **2500**, `SubFam`, then `cat consensus.fa input.clw` — **`.clw` NOT degapped** |
| D | `{sf}.evidence.txt` | **zero copies**: `sear <cons> <genome> 0.5 30 20` (relaxed) in a shared work dir, first 10 BED hits — so an absent family is recorded as *tested* |

**It pools firm + soft copies**, soft ones re-headed `>seqid|subfam|0|SOFT` from
`unassigned.fasta` via `Soft_Subfamily`. `extract_top100_rand100_subfam.sh` uses
`assigned.fasta` only — the two generators do not sample the same population.

`extract_subfam_only.sh` is Tier C alone (>400 copies, sample 2500, `.clw` not
degapped), resumable, for runs where A/B finished without SubFam available.

### `postprocess_flanks()` — the flank convention

The consensus row's first and last non-gap column define the body. Outside it,
per copy: **internal gaps deleted, bases lowercased, gaps padded on the outer
side** — left flank right-justified against the body, right flank left-justified.
Column count preserved, column correspondence gone.

This is the "no ragged tails" convention; `ALIGNMENT_DEPLOYMENT.md` confirms the
viewer "must show sequences WITH lowercase flanking DNA before/after the SINE".

**Therefore any column-wise flank statistic — island z-scores, flank identity
decay — is meaningless on `.core.fa`/`.best50.fa`.** It is valid on step8a's
`*_top100.aln.fa` / `*_rand100.aln.fa`, which do not get this treatment.

---

## Consensus building

**`conse <msa.fasta> [pct]`** — EMBOSS `cons` at a plurality of **pct % of the
number of sequences, default 35 %**, then **N → `-`**, one-line FASTA. Uncalled
columns become gaps, so the boundary emerges by itself. **A consensus threshold
is a fraction of the sequences, never a count.**

**`sine_consensus.sh input.fasta [subsample=100] [max_iters=50] [thresh=0.01]`**
— **positional only, no flags.** Bootstrap: subsample → `mafft --auto` →
mini-consensus → append to master → realign → converge on Hamming/len < 0.01,
minimum 5 iterations. Column rule: **>50 % gaps → `-`**, otherwise plurality
among non-gap bases, ties → `N`. *"No automatic trimming. Trim manually after
visual inspection."*

**`sine_consensus_smart.sh [options] input.fasta`** — same bootstrap plus
variance detection: keeps the last `VARIANCE_WINDOW=3` change values and, past
`MIN_ITERS=10`, stops early with `QUALITY: LOW_CONFIDENCE` if their mean exceeds
`VARIANCE_LIMIT=0.15` ("likely garbage/FP-contaminated data"). Uses
`esl-alistat` for per-iteration average identity when present. `-n` subsample,
`-m` max iters, `-s` stability, **`-t`** call threshold (50), **`-c`** min
coverage (30), `-k` keep trace.

**`sine_pairwise_consensus.sh [options] input.fasta`** — bootstrap a reference,
then `mafft --auto --op 5` each copy against the reference **alone** and project
onto reference columns; gap columns in the reference are recorded as insertions
`refpos:length`. Every projection is exactly `REF_LEN` long and this is checked.
Call per position: coverage gate `nongap ≥ n × MIN_COVERAGE`, then plurality with
`freq = max_count / n` — **gaps in the denominator** — kept if `freq ≥
CONS_THRESH`; ties → `N`. Writes `freqtable.tsv`
(`Pos Ref A C G T Gap Nongap Total Called Freq`) and an **insertion hotspot**
report for positions where >10 % of copies insert.

*Doc/code mismatch:* its help says `-t` defaults to 50 and `-c` to 30; the code
sets `CONS_THRESH=30`, `MIN_COVERAGE=10`. Pass them explicitly.

**"Trim the flanks" is `-t 70 -c 30`.** There is no `-c 70 -g 30`.

His rule on which: **approximate consensus → the bootstrap/pairwise method;
robust and data-supported → `conse` over the top 100 by bitscore.** Run that
final pass **once** — iterating creeps outward (203 → 252 → 269 bp on
hyd_SINE_0).

**Do not bolt on an automated consensus-refinement loop** (§6.1.5, explicit):
`Refiner`-style `refineUntil` or TEstrainer BEAT loops are unreliable alone —
they still need human judgement to catch a consensus that "converged"
confidently onto a chimeric or boundary-contaminated alignment.

---

## The manual subfamily step (§6.1) and its tools

1. Degap + realign `input.clw` (§6.1.1 — never cluster off the raw `.clw`).
2. Open `input.realigned.fasta` in MSA-viewer, identify groups sharing
   diagnostic columns and a consistent indel pattern (§6.1.2).
3. Extract each candidate group to its own FASTA.
4. Build a consensus per group (§6.1.3).
5. Re-run SINEderella with the new bank, or `step2_asSINEment.sh` directly
   against the existing `extracted.fasta` (§6.1.4).

**`step1b_cluster_subfamilies_assist.sh input.clw out_dir [--recurse]`** runs
`cluster_assist.js` → `subfam_cluster_lib.js`'s `SINEClusterer`. **An assist, not
a replacement**, with a measured limit: validated 2026-08-21 against Tal/saq's
9 manual subfamilies over 600 chunk consensuses, it recovered a true 2-subfamily
split correctly but collapsed the 9-subfamily case into **2 blobs (293 + 307)
plus one 6-seq cluster**; `--recurse` recovered nothing further.

### How `SINEClusterer` works — the peel loop, already written

- A candidate group is the **exact set of sequences sharing one (position,
  character)**. Gaps are never diagnostic. Columns where one base holds >80 % of
  the *whole* alignment are skipped (checked globally, not against the shrinking
  pool).
- Candidate sets fuzzy-merged at Jaccard ≥0.90 with size difference ≤5.
- Group size capped at **50 % of available**, with a relaxed-upper-bound retry.
- Features scored: `outside == 0` → **3** if every member matches, **2** at ≥80 %,
  **1.5** otherwise; else `qual = inside% − outside%` above a size-dependent
  threshold → **1**.
- Outliers pruned: members matching fewer than ~30 % of the group's features are
  dropped, then the group is re-validated.
- **Peel**: take the best group, remove it, repeat up to 20 iterations, with
  thresholds **relaxing** each round (`minPerfect` 4 → 1, quality down by up to
  50 points) and an ultra-relaxed rescue at ≤10 remaining.

### `getTrimBoundaries` — a second boundary tool

Sliding window of **15 columns** from each end; while the windowed gap fraction
exceeds the threshold, keep trimming. **Thresholds are asymmetric: 0.50 left,
0.80 right** — the 3′ end tolerates far more gap, which is what a poly-A tail
needs.

---

## SINE-de-novo-genome-scan — the stage before SINEderella

Finds SINEs with **no library**, on the premise that most SINEs share at least
faint similarity with some known SINE.

**1. Shred a SINE database** (`01_sanitize_shred.sh`,
`02_filter_lc_cluster_nr85.sh`). Deduplicated at 95 % with vsearch, then cut
region-aware because head and tail are diagnostic:

| region | span | step |
|---|---|---|
| 5′ | bases 1–150 | 10 bp |
| middle | 151 to L−100, only if L > 250 | 25 bp |
| 3′ | last 99 bp, only if L > 150 | 25 bp |

Headers `>SINE_ID|REGION:START-END`; IDs sanitised. Then low-complexity filtered
and clustered to `nr85`.

**2. Scan** (`sine_scan.sh`): per-fragment Smith-Waterman, filter on identity and
**query** coverage, merge hits within 500 bp into candidate loci, extract with
flanks → `query_summary.tsv`, `all_hits.filtered.m8`, `candidate_loci.bed`,
`candidates.fa` → feeds SubFam.

Two versions: the newer takes `-q -g -o` and parallelises per query; the scorpio
copy takes positional args and chunks the database into 100 Mb blocks.

**Consequence:** these loci are **anchors expanded to candidates**, not verified
elements. A de novo family scoring FRAGMENT_OF_LONGER is describing its input
honestly.

---

## `assembly_qc.sh`

Single-pass awk. Metrics: scaffolds, total bases, ACGT, N, `gap_pct`, N50, L50,
largest scaffold, `short_scaffolds_pct` (<1000 bp), and **`edge_N_pct`** — N
fraction in the first and last **500 bp** of every scaffold, which its comment
calls "the *Anilios bituberculatus* signature" for padding artefacts.

Flags: N50 ≥10 Mb ok / ≥1 Mb caution; gap <5 % ok / <15 % caution; edge N <5 % ok
/ <20 % caution; short scaffolds <5 % ok / <20 % caution. Verdict: any warn →
**WARN**; ≥2 cautions → **CAUTION**; else **SOLID** (the orchestrator's message
text expects `OK`).

## `SINEderella_multi`

Subcommands `full | add | exclude | plots | summary`. `genomes.tsv` is
tab-separated `Species_Name  /genome.fa  [consensus.fa]  [workdir]`, `#`
comments. Options `--project`, `--parallel N`, `--threads N`, `--resume`.

Output `multi_YYYYMMDD_HHMMSS/` with `project_manifest.tsv`, per-species dirs,
and `summary/cross_species_summary.tsv` — a wide table `Species | Total_Loci |
<sf>_firm | <sf>_total … | Unassigned_n | Unassigned_pct`. It detects old vs new
`summary.by_subfam.tsv` layouts by checking whether column 2 is `firm_assigned`.
Also calls `generate_sinekb_report.py` when present.

---

## Two competing models of divergence (RESEARCH_DIRECTIONS.md)

Everything in the pipeline assumes a **gradual, substitution-driven** model:
relatedness measured as identity. The competing model is **modular / "repeat
pangenome" (panconsensus)**: some repeats are composites from a library of
reusable, independently-mobile modules, so relatedness is a **presence/absence
and arrangement matrix**, not one number.

Evidence in hand: the `eri` e2-3 downstream flank is **bimodal** — about half of
sampled copies still share an exact motif 100–120 bp past the called SINE end
while the rest look like background at the same offset. The gradual framing
called that boundary under-calling and needed a mean → fraction-of-pairs fix to
detect at all.

The two predict different signatures: gradual under-calling correlates boundary
distance with subfamily age; a present/absent module gives a clean two-cluster
split with no gradient. Proposed as an **additive standalone step** (segment →
cluster blocks not copies → presence matrix → compare model fit), to be tried
first on `eri` e2-3. Idea only as of 2026-08-21, no code.

(`UnitFinder.sh` — "represent locus sequence as a graphical scheme of its parts",
using a pseudoconsensus that masks structural regions — may already be relevant
infrastructure.)

---

## Reading a run

`results/`: `assignment_full.tsv` (Sequence, Subfamily, Bitscore, Votes, Status,
Threshold, RunnerInfo), `summary.by_subfam.tsv`, `subfamilies/<sf>.fasta`
(element-only, headers `>chrom:start-end(strand)|subfam|bitscore`, `_` escaped as
`@U@`), `assigned.fasta`, `unassigned.fasta`, `unassigned.tsv`, `all_votes.tsv`,
`sim_scores.tsv`, `self_bits_real.tsv`, `consensuses.fa`, `plots/`.

- **No run has `manifest.txt` in the shape a `run_2026*/manifest.txt` glob
  expects** — assignment output is under `step2/step2_output/`.
- Cluster IDs `C<N>R` are **not stable** across runs.

## Environment

Needs `ssearch36`, `mafft`, `bedtools`, `samtools`, `seqkit`, `cons` (EMBOSS),
Python 3 + matplotlib + sklearn + pandas. On DRAGEN:

```bash
export PATH=/staging/conda/envs/bioinfo/bin:/staging/miniconda3/bin:/usr/local/bin:$PATH
```

`node` is **not** installed on DRAGEN or KIT, so `cluster_assist.js` cannot run
there as-is. `SubFam` has a `/bin/sh` shebang but bash syntax — on KIT run it as
`bash /data/V/toki/bin/SubFam`. Never write run data to `/home` on DRAGEN.

Run locations are in `DATA_LOCATIONS.md`. See also the
`sine-orth-loc-sinederella` skill for `SINE_orth_loc.bash`, PM/MP/SINE orthology,
and the three non-interchangeable subfamily labelling systems.

---

## Measured facts worth keeping

**The `sear` rebuild works.** saq, all 9 consensuses through `sear --slop 50,70
-b 100`, ~380 000 hits each:

| | blastn + 400 bp flanks | `sear --slop 50,70` |
|---|---|---|
| copy length (median) | 1054 bp | **366–381 bp** |
| `conse` at 35 % calls | **541 bases** | **281–321** |
| gaps in the consensus body | 47–73 % | **13–23 %** |

s8: 541 bases at 73 % → **307 at 14 %**, against a ~253 bp element.

**Peeling on diagnostic columns beats identity.** Timema, 596 chunk consensuses,
against his 8 groups: identity-based placed 111 with t3 never recovered;
feature-based placed **503 at 0.958 weighted purity**, recovering t1, t2, t3, t6,
t7, t8, with the re-alignment step firing twice. Against the v4 13-group
partition: **493 of 595 at 0.826**, ceiling 0.882 (the chunk labelling's own
purity), and **all four t1 subfamilies separate** — `t1_1`, `t1-2`, `t1-3`,
`t1-4`. Failures: `t3-1`/`t3-2` never separate (one 119-chunk cluster at 58 %) —
and his own `tim/subfam/` contains **`t345.al`**, so he groups t3+t4+t5 together
too.

---

## Traps, in one place

1. **Do not cluster subfamilies on pairwise identity** — §6.1.6 says that is the
   noise, not the signal.
2. **Do not measure column-wise flank statistics on `.core.fa` / `.best50.fa`** —
   `postprocess_flanks` has destroyed column correspondence there.
3. **Do not orient subfam chunk consensuses.** They inherit the input
   consensuses' orientation; flipping fused Timema t3 with t4 and split t1 in
   two. Orientation *is* needed for AnnoSINE seed candidates.
4. **Do not use single linkage.** It chained t3/t4/t5 into one 143-chunk cluster
   whose median within-cluster identity was 0.983 — the highest of any — while
   being 40 % pure.
5. **Do not treat `input.clw` as an alignment.** Degap and realign first.
6. **Do not dedup before SubFam.** The redundancy is what `-plurality 18` needs.
7. **Do not iterate the final `conse` pass.** Once.
8. **Do not apply an `undetermined` step7 boundary.** Fall back to the base flank.
9. **Check the step7 population** before reusing a boundary: top100 ≠ general.
10. **Keep `sear`'s query in the working directory** — it names outputs from the
    query path.
11. **Do not run step8b twice** — it appends before `</body>` when its
    placeholder is gone.
12. **Restore `@U@` → `_`** in anything you write out.
13. **Check the read/unread ledger below before claiming to know a script.**

---

## Read/unread ledger

**Read in full:** `SINEderella` (1036), `step1_search_extract.sh` (278),
`step1b` (50), `sear` (653), `sear_multi` (416), `SubFam` (58), `conse` (8),
`asSINEment` (440), `step2_asSINEment.sh` (900), `step3_postprocess.sh` (591),
`step4_plots.sh` (197), `step4_diagnostic.py` (1190),
`step5_align_subfamilies.sh` (231), `step5_direct.sh` (150),
`run_step5_wrapper.sh` (78), `step6_report.sh` (37),
`step7_boundary_refine.sh` (332), `step8a_extract_alignments.sh` (410),
`step8b_publish_report.sh` (206), `extract_alignments.sh` (577),
`extract_subfam_only.sh` (161), `extract_top100_rand100_subfam.sh` (462),
`run_subfam_per_sf.sh` (188), `assembly_qc.sh` (197),
`SINEderella_multi` (659, structure + summary builder),
`cluster_assist.js` (132), `subfam_cluster_lib.js` (504),
`sine_consensus.sh` (185), `sine_consensus_smart.sh` (288),
`sine_pairwise_consensus.sh` (686), `sine_scan.sh` (373),
`step6_report.py` — section structure, CSS, alignment section, both PCA builders,
SINEplot iframe.

**Also read:** `plot_subfamily.py` (283), `analyze_convergence.sh` (62),
`benchmark_sear.sh` (294), `import_squamata_run.py` (325), and the whole of
`step6_report.py`'s parsing, table and figure logic.

**Not read:** `generate_sinekb_report.py` and `sine-kb/models.py` (the
knowledge-base import path), `SINEplot.py`, `step4_diagnostic.sh` (95, a thin
wrapper around the python), and the MANUAL's 70 KB beyond §6 and §7.

**Do not describe anything in the second list as understood.**
