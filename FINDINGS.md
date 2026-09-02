# Findings — what his judgements have established

Separate from `HANDOFF.md` (working history) and `PLAN.md` (what to do next).
This file holds only **results that change the tool**, newest first. Updated
after every exchange.

*(Build-1 findings from 2026-08-31 moved to `FINDINGS_old_2026-08.md`.)*

---

## 2026-09-02 — the peel loop, run on Timema against his own partition

Input: 597 SubFam chunk consensuses from `run_20260821_132226` on DRAGEN.
Truth: his 8 groups, assigned per chunk by majority vote of its 50 members
against the run's own `assignment_full.tsv`. Mean chunk purity **0.953**, so a
chunk almost always belongs cleanly to one group — the chunk-level abstraction is
sound.

### It recovers five of his eight groups cleanly, in his order

Ranked by isolation — within-group identity minus the best identity to anything
outside, which is his "cleanly isolated from alignment":

| isolation | chunks | his group | purity |
|---|---|---|---|
| +0.414 | 30 | **t2** | **100 %** |
| +0.409 | 69 | t3 | 45 % |
| +0.234 | 16 | **t6** | **100 %** |
| +0.221 | 74 | t4 | 54 % |
| +0.161 | 55 | **t8** | **100 %** |
| +0.135 | 297 | **t1** | **95 %** (his t1 = 288) |
| +0.117 | 35 | **t7** | **97 %** |

Residue 20 chunks. So `t1, t2, t6, t7, t8` come out at 95–100 % pure, and only
the `t3`/`t4` pair is wrong.

### His correction, which changed the result

> subfam after sinederella usually doesnt need flip because consensuses used are
> expected to be properly oriantated

Correct, and my orientation pass was doing real damage. Compared directly:

| | with my flipping | as SINEderella left them |
|---|---|---|
| `t1` | split into 210 + 89 | **one cluster, 297, 95 % pure** |
| `t3` / `t4` | **fused**, one 143-chunk blob at 40 % | **separated**, 69 and 74 |

Flipping fused t3 with t4 and split t1 in two. **Do not orient subfam chunk
consensuses** — they inherit the orientation of the input consensuses.

Worth keeping separate: orientation *is* needed for AnnoSINE seed candidates,
where about half are reverse-complemented. It is not needed after SINEderella.
Two different situations, one of which I applied to the other.

### Single linkage chains; average linkage does not

Single-link produced the 143-chunk blob whose **median within-cluster identity
was 0.983** — the highest of any cluster — while being 40 % pure. A chain of
groups linked through intermediates has high pairwise identity at every link, so
no within-cluster statistic can see the join. Average linkage cannot make that
merge.

### Where it disagrees with him — and the evidence

Mean pairwise identity between his groups, at chunk level, random-sampled:

|  | t1 | t2 | t3 | t4 | t5 | t6 | t7 | t8 |
|---|---|---|---|---|---|---|---|---|
| **t1** | **0.878** | 0.467 | 0.448 | 0.421 | 0.481 | 0.447 | 0.473 | 0.448 |
| **t2** | | **0.950** | 0.475 | 0.398 | 0.550 | 0.350 | 0.369 | 0.410 |
| **t3** | | | **0.685** | 0.685 | **0.727** | 0.393 | 0.401 | 0.413 |
| **t4** | | | 0.695 | **0.920** | 0.466 | 0.418 | 0.392 | 0.402 |
| **t5** | | | | | **0.979** | 0.364 | 0.413 | 0.420 |
| **t6** | | | | | | **0.846** | 0.594 | 0.594 |
| **t7** | | | | | | | **0.874** | 0.665 |
| **t8** | | | | | | | | **0.876** |

Every group is internally tight — 0.85 to 0.98 — **except t3, at 0.685**. And
t3's members are on average *more* similar to t5 (0.727) than to each other, and
exactly as similar to t4 (0.685) as to each other.

So the reason the loop cannot recover t3 is not that the loop is blind: **t3 as
labelled does not behave like one group at chunk level.** That is a question for
him, offered as a question — he said he may have made errors, and this is the
one place the measurement disagrees. It is equally possible that what separates
t3 from t4 lives in the members and is smoothed away by the consensi-of-50s,
which is his own reason for building group consensuses from members.

Also visible: t6, t7 and t8 are a related triplet (0.59–0.67 between them) while
each is internally distinct, and t2 and t5 are the tightest things in the set.

### What is still not exercised

Everything peels in round one, so the re-alignment he described — the reason for
peeling at all — never gets tested. The intrinsic defer test (re-cluster a
candidate at a stricter threshold and hold it back if it splits) fires on
nothing, because t3's cluster is homogeneous at 0.83 despite being impure. A
better hold-back rule is needed before the re-alignment effect can be measured.

---

## 2026-09-02 — two discovery routes, one genome, opposite answers

*Centruroides vittatus* has now been through both routes, so they can be compared
on identical DNA.

| route | candidates | scored | at 50+ | dominant reason |
|---|---|---|---|---|
| **AnnoSINE_v2** (structure: TSD, poly-A, pol-III boxes) | 15 | 12 | **7 clean at 100.0** | none - six are flagless |
| **his de novo fragment scan** (SINEBase fragments, 65 % over 50 bp) | 47 | 19 | 3 | **FRAGMENT_OF_LONGER, 16 of 19** |

The structural route gives `sco_SINE_1, 2, 3, 4, 6, 7` at **100.0 with no flags
at all** - identity 0.69-0.91, 74-98 of 100 copies in the core, flank background
0.27-0.31, which is the signature of independent insertions.

The fragment route gives the opposite signature on the same genome: 16 of 19
families are read as **fragments of something longer**, 14 carry flank islands,
and flank backgrounds run 0.40-0.98. Extending every one of them by its decay
distance and re-searching does not rescue them - only `CQ4` (0.0 -> 64.7) and
`AFC3` (52.2 -> 70.6) improve, while `Mar1`, `MARE3` and `NV1` get worse.

### What that means

A cross-species SINE fragment matched at 65 % identity over 50 bp is weak
evidence on its own, and in a repeat-rich scorpion genome most of those matches
land inside longer repeats rather than at independent SINE insertions. The
discriminator says so consistently and for a stated reason, which is the useful
part: **the fragment sweep needs a structural or flank-independence filter
behind it, and the discriminator is that filter.**

Two details worth keeping:

- `sco_OK`, built from the single most abundant query family (357 loci), fires
  **MICROSATELLITE_ELEMENT** - the most common thing the sweep found in this
  genome is simple repeat.
- `sco_SINE_5` from the structural route is rejected at 0.0 with flank
  background **0.885**, so the structural route is not immune either; it is
  caught by the same rule.

Neither route has been looked at by eye yet. The alignments are published and
the seven clean AnnoSINE families are the ones worth his time first.

---

## 2026-09-02 — his de novo scan, found and checked

*"sine de novo scan (not jargon, my script, find it) - find these assemblies and
check run results."*

**The script** is `/data/W/toki/scorpio/denovo_scan/sine_scan.sh`:
`sine_scan.sh <search_db.fa> <queries.fa>`. A chunked ssearch36 fragment search
run in parallel, ≥65 % identity and ≥0.90 query coverage, hits normalised back
from composite chunk names to original scaffold coordinates, then merged with
50 bp of flank. It is the stage the SINEderella README points at when it says
"see the separate SINE-de-novo-genome-scan tool for that earlier stage".

**The run completed 2026-08-14 and is sound:**

| | |
|---|---|
| queries with hits | 792 of 958 |
| total hits | 4,574 |
| merged loci | **3,408** |
| locus length | 141–304 bp, median **148** |
| target | `NW027120657.1` — *Centruroides vittatus* |

Queries are SINEBase family fragments, and every locus header names the query
that found it. Loci per family: OK 357, CQ4 67, BbeA2 55, CQ1 50, Mar1 49, SP1
43, Pca1 39, NV1 38, SacSINE1 38, SINE3 37 and a long tail.

### The assemblies — two, not one per species

Searched `/data/W/toki`, `/data/V/toki` and `/staging`. Exactly two scorpion
genomes exist:

- ***Centruroides vittatus*** — `/data/W/toki/scorpio/sco.bnk`, NCBI ASM3068694
  (= GCF_030686945.1), 770 MB. This is the scan's target.
- ***Buthus sindicus*** — `Genomes/lower/Arthropoda/Scorpions/unp/buthus_sindicus_contigs.fasta`,
  SPAdes contigs (`NODE_1_length_582274_cov_6.5`), **not** from NCBI.

If more species were downloaded they are not in those trees.

### Grouping the loci: by query, not by SubFam chunk

SubFam on the 3,408 loci gives 69 chunk consensuses whose **median pairwise
identity is 0.205** — near random. That is expected and not a failure: SubFam
chunks by input order after a rough reorder, so with loci drawn from 792
different query families each 50-locus chunk spans many families. Chunk-level
grouping only works when the input is already one family, which is how
SINEderella uses it.

For a heterogeneous de novo set the grouping is already in the data: **the query
that found each locus.** 47 families with 20 or more loci, each consensus built
by column majority over a MAFFT alignment of its own loci, all SINE-sized
(117–172 bp).

### My own bug: `cons -plurality 18` on a 20-sequence family

The first attempt used SubFam's own consensus call, `cons -plurality 18`, and
produced 11 consensuses of which **6 got zero blastn hits against the genome the
loci came from** and none scored above 50.

`-plurality 18` asks for 18 of the sequences to agree. In SubFam that runs on
chunks of ~50, where 18 is a little over a third. Applied to a family of 20 it
demands near-unanimity, so most columns come out N — and stripping the N's
splices non-adjacent columns together into a chimera that matches nothing.

**A consensus threshold is a fraction of the sequences present, never a count.**
Rebuilt as a plain column majority over columns where at least half the copies
have a base, the same families give 47 usable consensuses.

---

## 2026-09-02 — re-evaluating the old datasets, and three bugs it flushed out

He asked for the old sets to be re-run with the new criteria: *"the old tal sine,
ere, zebrafush, timema, scorpions (need full re-run - find where they are),
human also need to be re-evaluated using new skills in discrimination."*

The talpid sets (`saq`, `ccr`, `teu`, `dmo`) and hedgehog (`eri`) **are** the
labelled corpus. Human Alu (64 sets) and Timema (55) are in `bench_in`.
Scorpion is *Centruroides vittatus* on KIT at `/data/W/toki/scorpio/sco.bnk`
with three consensuses; zebrafish is GRCz11. Both are downloading on DRAGEN
(GCF_030686945.1, GCF_000002035.6) with AnnoSINE queued behind.

### The new criteria disturb nothing they should not

Re-scored against the previous run:

| dataset | sets | scores moved | verdict flips | new reasons |
|---|---|---|---|---|
| POS, MIXED10/30, NEGJITTER, NEGSPLICE, NEGTRUNC5, NEGRAND, NEGCHIM, MIXSUBFAM, NEGMOSAIC, all SIM*/MOSAIC* | 213 | **0** | **0** | none |
| ERI (hedgehog) | 8 | 3 | 1 | FLANK_ISLANDS x4, SMALL_CORE x2 |
| NEGLINEORF | 4 | 0 | 0 | FLANK_ISLANDS x4, SMALL_CORE x1 |
| NEGSEGDUP | 3 | 2 | 0 | FLANK_ISLANDS x3, SMALL_CORE x1 |
| NEGSAT | 3 | 2 | 0 | FLANK_ISLANDS x3, SMALL_CORE x1 |

Every new reason lands in a flank-related or messy class and nowhere else.

### Human Alu is the cleanest validation the project has

| | sets | accepted | not assessable | rejected |
|---|---|---|---|---|
| **HUM** (Alu subfamilies) | 64 | **55** | 9 | **0** |
| **TIMB** (Timema) | 55 | 48 | 2 | 5 |

**Not one assessable human Alu is rejected.** The nine set aside are alignments
with under 30 copies or no genomic flanks - the tool declining to answer rather
than answering wrongly. The five Timema rejections each name their reasons
(FLANK_ISLANDS + SMALL_CORE on four of them).

### Bug: capping on INSUFFICIENT_COPIES was scoring absence as guilt

`HUM__hum__AluYb9` - which he calls *"no problems good sine"* - came out at
**45.0**, because I had put `INSUFFICIENT_COPIES` in the list of reasons that cap
the score. That is the rule already written down in note 2 of this file and
broken again: **an absent measurement must never be scored as guilt.**

`MICROSATELLITE_ELEMENT` and `SMALL_CORE` are evidence *against*.
`INSUFFICIENT_COPIES` and `NO_FLANKS_PRESENT` are evidence *missing*. The scorer
now carries a third state, `assessable: false`, instead of a low score. Human
rejections went from 1 to 0 and hydra's `hyd_SINE_2` moved from a forced
rejection to "not assessable", which is exactly his reason for it: *"too short,
few copies"*.

### Three bugs in flank decay, all the same mistake

Regenerating the decay profile over the whole corpus - it had only ever been
computed for part of it - dropped the minimum POS score from 90.1 to **0.2**.
Ten sets flipped to rejection including a positive and `NEGCHIM__dmo__d1`, which
he calls a good sine. Three separate defects, each an unmeasured value treated
as measured:

1. **SATELLITE_OR_DUPLICATION was called on distance alone.** No edge guard, so
   a set with flank identity **0.255** - pure background - was called a
   satellite. Real satellites here sit at 0.73-0.95 at the edge; the ten
   casualties were 0.26-0.55. Now requires `edge_max > 0.60`.
2. **An unmeasurable distance defaulted to 400 bp**, the maximum, which then
   collapsed the uniqueness term to zero regardless of the classification. If
   the edge is already at background the element is isolated from its first base
   out: that is a distance of 0, not 400.
3. **The distance was claimed past where the profile was readable.** These
   copies mostly do not carry 400 bp of flank, so the profile goes NaN after two
   or three offsets. `POS__saq__s4_60seqs` is measurable to 50 bp on its right
   flank and was being credited with 400 bp of continuous similarity. It now
   claims only as far as it actually looked.

After all three: **true positives n=132 min 81.5, true negatives n=30 max 45.0** -
still cleanly separated, with every set now carrying real decay data instead of
a neutral fallback.

Two sets remain moved and both are defensible. `NEGTRUNC5__teu__t4` drifts
54.9 -> 49.0, a marginal negative moving the right way. `ERI__eri__e2-3` goes
92.6 -> 45.0 on SMALL_CORE, which converts his *"should be regarded as SINE but
with caution and need manual reinspection"* into the tool declining to assert -
closer to his reading than 92.6 was, though it does cross the line on a set he
leans SINE on. Flagged, not hidden.

### hyd_SINE_5, looked at again

His second reading, with a screenshot of the left end: *"looks like mosaic, i
still cant decide but its not a good sine if at all."*

**Fifth attempt, and the first that measures the right thing.**
`region_partition.py` splits the element into thirds, splits the copies at the
median identity in each third, and asks with an adjusted Rand index whether the
*same* copies group together everywhere. Two subfamilies split the same way
along the element; a mosaic shuffles membership.

On the labelled corpus it does exactly that:

| class | median consistency |
|---|---|
| SIMSUBFAM | **0.560** |
| MIXSUBFAM | 0.064 (max **0.870**) |
| MIXED30 | 0.202 |
| POS | 0.010 |
| NEGMOSAIC / SIMMOSAIC | **0.010 / 0.012** |
| MOSAICKALEID / MOSAICDEL | -0.003 / -0.002 |

`hyd_SINE_5` reads **0.037** with a detected split - the split does not hold
along the element.

**But it cannot be used to reject.** POS is also ~0.010, because a clean family
has no real split to be consistent about, and MIXSUBFAM's median is only 0.064 -
so any threshold catching `hyd_SINE_5` would also catch subfamily mixtures he
considers fine. What the measure gives is *conditional* information: when a
split has been detected, it says whether the split is two families or patchwork.

So it goes into the subfamily note rather than the score. `hyd_SINE_5` now
reads: *"group membership shuffles between the 5' end, the middle and the 3' end
(consistency 0.04, against 0.56 for a real subfamily split): the copies are not
two families so much as patchwork, which is worth looking at by eye."* That is
what he sees, said back to him, without a verdict the evidence cannot carry -
which matches his own position, *"i still cant decide"*.

The score stays 96.6 and it stays the one disagreement.

Hydra now stands at **20 agree, 1 disagree, 1 not assessable** out of 22.

---

## 2026-09-02 — his 22 hydra judgements, and what they changed

He judged every hydra candidate against RepBase, top100 view only, and named two
faults the tool had no measure for. Agreement went from **16 of 22 to 20 of 22**,
then 21 after the consensus repair below.

### RepBase settles the one ambiguity without asking him

He wrote `hyd_SINE_7` twice with two different readings. Rather than ask, all 71
candidate consensuses went against RepBase on KIT. `hyd_SINE_7` is **L2-10_Hmel,
a 2,378 bp LINE**, so the *"combination of microsatellites to me, not sine"*
reading is the one that belongs to it. RepBase confirms his assignment on 18 of
21 and disagrees on `hyd_SINE_4` and `hyd_SINE_10` (Sola2-7_HMa, a DNA
transposon, where he read SINE2-1_HM) — recorded, not resolved.

### Microsatellite content: the criterion he asked for by name

*"looks like combination of microsatellites to me, not sine - not caught by your
filters! - needs new criteria on microsatellite content?"* and, on
`hyd_SINE_17`, *"surrounded by long 2-nt microsatellites"*. Two places, so two
numbers: fraction of the **element** inside a tandem repeat of period 1-6, and
fraction of the **flank**. Median over copies, ungapped, on the raw long flanks.

| candidate | msat element | msat flank | his reading |
|---|---|---|---|
| `hyd_SINE_7` | **0.517** | 0.000 | not a sine, microsatellites |
| `hyd_SINE_17` | **0.229** | **0.253** | surrounded by 2-nt microsatellites |
| every hydra SINE he accepted | <= 0.092 | 0.000 | |

**It is a genuinely new axis.** Over the 673-set labelled corpus, *every* class
has a median element-microsatellite content of **0.000** — POS, NEGSAT,
NEGSEGDUP and NEGLINEORF alike, highest 90th percentile anywhere 0.076. No
existing class exercises it, so no threshold could have been fitted from them;
the bands come from his hydra calls. Thresholds 0.15 (element) and 0.20 (flank),
the latter set above the 0.12–0.17 the synthetic SIM* sets carry in their
generated flanks.

### The background was being measured on the element itself

`hyd_SINE_0`, which he calls a good SINE and RepBase calls SINE2-2B_HM, scored
**0.0 NO_ELEMENT**: *"looks like false alarm on your side"*. Cause: its flank
identity came out at **0.614, above its own element identity of 0.560**, so the
cliff went negative and the support threshold excluded all 100 copies.

The flank was not shared — the consensus covers **79 bp of a 208 bp element**, so
~150 bp of the SINE sat in what the code called flank, and the element was being
measured against itself.

Fix, using data already computed: the decay profile walks outward in 25 bp steps,
so its far end is unrelated DNA wherever the consensus ends. Both the cliff and
the support threshold now use that far-field value. On the labelled corpus this
moved **2 of 273 scores**, both trivially and in the right direction.

### CONSENSUS_UNDEREXTENDED, the mirror that never existed

`CONSENSUS_OVEREXTENDED` had no counterpart. It fires when similarity carries on
past the consensus edge but **ends within ~200 bp** and the far background is
ordinary — a LINE fragment runs on for kilobases, a short consensus does not.

The first version also fired on 4 NEGLINEORF sets and a segmental duplication,
whose consensuses are 260–300 bp and whose decay also reads under 200. What
separates them is the rest of the flank: an under-extended consensus has an
ordinary flank past the missing piece (`hyd_SINE_0`: 6 % in islands), a LINE ORF
has one shared all the way out (37–82 %). With that guard it fires on **0 of
273** labelled sets.

### And then the tool repairs it by itself

`extend_consensus.py` acts on the flag instead of only reporting it:

| stage | length | score | reading |
|---|---|---|---|
| AnnoSINE seed | 110 bp | **0.0** | NO_ELEMENT, consensus too short |
| extended by the decay distance | 310 bp | 0.2 | 88/100 in core at 0.946, now flagged **too long** |
| trimmed to the window the copies support | **211 bp** | **100.0** | clean, 93 of 100 in core |

**211 bp against RepBase's 208 bp for SINE2-2B_HM.** The tool diagnosed the
fault, fixed it, overshot, detected the overshoot with a rule that already
existed, and converged. His three starting situations include "a library needing
adjustment"; this is the tool doing the adjustment.

### Reasons now bind on the score

Three reasons could fire and change nothing: `hyd_SINE_7` was 52 % simple repeat
and scored 77; `hyd_SINE_6` rested on 25 of 100 copies and scored 93;
`hyd_SINE_2` had 8 copies and scored 95 while its own flag text said *"this score
is an impression rather than a measurement"*. He rejected all three.

`MICROSATELLITE_ELEMENT`, `SMALL_CORE` and `INSUFFICIENT_COPIES` now **cap** the
score at 45, just under the acceptance line. Not a subtraction — there is no
principled amount — a statement that on this evidence the answer cannot be yes.
`SMALL_CORE` (core under 35 % of copies) fires on 4 of 273 labelled sets, all in
ERI/NEGSEGDUP/NEGLINEORF, none in POS or any MIXED/NEG class.

Also: flank decay had **never been run on the new species**, so 34 of 52
alignments carried FLANKS_UNMEASURED for no reason. Now measured.

### Where it still disagrees

`hyd_SINE_5` — *"not looking like sine, weak right end, mosaic left end"* —
scores 96.6 with only a subfamily note. Core fraction 0.68, identity 0.641,
nothing else fires.

Two measures were built for it and **neither works**, which is worth recording
rather than hiding:

- **edge sharpness** (`edge_quality.py`: identity over the first/last 25 element
  columns against the middle, and the drop into the flank just outside).
  `hyd_SINE_5` reads ratio5 0.801 / ratio3 0.735, mid-range; `hyd_SINE_6`, also
  rejected, reads ratio3 0.978, near the top. No separation.
- **regional mosaic** (`regional_mosaic.py`: patch2d over the first, middle and
  last third, because he said mosaic *left* end). `hyd_SINE_5` left 0.115 against
  `hyd_SINE_12` left 0.124 and `hyd_SINE_16` left 0.250, both of which he treats
  differently. No separation.

Whole-element patch2d comes closest — `hyd_SINE_5` 0.136 — but `hyd_SINE_13`,
which he accepts, is **0.169**. A threshold that catches `hyd_SINE_5` rejects a
family he accepted, so no threshold was set.

So the property is real and still unmeasured. What is missing is probably not a
sharper edge statistic but a boundary that is defined by the copies rather than
by the consensus: where does each copy actually stop matching, and how much do
those stopping points agree?

**A third attempt at exactly that also failed** (`copy_boundary.py`): slide a
15 bp window along each copy, find the outermost position on each side where it
still matches above chance, then report the spread of those positions. It
returns nonsense - offsets of +304 and -359 on a ~300 bp element, and a spread
of 0.0 for almost every set. The cause is that it walks along **alignment
columns**: MAFFT leaves long gap runs, the window drops below half-valid, the
run of "still element" breaks immediately, and every copy collapses to the
element midpoint.

The fix is to walk each copy's **own ungapped sequence** and map back to columns
only at the end. Not done - recorded so the next attempt starts from the right
place rather than repeating this one.

---

## 2026-09-02 — SubFam, and what a subfamily alignment actually is

He had to correct this three times, so it is worth stating plainly.

**The consensuses fed to a run are families, not subfamilies.** Giving
SINEderella two snail consensuses produces two families with copies assigned to
them. Taking top100/rand100 of those assigned loci — which is what I built first
— is the family-level view that already existed.

**A subfamily alignment is an alignment of SubFam chunk-consensus rows.**
`SubFam` splits the bank into chunks of ~50, makes **one consensus per chunk**,
and aligns those to each other. Collapsing thousands of copies into tens of
consensus rows is what makes subfamily structure visible in one screen.

**And the orchestrator does not run every step.** `SINEderella` runs steps 1–4
and 6; steps **5, 7, 8a, 8b are standalone**. The snail run had no subfamily
alignments because step 5 had never been invoked — I should have checked for a
`step5.*.log` before calling the run complete.

Three more things the docs say that I had not read:

- `ALIGNMENT_DEPLOYMENT.md`: **"CRITICAL: Use extract_alignments.sh, NOT
  step5_align_subfamilies.sh"** — step 5's output has no flanks, and flanks are
  mandatory for boundary inspection. `extract_alignments.sh` writes `.core.fa`
  (50 random copies), `.best50.fa` (top 50 by bitscore) and `.subfam.fa` (the
  SubFam rows, only at >= 400 copies), all with 30 bp left and 70 bp right
  flanks.
- MANUAL §6.1.1: **`input.clw` is not a trustworthy alignment** — independently
  derived batch consensuses that ended up discordant. Degap and realign before
  using it. That is why my hand-run of step 5 produced empty files; the failure
  was mine, not the pipeline's.
- MANUAL §6.1 is the **subfamily discovery** workflow for when the consensus set
  is too coarse: realign the chunk consensuses, group them by eye, build a
  consensus per group with `sine_consensus.sh`, re-run.

Done for the snail: all three alignment tiers for both families, and the 178
chunk consensuses degapped and properly realigned
(`POM__subfam_chunks.realigned.fa`) so the grouping can be done by eye. Written
up as a `sinederella` skill.

---

## 2026-09-02 — the flank islands, tested against his own judgements

One candidate is not a finding, so the island scan was run over the whole
labelled corpus on KIT (`islands_corpus.py`, 673 of 674 sets), measuring the
**fraction** of the flank that sits inside an island rather than the raw count -
a set with a 2,400 bp flank has more room for islands than one with 400.

### It separates exactly the classes whose problem IS the flank

| class | n | median island fraction | median flank_bg |
|---|---|---|---|
| NEGSEGDUP | 3 | **0.871** | 0.793 |
| NEGLINEORF | 4 | **0.819** | 0.493 |
| NEGSAT | 3 | **0.523** | 0.877 |
| ERI (hedgehog, real but messy) | 8 | 0.115 | 0.283 |
| NEGLINE / SIMNEST / NEGLINE2 | 12 | 0.055-0.061 | ~0.32 |
| **POS** | 28 | **0.025** | 0.309 |
| NEGRAND | 20 | 0.000 | 0.264 |

Segmental duplications, LINE ORFs and satellites: the three classes defined by
copies not being in independent places. Everything else sits at 0.02-0.06.

### It never misses what the average catches, and catches three more

All **11** sets with `flank_bg >= 0.40` also have island fraction >= 0.15, so
nothing is lost by using it. Three sets go the other way - high islands with an
ordinary flank average, which is the `aca_SINE_0` pattern:

| set | island frac | flank_bg | score |
|---|---|---|---|
| `NEGSEGDUP__eri__r00` | 0.368 | **0.326** | 47.0 |
| `ERI__eri__e2-4` | 0.195 | **0.281** | 97.4 |
| `ERI__eri__e2-3` | 0.194 | **0.371** | 92.6 |

`NEGSEGDUP__eri__r00` is a **labelled segmental duplication whose flank average
is 0.326** - indistinguishable from a clean family. Only the islands see it.

### And all three are sets he had already flagged for flank checks

This was not fitted; it fell out. His verbatim calls on exactly those three:

- `ERI__eri__e2-4` - *"nested copies, i agree, but there are enough uniq left
  flanks to lean towards SINE - requires post-processing with **proving uniqness
  of at least some flanks on whole-genome level**"*
- `ERI__eri__e2-3` - *"can be SINE, but very turbulent at manu places, should be
  regarded as SINe but **with caution and need manual reinspection**"*
- `NEGSEGDUP__eri__r00` - *"**unclear situation with left end**, could be just
  not extended enough, more like grey zone or technically badly prepared"*

Across all 55 judged sets that have an island measurement, the fraction falls
into three bands with wide gaps between them:

| his call | n | island fraction |
|---|---|---|
| plain `SINE` | 31 | **<= 0.067** |
| the three flank-caution calls | 3 | **0.19 - 0.37** |
| `NOT_SINE` / `UNUSABLE` / `MOSAIC` | 7 | **>= 0.50** |

No plain SINE reaches 0.07. No caution call is below 0.19. No rejection is
below 0.50.

### So it goes in as a reported reason, not a score penalty

He still leans SINE on two of the three. Subtracting for islands would turn
`e2-4` (97.4, which he accepts with a caveat) into a rejection, which would be
wrong. `verdict.py` now emits **`FLANK_ISLANDS`** at fraction >= 0.10, with
different wording above 0.45, and says explicitly when the flank average is
ordinary - because that is the case a person cannot get any other way.

Thresholds 0.10 and 0.45 sit in the empty gaps, not on top of any judged set.

### Where the new candidates land on his ladder

| candidate | island fraction | flank_bg | band |
|---|---|---|---|
| `hyd_SINE_9` | 0.331 | 0.914 | caution, but the average already rejects it |
| `aca_SINE_1` | 0.217 | 0.450 | caution - a **fourth** independent fault |
| **`aca_SINE_0`** | **0.171** | **0.271** | **2.5x the highest plain SINE he judged** |
| `hyd_SINE_0`, `hyd_SINE_17` | 0.088-0.091 | 0.61 / 0.54 | slightly raised |
| the other 21 | <= 0.029 | ~0.31 | clean |

`aca_SINE_0` sits between the plain-SINE ceiling (0.067) and the caution floor
(0.19), nearer the caution end - which is what his eye reported.

**And the signal is in the top hits only.** `aca_SINE_0__top100` is 0.171;
`aca_SINE_0__rand100` is under 0.03. `hyd_SINE_9` is 0.331 on both views. So
aca_SINE_0 is not a family sitting in a shared context - it is a family where
the **most similar copies** share one, which is what a subset inside a larger
duplication looks like. That is a concrete next check: pull those flanks and
search them against the starfish genome.

---

## 2026-09-02 — all three genomes rebuilt, and where the islands actually live

### The five hydra candidates that were never scored

Hydra's alignments finished at 13:34; the scoring ran at 13:27. So
`hyd_SINE_17` through `hyd_SINE_21` had alignments and no verdicts, and the
results page showed 17 of 22 candidates without saying so. Everything is now
rebuilt in one pass (`newsp_all.py` on DRAGEN): raw alignment, flanks
justified, flank width trimmed, scored, islands measured. **52 alignments over
26 candidates**, all published.

Hydra: 19 of 22 candidates score 50 or above on their best view. The three that
do not are `hyd_SINE_0`, `hyd_SINE_9` and `hyd_SINE_16`, all with NO_ELEMENT;
the first two also with SHARED_FLANKS.

### Trimming the flanks destroys the island signal

This one matters for how the tool is put together. Measured on the **trimmed**
alignment, `aca_SINE_0` has **6** island columns. Measured on the **full 400 bp
flank**, it has **426**.

Trimming cuts the flank panel to the width the copies actually fill, which is
right for looking at a boundary and is what fixed his complaint about walls of
dashes. But it removes exactly the far-out columns he said the islands sit in
(*"these islands can be far off on the flanks (both directions)"*).

**So the island scan must never run on the display alignment.** Same rule as
flank decay: short flanks to look at, long flanks to measure. Recorded in
`ALGORITHM_NOTES.md` as part of note 10.

### Rescanned uniformly, only ONE candidate is invisible to the score

All 50 measurable alignments scanned with one null (per-column pairwise
identity against that set's own flank composition, SD from the number of pairs
actually present, columns with under 8 copies not scored):

| alignment | patches | max z | island cols | flank bg | score | already caught? |
|---|---|---|---|---|---|---|
| `hyd_SINE_9` top100 | 50 | 112.2 | **739** | 0.914 | 0.0 | yes, SHARED_FLANKS |
| `hyd_SINE_9` rand100 | 46 | 110.8 | 722 | 0.800 | 0.0 | yes |
| **`aca_SINE_0` top100** | 49 | 69.1 | **426** | **0.271** | **100.0** | **no** |
| `aca_SINE_1` all | 38 | 74.9 | 424 | 0.450 | 26.3 | partly, bg raised |
| `hyd_SINE_17` top100 | 28 | 88.1 | 303 | 0.543 | 100.0 | yes |
| `hyd_SINE_0` top100 | 10 | 109.9 | 139 | 0.614 | 0.0 | yes |
| everything else | <= 10 | | <= 75 | ~0.30 | | too little to matter |

Taking "invisible" to mean over 150 island columns with a flank average under
0.40, **`aca_SINE_0` is the only one in the entire corpus.** Every other set
with a big island count has a raised flank average too, so the existing
SHARED_FLANKS rule already fires on it.

That is the whole case for the measure: it is not that averaging misses
similarity in general, it is that averaging misses *localised* similarity, and
exactly one candidate here has localised similarity without the global kind.

**Earlier numbers superseded.** The 839 / 754 / 713 counts I reported before
came from a scan with a different null; these come from one method applied to
all 50 sets and are the ones to use. The ordering and the conclusion are
unchanged.

### Correction: hyd_SINE_9 is not a second blind spot

I wrote that `hyd_SINE_9` (713 island columns) was a second case the score
could not see, and predicted it would look like `aca_SINE_0` by eye. That was
wrong. Its flank background is **0.914** — the score rejects it outright with
SHARED_FLANKS and NO_ELEMENT. It should look *different* under the eye: similar
all the way along the flank rather than in patches. `aca_SINE_0` is the one to
check.

### The results page

`site/newspecies.html` rebuilt from the data files with no hand-typed numbers:
all three species, all 26 candidates, every view on its own row, links to the
**justified and trimmed** alignments, a flank-island column with a magnitude
bar, the full SINEderella table for the snail, and hover text on every flag and
every column heading. Five tables, all fitting without horizontal scroll.

---

## 2026-09-02 — SINEderella on the snail candidates, and aca_SINE_1's full reject case

### Full SINEderella run on the two refined snail consensuses

Consensuses refined first from their 100 aligned copies rather than the
single-locus AnnoSINE seeds (174 bp and 261 bp). Run
`/staging/tmp/newsp/pom/run_20260901_131644`, completed in about 4 minutes.

| | firm | soft | **total assigned** | **leak** | conflicts | sim_mean |
|---|---|---|---|---|---|---|
| `pomSINE0` | 4,712 | 1,304 | **6,016** | **0.00 %** | 3 | 0.737 |
| `pomSINE1` | 1,249 | 1,599 | **2,848** | **0.00 %** | 0 | 0.671 |

**AnnoSINE had estimated 537 and 339 copies — SINEderella finds 11x and 8x
more.** A discovery tool's copy count is a floor, not an estimate.

**Zero leak on both.** In Timema the candidates he judged real ran 0.00–0.18 %
leak and the noisy ones 65–98 %, so both snail families behave like his real
ones. `pomSINE1` is again the weaker: firm assignment 14 % against 53 %, and
sim_mean 0.671 against 0.737.

### aca_SINE_1: all three of his additional reject reasons confirmed

He said the rejection is *"not only mosaic, but also weak identity, too long
identity with bad flanks"*. Measured:

| | cons bp | id_all | flank bg | cliff | span/cons | patch2d |
|---|---|---|---|---|---|---|
| `aca_SINE_0` (he: legit) | 281 | 0.723 | 0.271 | 0.452 | 2.76 | 0.097 |
| **`aca_SINE_1`** | 243 | **0.552** | **0.450** | **0.102** | 2.52 | **0.143** |
| `pom_SINE_0` | 175 | 0.834 | 0.307 | 0.527 | **1.07** | 0.061 |
| `pom_SINE_1` | 262 | 0.730 | 0.279 | 0.451 | 2.26 | 0.150 |

- **weak identity** — 0.552 against 0.72–0.83. Confirmed.
- **bad flanks** — flank background **0.450**, nearly double the 0.27–0.31 of the
  others. Its copies share flanking sequence and are therefore not independent
  insertions. Confirmed, and stronger than he stated.
- **almost no boundary** — cliff 0.102 against 0.45–0.53. Confirmed.
- **mosaic** — patch2d 0.143. Confirmed.

**The tool scored it 26.3 with flags "clean".** Right number, no reason. Four
independent faults and it reported none of them. Recorded as the central defect
in `ALGORITHM_NOTES.md` §1.

Note `span/cons` — how far the consensus is smeared across the alignment.
`pom_SINE_0`, the strongest candidate, is 1.07; everything else is 2.26–2.76.
That may be worth a measurement of its own.

---

## 2026-09-02 — R12. The mosaic he means, caught in the wild — and a measurement for it

He said months of examples earlier: *"What i mean with mosaic looks truly
different and maybe we can get back to it when we'll catch it in the wild."*
None of the six synthetic MOSAIC* constructions was it — he called all six SINEs.

**`aca_SINE_1` (starfish) is it.** His words: *"sine_1 is bad and it includes
mosaicism - inconsistent mixture of conserved and discordant columns/spots"*.
The tool had already rejected it at 26.3, but for the wrong reason (weak overall
identity), not for being a mosaic.

### What does NOT measure it

- **Column jumpiness** — `aca_SINE_1` 0.261 vs a legitimate `pom_SINE_0` at
  0.239. No separation.
- **Block structure** (runs test against chance): `aca_SINE_1` blockiness 1.17,
  while the legitimate `aca_SINE_0` is the second-blockiest at 1.30. Backwards.
- **Mean identity** — just says "less conserved", not "mosaic".

### What does: per-copy variation across regions (`patch2d`)

For each copy, its concordance with the consensus in each 20 bp window; then how
much that varies **within** a copy. A clean family: every copy is uniformly
concordant. A mosaic: a copy is concordant in some regions and discordant in
others, and which regions differ between copies.

| candidate | his judgement | patch2d |
|---|---|---|
| `pom_SINE_0` | strongest — 86 % TSD, perfect B-box | **0.061** |
| `aca_SINE_0` | "sine_0 looks legit" | **0.097** |
| `aca_SINE_1` | **"bad ... includes mosaicism"** | **0.143** |
| `pom_SINE_1` | "possible sine too" — only 35 % TSD | **0.150** |

**It orders all four exactly as he does.** Four points is not a threshold, but it
is the first measurement that responds to the thing he has been describing since
the mosaic discussion began, and it is row-wise-within-column — neither purely
per-column nor purely per-copy, which is why every earlier attempt missed it.

Note `pom_SINE_1` scores highest of all, above the acknowledged mosaic. That is
consistent with its other weaknesses (A-box 3 mismatches, TSD 35 %, no tail) and
suggests it deserves the same scrutiny.

**Needs more of his judgements to become a threshold.** The obvious source is
the 22 Hydra candidates, which are unjudged and from a fourth phylum.

---

## 2026-09-02 — pom_SINE_1 structure, and both candidates vs RepBase

Alignments rebuilt at 100 bp flanks, justified and degapped, published as
`NEW__pom100_*__*.degap.aln.fa`. Flank gaps:

| | raw | justified | degapped |
|---|---|---|---|
| SINE_0 top100 | 0.69 / 0.68 | 0.07 / 0.04 | **0.07 / 0.03** |
| SINE_1 top100 | 0.74 / 0.79 | 0.19 / 0.45 | **0.19 / 0.22** |
| SINE_1 rand100 | 0.78 / 0.87 | 0.38 / 0.70 | **0.23 / 0.30** |

### The two candidates side by side

| | pom_SINE_0 | pom_SINE_1 |
|---|---|---|
| length | 175 bp | 262 bp |
| AnnoSINE source | tRNA | 7SL RNA |
| genomic copies | 537 | 339 |
| **A-box** | pos 2-12, **2 mismatches** | pos 0-10, **3 mismatches** |
| **B-box** | pos 54-62 `GTTCAACTC`, **0 mm** | pos 50-58 `GTTCAAGAC`, **0 mm** |
| **TSD** | **86 % of copies**, median 13 bp | **35 % of copies**, median 8 bp |
| Pol III terminator | 24 % strong, 8 % moderate | 16 % strong, 9 % moderate |
| 3' tail | poly-C, 8 units | none detected |

**Both have a perfect B-box** — an intact RNA-pol-III internal promoter, which is
the hardest part to acquire by chance.

**SINE_1 is the weaker of the two**: its A-box has 3 mismatches rather than 2,
only 35 % of copies carry a TSD against 86 %, and no 3' tail repeat is detected.
Real but less convincing, and worth his eye on the degapped alignment.

### Neither has a relative or a LINE partner

Searched at `word_size 7`, no filtering, against **49,007 RepBase entries**,
of which **1,784 L1** and **7,117 other LINE families** (RTE, CR1, Jockey, L2,
Rex, R1/R2, Nimb, Penelope):

| query | whole RepBase | best LINE hit | where on the LINE |
|---|---|---|---|
| pom_SINE_0 | 3 hits, e 0.25-0.89 | 28 bp, e 0.14 | **1,246 bp** from the 3' end |
| pom_SINE_1 | 1 hit, 19 bp, e 0.39 | 28 bp, e 0.22 | **3,000-4,600 bp** from the 3' end |

**The "similarity to L1 tails" does not survive a targeted search.** Every hit is
14-28 bp with e >= 0.14, and none is near a LINE 3' terminus. A genuine SINE-LINE
partnership puts the SINE's tail on the LINE's last ~50 bp; these sit thousands
of bases inside. It is background.

Combined with the SINEBase survey above — 23-55 % of that database is isolated —
being unrelated to known elements is a normal property of a real SINE, not
evidence against one.

**Open question, unresolved by sequence alone:** if neither has a LINE partner in
RepBase, what mobilises them? Molluscan LINEs are poorly represented there, so
the partner may simply be undescribed — findable by searching the snail genome
itself for a LINE whose 3' end matches these SINEs' tails.

---

## 2026-09-02 — pom_SINE_0 structure, and the SINEBase isolation survey

### pom_SINE_0 has a complete SINE structure

He asked whether it has TSDs or boxes. It has all of them:

| feature | detail |
|---|---|
| **A-box** | positions 2-12, `TGGCGCAACAC`, 2 mismatches |
| **B-box** | positions 54-62, `GTTCAACTC`, **0 mismatches** |
| **TSD** | present in **86 % of copies**, median **13 bp** |
| also detected | tRNA region, conserved core, terminator, tail repeat |

A perfect B-box and TSDs on 86 % of copies means an intact RNA-pol-III internal
promoter and direct repeats flanking each insertion. That is strong internal
evidence, and it is consistent with his reading: a real SINE with no relationship
to known families. 537 genomic copies.

Rebuilt with 100 bp flanks (was 400) on his instruction — published as
`NEW__pom100_SINE_0__*`.

### SINEBase all-vs-all: how much of it is isolated?

The only local SINEBase was 958 x 50 bp fragments already reduced to 85 %
non-redundant — useless for this question. Downloaded the real bank from
sines.eimb.ru: **231 entries, 230 unique ids** (one duplicate header).

blastn all-vs-all, `word_size 7`, self-hits removed:

| significance cut | isolated | share |
|---|---|---|
| e <= 1e-20 | **127** | 55 % |
| e <= 1e-10 | 94 | 41 % |
| **e <= 1e-5** | **53** | **23 %** |
| e <= 1e-3 | 16 | 7 % |
| e <= 0.01 | 5 | 2 % |

**The threshold decides the answer**, which is itself the finding. At a strict
cut more than half of SINEBase has no relative in it; at a loose cut almost
everything is connected to something.

Connectivity of the rest: **median 2 neighbours, maximum 15**. So SINEBase is
not one connected family — it is mostly small clusters plus a large isolated
tail.

The 53 isolated at e<=1e-5 are in `sinebase_isolated.txt`: `ACar-1`, `AFC-3`,
`Asu1`, `Au`, `Bm1`, `BmSE`, `BraS-I`, `CQ-2`, `Cre-3`, `Cry`, `CucuS-II`,
`DR-2`, `ERI-2`, `EuphS-I`, `FabaS-I/III/V/VI/VII/IX`, `Feilai`, `IscNinDC`,
`Lj-2`, `LF`, `Lm1`, `Mad-1`, `Mad-2`, `NV-2`, `NymS-I`, `OK`, `p-SINE1`,
`PinS-I`, `PoaS-I`, `PoaS-II`, `Ruka`, `SaliS-III`, `SB1`, `SB3`, `SB4`,
`SB12`, and 15 more.

**Bearing on pom_SINE_0:** an unrelated SINE is not unusual. Somewhere between a
quarter and a half of the established database has no significant similarity to
anything else in it, so "no similarity to known SINEs" is a normal property of a
real SINE family, not a reason for doubt.

---

## 2026-09-02 — First prospective results: candidates with no answer key

`candidate_to_aln.py` is new and is the step that never existed: candidate
consensus + genome -> ~100-locus alignments with flanks, in the three views he
uses (top hits / random sample / everything). Plan section 6.1, working.

### AnnoSINE_v2 on three genomes with no SINE library

| genome | phylum | candidates found |
|---|---|---|
| *Pomacea canaliculata* (snail) | Mollusca | **2** |
| *Acanthaster planci* (starfish) | Echinodermata | **2** |
| *Hydra vulgaris* | Cnidaria | **22** (20 tRNA-derived, 1 5S, 1 unknown) |

Against Timema's 55. He read the low counts correctly: snail and starfish are
genuinely SINE-poor, not a failed run — all seven AnnoSINE steps completed and
the candidates found are substantial, not junk.

### Verdicts on the four snail/starfish candidates

| candidate | derived from | genomic copies | top100 | rand100 |
|---|---|---|---|---|
| `pom_SINE_0` | tRNA | **537** | **100.0** clean | 93.7 TRUNCATED_COPIES |
| `pom_SINE_1` | 7SL RNA | **339** | **100.0** clean | 86.2 TRUNCATED_COPIES |
| `aca_SINE_0` | tRNA | 176 | **100.0** | 95.9 SUBFAMILY_NOTE |
| `aca_SINE_1` | 5S rRNA | 11 | **26.3** | 26.3 |

Three accepted, one rejected. **Nothing can check this except his eye** — there
is no curated library, no prior annotation, no RepeatMasker track for these
genomes. Published at `newspecies.html`.

Note the pattern in both snail candidates: **top100 clean at 100.0, rand100
flagged TRUNCATED_COPIES**. The best-scoring hits are full-length; a random
sample includes fragments. That is expected for a real family and is itself
weak evidence these are genuine.

### Two more genomes running

Chosen on his instruction — *"find other genome presumably containing sines but
make sure its not properly annotated/studied"*:

- *Acipenser ruthenus* (sterlet sturgeon) 1.9 Gb
- *Amblyraja radiata* (thorny skate) 2.6 Gb

Both are groups where SINEs are documented in the literature while these
assemblies carry no repeat library.

---

## 2026-09-01 — R11. Flanks in the aligner: a hypothesis of mine, DISPROVED

He asked: *"why you insist on aligned flanks in your examples?"*

The two corpora were built differently. `kit/extend.py` (400 bp) aligns the
element with the consensus and attaches raw flanks afterwards — its own comment
says *"no flank in the aligner at all"*. But `aln_v2`, which every review
alignment came from, put element and flanks through MAFFT together and de-gapped
the flanks afterwards.

**My hypothesis:** letting unalignable flanks into the aligner distorts the
ELEMENT alignment, because MAFFT optimises one score across the whole row.

**Tested and disproved.** Re-aligning the element alone and attaching raw flanks
changes the element alignment barely at all:

| alignment | element gaps, flanks in aligner | element only |
|---|---|---|
| `NEGCHIM__ccr__g1_180seqs` | 0.752 | 0.789 |
| `NEGCHIM__ccr__g3_71seqs` | 0.765 | 0.769 |
| `MIXSUBFAM__teu__t1_45seqs_t5_31seqs` | 0.451 | 0.434 |
| `NEGTRUNC5__saq__s2_38seqs` | 0.663 | 0.654 |
| `NEGCHIM__saq__s8_225seqs` | 0.738 | 0.747 |

Two go slightly worse, three slightly better - noise. MAFFT was already
effectively ignoring the flanks. **So this is not the cause of what he is
seeing**, and the question remains open: after justification each copy's flank
is raw sequence butted against the element with padding on the outer side only,
because a FASTA alignment must be rectangular. Whether that padding is what he
reads as "aligned" is unresolved and was asked back to him rather than guessed
at a third time.

`realign.py` is kept on DRAGEN in case element-only rebuilding is wanted for
another reason, but it does not fix this.

---

## 2026-09-01 — R10. The mixture test rebuilt: 5/9 → 9/9, but it now fires 3x more often

### The old measure was invalid, not miscalibrated

Group **length** difference, tested against nine of his judgements:

| | rel_len |
|---|---|
| his mixtures | 0.122, 0.133, 0.138, 0.309 |
| his non-mixtures | 0.089, 0.119, **0.228** |

`s8_225seqs` is "not mixture" at 0.228 — above three of the four real mixtures.
The ranges overlap completely, so **no threshold exists**. Re-fitting was the
wrong response; the measure had to be replaced.

### What replaced it

Group **identity** difference — which is what he actually describes: *"top proper
about half longer similarity sequences ... bottom more discordant shorter ones"*.

| | d_ident |
|---|---|
| his mixtures | 0.128 – 0.215 |
| his non-mixtures | 0.012 – 0.124 |

No overlap. Threshold 0.126.

### The grouping mattered as much as the measure

A first attempt applied that threshold to `subfamily_split`'s grouping and got
7/9. Under that grouping CAS and AluYh9 — both of which he calls mixtures —
score **0.013 and 0.012**, the lowest of all nine. `subfamily_split` is driven by
sequence structure and finds the wrong halves.

Splitting instead on each copy's own **identity, length and coverage** (the SVD
split in `row_groups`) recovers 0.209 and 0.128. **9 of 9.**

### Regression, and one number that needs his eye

| | n | mean | extreme | crossing 50 |
|---|---|---|---|---|
| true negatives | 30 | 1.59 | max 47.0 | 0 |
| true positives | 432 | 99.96 | min 90.1 | 0 |
| curated `tim/` | 42 | 98.29 | min 71.3 | 0 |

Unchanged. **But the flag now fires on 120 of 834 sets (14%), up from 39 (4.7%).**

| class | fires | of | |
|---|---|---|---|
| MIXED10 | 28 | 28 | 100% |
| MIXED30 | 28 | 28 | 100% |
| NEGTRUNC5 | 28 | 28 | 100% |
| HUM | 14 | 64 | 22% |
| TIM | 13 | 42 | 31% |

MIXED10/MIXED30 are *deliberately contaminated*, so firing there is arguably
right — contamination is two things in one set. NEGTRUNC5 at 100% is more
doubtful: he called those *"more like a SINE"*, not mixtures.

**This flag withholds the verdict**, so 120 deferred sets is a large behavioural
change resting on nine judgements. Needs his eye on a NEGTRUNC5 and a MIXED30
before it can be trusted at this rate.

---

## 2026-09-01 — Review batch 2 opened; two failures he caught

### R7. `CONSENSUS_OVEREXTENDED` has a false negative

He judged `NEGCHIM__ccr__g1_180seqs`: *"overextended consensus and no degapped
right flank"*. The tool flags `FLANKS_UNMEASURED, RECOVERABLE_CORE` — **the
overextension flag does not fire**, although it fires correctly on `g3` and `g5`
which he confirmed.

So the flag is right when it fires (4 of 4 in batch 1) but misses cases. Its
condition — a core window beating the whole span by ≥0.15 with the remainder at
background — is too strict. Needs re-fitting against `g1` as a positive.

### R8. I asked about an alignment whose answer he had already given

> "it feels like you havent tried to learn anything from my previous answers
> because it perfectly matches what i have already answered"

He had judged **eight** NEGCHIM sets, five of them *"consensus too long but good
sine"*. I put a ninth in the next batch. He answered it identically.

Fixed by `already_answered.py`: before asking, check whether he has answered
enough alignments of the same construction, consistently, and if so predict
rather than ask.

Grouping by exact label was too literal — his nine NEGCHIM answers use six
different labels and the top one is only 44%. Grouped by **substance** (is it a
SINE at all, ignoring which extra property is named) it is **89% "a SINE"**.
That is the level at which his answers agree.

Constructions now predictable without asking:

| construction | judged | agreement | answer |
|---|---|---|---|
| NEGSPLICE | 20 | 100% | a SINE |
| NEGSAT | 3 | 100% | not a SINE |
| NEGCHIM | 9 | 89% | a SINE |
| ERI | 4 | 75% | a SINE |
| HUM | 26 | 62% | a SINE |

### R9. The flank trim needed two more fixes

His *"no degapped right flank"* on `g1_180seqs` was still true after the first
fix. Two causes:

1. A 75th-percentile cap is not enough when flank lengths are skewed — right
   flank median 20 bases but p75 58, leaving the panel 53% gaps. Now the cap is
   chosen as the widest width whose panel is at most 25% gaps.
2. A constant floor of 25 columns cannot help when copies carry ~10 bases. The
   floor now follows the copies' own median.

| alignment | right-flank gaps: raw → first fix → now |
|---|---|
| `NEGCHIM__ccr__g1_180seqs` | 0.78 → 0.53 → **0.36** |
| `NEGCHIM__ccr__g3_71seqs` | 0.85 → 0.55 → **0.42** |
| `MIXED10__ccr__g4_77seqs` | 0.29 → 0.25 → **0.25** (left 0.80 → **0.23**) |

The residue is intrinsic: those copies genuinely carry only 10–20 bases of right
flank, so any panel showing them has gaps for the shorter ones.

---

## 2026-09-01 — Review batch 1: 16 judgements, six results

Recorded verbatim in `calls.tsv` (now 87 judgements).

### R1. The mixture test is wrong 3 times out of 4 — must be re-fitted

| alignment | what he said |
|---|---|
| `NEGCHIM__saq__s5_5seqs` | "not mixture" |
| `NEGCHIM__saq__s8_225seqs` | "not mixture" |
| `NEGCHIM__ccr__g5_7seqs` | "i dont see mixture" |
| `MIXSUBFAM__teu__t1_45seqs_t5_31seqs` | "mixture of 2 subfamilies" — the only true one |

`HETEROGENEOUS_SELECTION` fires on all four. **75 % false positive rate.**

Its threshold (relative group length difference ≥ 0.12) was fitted to three
examples and was already flagged as overfit in `HANDOFF.md` §53a. This confirms
it. And the one true positive is a *subfamily* mixture — which by his own rule
does not affect the sine/not-sine verdict anyway.

**Action:** re-fit against these four; reconsider whether a subfamily mixture
should trigger it at all.

### R2. The over-long-consensus test is right 4 of 4 — but it is a NOTE, not a fault

| alignment | what he said |
|---|---|
| `NEGCHIM__ccr__g3_71seqs` | "yes consensus is too long, but its a good sine" |
| `NEGCHIM__ccr__g6_58seqs` | same |
| `NEGCHIM__dmo__d1_16seqs` | same |
| `NEGCHIM__ccr__g5_7seqs` | "consensus too long in 3' end" |

`CONSENSUS_OVEREXTENDED` is correct every time, **but he calls all of them good
SINEs**. It must never lower the verdict — it is an instruction to trim, not
evidence against the family. He also localises it further than the tool does
("3' end"): the flag should say which end.

### R3. The copy-count threshold is bounded, and 30 is too high

| alignment | n | what he said |
|---|---|---|
| `HUM__hum__AluYb9` | 21 | "no problems good sine" |
| `TIMB__timb__SINE_46` | 11 | "not enough sequences but left flank is probably bad or incomplete" |

`INSUFFICIENT_COPIES` fires below 30, so it wrongly fires on `AluYb9`.
**The line sits between 11 and 21.** Needs a judgement near n≈15 to place it.

### R4. Display bug — flank width set by one outlier copy *(FIXED)*

> "flanks are badly degapped and it hinders my estimates"

The flanks **were** justified — every copy sits hard against the element,
gap-to-element 0 for all. The fault was the column **width**, set by the single
longest copy while the median copy has far fewer bases.

`NEGCHIM__ccr__g3_71seqs`: right flank 125 columns, median copy has 10 bases in
it, panel 85 % gaps. One outlier padding seventy copies with whitespace.

Fixed by `trim_flanks.py` — cap each flank at the 75th percentile of the copies'
own flank lengths:

| alignment | flank gaps before | after |
|---|---|---|
| `NEGCHIM__ccr__g3_71seqs` | 0.56 | **0.16** |
| `NEGSAT__eri__r02` | 0.68 | **0.26** |
| `MIXSUBFAM__teu__t1_45seqs_t5_31seqs` | 0.67 | **0.07** |
| `NEGCHIM__saq__s8_225seqs` | 0.45 | **0.23** |

All review alignments republished as `.trim.aln.fa`.

### R5. Some published alignments are too poor to judge — confirmed

> "you provided it unaligned; too long and no flanks"

Technically aligned (all rows equal length) but **59–70 % gaps, smeared to 4–6×
the element length**, against 33 % and 2.3× for a good one. Any alignment past
roughly 3× its consensus length should not be shown for judgement without
saying so.

### R6. Batch size

> "actually i asked by 6 not to overload me"

**Six alignments per review batch.** Batch 1 had 28.

---

## 2026-09-01 — Tests 1–3

### T1. Three existing measurements carry what the score discards

- **`edge_drop`** (how far copy-to-copy similarity falls at the element edge):
  clean 0.52, grey 0.15, badly-presented 0.18, wrong-consensus 0.13 — all of
  which score 89–100 and are identical by score.
- **`pair_bg`** explains `ERI e1-4`: flanks 0.57 similar against a normal 0.25.
  The copies share flanking sequence — a different statement from "no element".
- **`cov_min` = 0.000** marks every badly-presented alignment.

### T2. The composition-aware 3′ edge rule WORKS — the plan's go/no-go, passed

Identity alone fails: it collapses in the poly-A tail because A-run length
varies between copies (AluY 0.93 in the body → 0.27 in the last 90 bp, while
A+T rises 0.26 → 0.63). An identity-only rule amputates the tail from every SINE.

Watching composition as well as identity, the tail correctly kept:

| what he said | tail rescued |
|---|---|
| SINE | **53.9 bp** |
| SINE_MINOR_NOTE | 42.5 bp |
| MIXTURE_SPLIT_FIRST | 56.0 bp |
| NOT_SINE | **10.6 bp** |
| ANCIENT_BADLY_PRESENTED | 9.3 bp |
| SINE_CONSENSUS_WRONG | **0.0 bp** |

Real SINEs carry ~50 bp of A/T-rich tail; things he called not-SINEs carry
almost none. So it is both a boundary rule and a discriminating measurement.

### T3. `rowsplit_dident` earns its place; local instability does not

Of 10 comparable pairs, all separate. **`rowsplit_dident`** — how much two
groups of *copies* differ in identity — is the best separator for both mixture
comparisons, a property nothing measured before (everything else is per-column;
he was describing groups of rows).

**`dip_depth`** (local instability) separates nothing: 0.55 for clean SINEs vs
0.61 for the "island of instability" case. Dropped.

### T3-bug. A zero-variance artifact, caught before reporting it as success

The first Test 3 run reported separations of ~1e11 and "0 of 120 pairs not told
apart". Nine of his sixteen properties have a single example, whose variance is
zero, so the pooled spread collapses to the epsilon. Now requires 3 examples per
side — leaving **only 5 of 16 properties comparable at all**. That is why Test 4
exists.

---

## Standing conclusions

- **He names a property and attaches an action; the tool emits a number.** Four
  of his properties all score exactly 100.0. No threshold tuning fixes this.
- **Human Dfam is not ground truth.** His judgements and the curated `tim/` set
  are.
- **Two thresholds are fitted to almost nothing** and must be re-fitted: the
  mixture test (now measured at 75 % false positive) and the over-long-consensus
  tail bound.
