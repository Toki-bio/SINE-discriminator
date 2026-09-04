# The two schooling problems: what exists, how ready it is, and whether it is doable

Written 2026-09-04. Every number here is measured. Where a number was obtained by
tuning against the answer, it says so.

This is not an inventory. Four experiments were run for it, three of which changed
my view of whether the work is doable at all.

---

## Executive summary

| | Track 1 — SINE or not | Track 2 — subfamily separation |
|---|---|---|
| Expert labels | **113 verdicts** with verbatim reasoning | 6 blocks on saq, 2 on teu |
| Best measured performance | 0.98 on the binary call, on **6 negatives** | **1 of 6** finding boundaries unaided |
| Biggest surprise | **the negative corpus is 81 % mislabelled** | **the criterion has a formulation flaw, and it is fixable** |
| Readiness today | triage filter, human-in-the-loop | assistant, not classifier |
| Doable? | binary yes; the 28-way verdict, partly | **yes — better than I thought a day ago** |
| Real blocker | no genuine negatives except from the de novo scan | thresholds fitted on n=6 |

**The headline for track 1** is that the thing everyone assumed was the asset — 42
synthetically-degraded "negative" alignment sets — is not an asset. He called 34 of
them positive. The real negatives number about ten, and six of them came from the
de novo scan.

**The headline for track 2** is the opposite. Yesterday's blind test read as "the
method finds clades, not subfamilies, and there is no way to pick the level." That
diagnosis was incomplete. The level problem is now measured, one of its two halves
is solved, and the reason the other half failed turns out to be a fixable defect in
how the criterion is written.

---

# Track 1 — "is this a SINE?"

## What exists

`calls.tsv` — **113 verdicts on real alignment sets, each with the sentence he
wrote while deciding.** This is the most valuable artefact in the project, because
it records the reasoning and not just the label.

`verdict.py` — a weighted score over four evidence groups: element support 0.45,
homogeneity 0.25, uniqueness 0.20, insertion evidence 0.10. **The weights are
stated, not fitted**, deliberately, to avoid fitting the synthetic negatives.
`calls_scored.tsv` holds 72 sets scored.

## Experiment 1 — do the synthetic negatives work?

42 sets were built to be negative. **He called 8 of them negative. 34 he called
positive.**

| construction | n | he called negative | what he actually said |
|---|---|---|---|
| `NEGSPLICE` | 20 | **0** | *"very clear SINE, almost perfect"* — all twenty |
| `NEGCHIM` | 9 | **0** | *"yes consensus is too long, but its a good sine"* |
| `NEGTRUNC5` | 2 | **0** | *"indeed grayish on the left edge but still more like a SINE"* |
| `NEGMOSAIC` | 1 | **0** | *"absolutely clear SINE"* |
| `NEGLINE` | 1 | **0** | *"these are clear SINEs, what made you think that it is line?"* |
| `NEGSAT` | 3 | 3 | |
| `NEGLOWBIT` | 2 | 2 | *"not informative too obvious negative control"* |
| `NEGSEGDUP` | 2 | 1 | |
| `NEGRAND`, `NEGLINEORF` | 2 | 2 | *"no questions negative, good work"* |

**19 % of the synthetic negatives are actually negative.** The generation strategy
encoded a theory of what breaks a SINE — splice a family in half, chimerise two
subfamilies, truncate the 5′ end — and that theory is wrong. Splicing two
subfamilies together produces something he reads as a perfectly good SINE, because
*it is one*; the copies are real SINE copies and the element is still there.

This has a direct consequence: **the 0.98 binary accuracy is computed against six
negatives, of which two he dismissed as trivial.** Effectively four informative
negatives. That is not a characterised filter.

## Experiment 2 — where do the real negatives come from?

Of the 14 `NOT_SINE*` verdicts, **six come from the de novo scan corpus** (`NEW__hyd_SINE_*`):

- `hyd_SINE_9` — *"CR1-1_HM CR1 Hydra vulgaris (4692 nt) - definitely not sine"*
- `hyd_SINE_7` — *"looks like combination of microsatellites to me, not sine — **not caught by your filters!** — needs new criteria on microsatellite content?"*
- `hyd_SINE_6` — *"mosaic columns in flanks, few copies, not firm edges, unknown middle part"*
- `hyd_SINE_2` — *"too short, few copies"*
- `hyd_SINE_16`, `hyd_SINE_17` — `NOT_SINE_NEEDS_LOOK`

**The de novo scan is the only generator of genuine negatives this project has.**
That is a structural argument for keeping the scan close rather than treating it as
a separate concern: it is the source of the training data track 1 needs.

## Experiment 3 — what does he actually look at?

Term frequency across all 113 verbatims:

| what he mentions | frequency |
|---|---|
| **flank handling — degapping, length** | **21 / 113 (19 %)** |
| mixture of subfamilies / needs splitting | 14 (12 %) |
| alignment quality, presentation | 13 (12 %) |
| length — too long, too short | 12 (11 %) |
| edge / end definition | 10 (9 %) |
| divergence, age | 7 (6 %) |
| consensus wrong or too long | 6 (5 %) |
| too few copies | 6 (5 %) |
| mosaicism | 6 (5 %) |
| nested elements, LINEs | 5 (4 %) |
| microsatellite / tandem content | 3 (3 %) |
| flank uniqueness (needs the genome) | 2 (2 %) |

**The single most common subject of his verdicts is not biology. It is flank
handling.** *"in all cases flanks should have been degapped"*, *"needs much longer
or properly handled (degapped) left flank"*, *"not presented properly with aligned
gappy flanks"*.

## Experiment 4 — how much of the problem is the input, not the element?

**24 of 113 verdicts (21 %) are wholly or partly a complaint about how the set was
prepared**, not about whether the element is a SINE:

`ANCIENT_BADLY_PRESENTED` 4 · `SINE_CONSENSUS_TOO_LONG` 4 · `MIXTURE_SPLIT_FIRST` 3 ·
`UNUSABLE_UNALIGNED` 2 · `NO_FLANKS` 2 · `SINE_CONSENSUS_WRONG` 1 ·
`BADLY_PRESENTED` 1 · `GREY_OR_BADLY_PREPARED` 1 · `SINE_NEEDS_SUBGROUPING` 1 ·
`TOO_FEW_AND_BAD_LEFT_FLANK` 1 · `UNUSABLE` 1 · `SINE_TOO_FEW` 1 ·
`MIXTURE_TWO_SUBFAMILIES` 1 · `SINE_CONSENSUS_TOO_LONG_NOT_MIXTURE` 1

**This is the most actionable finding in track 1.** A fifth of the verdict problem
is a data-preparation QC problem with deterministic answers — degap the flanks,
don't over-extend the consensus, count the copies, check for tandem repeats. It
does not require learning his judgement at all. Fix the preparation and those 24
verdicts change or disappear, which both improves the pipeline and makes the
remaining 89 a cleaner learning problem.

## Honest assessment of track 1

**Usable today:** as a triage filter, human-in-the-loop, with the explicit
understanding that its false-negative rate is unknown.

**Not usable:** as an autonomous accept/reject. Four informative negatives cannot
support that claim, and the threshold of 57 was chosen on the same 49 sets it is
scored on.

**Is it doable?** The binary call — yes, and the gap is data, not ideas. The
28-category verdict is *partly* doable: several categories are deterministic input
checks (see above), several are genuinely visual (`MOSAIC`, `ANCIENT_BADLY_PRESENTED`),
and at least one is **not decidable from the alignment at all** —
`SINE_NEEDS_GENOME_CHECK` explicitly requires proving flank uniqueness genome-wide.
Any honest system has to be able to return "I need the genome for this."

**What would move it fastest:** 30–50 genuinely negative sets, chosen by him or
harvested from de novo scan false positives. Nothing else comes close in value.

---

# Track 2 — "are these one subfamily or two?"

## What exists

The criterion, arrived at by walking six of his peels on saq:

> In the current alignment, a block is a subfamily if it has ≥ 1 column where ≥ 90 %
> of members share a state that ≤ 2 % of everything else carries, and its
> within-minus-between identity gap is clearly positive. Membership of an individual
> chunk is **carriage of those columns — never an identity threshold.**

Supporting code, all ported from his tools rather than invented: `viewer_cons.py`
(MSA-viewer's Copy Consensus, exact), `mafft_dist.py` (MAFFT 6-mer distance,
self-distance 0.0), plus the peel drivers.

## What it does well

Given his boundaries, it reproduces his curation: **three byte-exact consensus
reproductions** (s8, s6, s5 at identity 1.0000) and two near-exact (s2 0.9960,
s3 0.9962), from chunks alone, of consensuses he built months earlier.

It settles that **membership is not similarity**: `input_296` at 0.930 identity to
its own group and 0.960 to the neighbouring group **belongs to the first**;
`input_342` at 0.918 does **not** belong. No identity threshold reproduces both;
diagnostic carriage reproduces both exactly.

## What it does badly

Asked to find boundaries itself on an unseen species, committed in writing before
comparison: **0 of 6** flat, **1 of 6** with recursion added afterwards. Its first
pick was a 214-chunk block spanning three subfamilies whose labels sum to 212.

## Experiment 5 — is the level problem solvable at all?

Test bed: the labelled teu benchmark, six groups of 200 chunks with membership known
by construction. Three classes of candidate block were built —
**SUBFAMILY** (the 6 true groups), **CLADE** (unions of 2 and 3), **FRAGMENT**
(contiguous halves and quarters of a true group) — and each of the four proposed
level-criteria was scored on all 65.

A usable criterion must separate SUBFAMILY from **both** others.

| statistic | non-subfamily blocks falling inside the subfamily range |
|---|---|
| diagnostic count | 97 % |
| diagnostic density per member | 100 % |
| stability under resampling | 100 % — every value zero, the statistic is inert |
| bimodality of the cohesion distribution | 51 % |
| cohesion gap | 42 % |
| **threshold-sensitivity (his own idea)** | **36 %** |

All four fail, but **not symmetrically**, and that is the useful part:

| | clades | subfamilies | fragments |
|---|---|---|---|
| threshold-sensitivity | 0.017 – 0.286 | 0.003 – **0.024** | 0.000 – 0.045 |
| cohesion gap | −0.025 – 0.210 | **0.116** – 0.222 | 0.122 – 0.242 |

**Both detect "too coarse" and both are blind to "too fine."** Fragments are
indistinguishable from subfamilies by every statistic computed on the block —
necessarily so, because half of a homogeneous group is a homogeneous group.

## Experiment 6 — the missing half

If the information is not in the block, it is in what sits nearest to it. A
fragment's nearest neighbours are its own siblings, which lie outside it. A
subfamily's nearest neighbours are its own members.

**Nearest-neighbour containment** — the fraction of members whose single nearest
neighbour in the whole alignment is also a member:

| class | min | median | max |
|---|---|---|---|
| SUBFAMILY | **0.945** | 0.988 | 0.995 |
| CLADE2 / CLADE3 | 0.968 | 0.988 | 1.000 |
| FRAG_half | 0.360 | 0.830 | 1.000 |
| FRAG_quarter | 0.000 | 0.870 | 0.920 |

It separates fragments, which nothing else did. Combined with his
threshold-sensitivity, which separates clades:

> **subfamily ⟺ NN-containment ≥ 0.94 and threshold-sensitivity ≤ 0.025**
>
> **6 of 6 subfamilies accepted. 5 of 59 non-subfamily blocks wrongly accepted.
> Accuracy 0.923.**

The five errors are three clades that slip under the sensitivity bar (`t4+t6`,
`t5+t6`, `t1+t2+t4`) and two fragments whose first half happens to be
self-contained (`t1/2a`, `t2/2a`).

**Both thresholds were fitted on this data.** It is a demonstration that the level
is recoverable, not a validated rule.

## Experiment 7 — why the blind test really failed

In the benchmark, t1 and t2 each have **zero** diagnostics while their union has
**16**. Read naively that says the union is the real unit — that the ground truth is
wrong. It isn't. It is an artefact of the denominator:

| group | nearest sister | diagnostics vs **everything** | vs everything **minus sister** |
|---|---|---|---|
| t1 | t2 (identity 0.886) | **0** | **14** |
| t2 | t1 (0.886) | **0** | **16** |
| t3 | t5 (0.793) | 8 | 20 |
| t4 | t2 (0.880) | 1 | 0 |
| t5 | t6 (**0.922**) | 0 | 0 |
| t6 | t5 (0.922) | 0 | 0 |

**t1 and t2 each have about fifteen private diagnostics that are invisible only
because the other is in the denominator.** The criterion asks *"does this block
differ from everything?"* when a subfamily boundary means *"does this block differ
from its nearest competitor?"* Two sisters therefore cancel each other out and their
union scores high — which is exactly the merge the blind test made.

**This is a formulation defect, not a conceptual dead end.** And it has a
label-free fix: measure `f_out` against the nearest 20 % of non-members rather than
against all of them.

| criterion | subfamilies with ≥ 1 diagnostic | fragments with ≥ 1 |
|---|---|---|
| `f_out` vs everything (current) | **2 of 6** | 0 of 24 |
| `f_out` vs nearest 20 % (proposed) | **4 of 6** | 0 of 24 |

It doubles subfamily detection and still rejects every fragment.

The two it cannot recover are **t5 and t6, whose consensuses are 0.922 identical**
and which have no private diagnostics under any formulation. That is consistent
with the operating range measured on saq: pairs above roughly 0.95 consensus
identity were exactly the ones the method could not separate.

## Honest assessment of track 2

**Usable today:** as an assistant. Given a boundary it will say whether the block
holds up, which chunks fail it, what the diagnostic columns are, and it reproduces
consensuses exactly. That is real value in a human-in-the-loop workflow.

**Not usable:** unsupervised, on a new genome.

**Is it doable? Yes — and I am more confident than I was a day ago**, for three
specific reasons rather than optimism:

1. The failure was **diagnosed to a named defect** with a label-free fix, and the
   fix doubles detection on the benchmark.
2. The half of the problem that looked hardest — telling a subfamily from a fragment
   of itself — turned out to be **solved by a statistic that takes ten lines**
   (nearest-neighbour containment), and it is the half no statistic in the original
   four addressed.
3. A two-statistic rule reaches **6 of 6 with 8.5 % false accepts** on labelled
   data. Fitted, but the signal exists to be fitted to.

**What is genuinely hard, and may not be solvable:** pairs whose consensuses are
≥ 0.92 identical (t5/t6 here; s7g/s1/s9 on saq at 0.951–0.973) have no private
diagnostics under any formulation tried. For those, either a different kind of
evidence is needed — copy age distributions, insertion sites, genome context — or
they are genuinely one subfamily and his split is a convention. **Nothing measured
so far distinguishes those two possibilities**, and that is the honest limit.

## What would move track 2 fastest

1. **Re-run the blind test with local `f_out` and the two-statistic level rule.**
   Everything needed exists; it is a day's work and it is a real test, since the
   rule was fitted on the per-subfamily file and the blind test uses the
   genome-wide one.
2. **A second walked peel on one species**, so the level rule has more than n = 6 to
   stand on. Keep `ccr` or `dmo` untouched for the blind test afterwards.
3. **Decide the t5/t6 question deliberately** — take one pair he considers separate
   whose consensuses are ≥ 0.92 identical, and ask what evidence *he* uses. If it is
   not in the alignment, the method needs another input and we should know that now.

---

# Bearing on the manuscript

Neither track is ready to appear as a **method**. The post-processing pipeline —
known consensi to a full genomic repertoire with per-copy description — is validated
across six talpids and nine scorpion genomes and stands alone.

These two schools belong in that paper, if at all, as a stated limitation and a
roadmap: **subfamily boundaries are currently expert-curated.** That is accurate,
normal for a tool paper, and much safer than implying discovery is automated.

One caveat that cuts the other way, and is worth weighing: **the de novo scan is the
only source of genuine negatives track 1 has.** If the scan is cut out of the
manuscript entirely, the natural pipeline for generating the data both schools need
gets cut with it. That is an argument for at least citing it as a companion rather
than severing it.

---

# The single biggest weakness, across both tracks

Every number in this project except one was obtained on data that helped produce it.
The one exception is the teu blind test, and it scored 1 of 6.

That is fixable only by holding data back — and there are two untouched species left
to hold back with. Spending them on validation rather than on more calibration is,
I think, the most important methodological decision facing this work.
