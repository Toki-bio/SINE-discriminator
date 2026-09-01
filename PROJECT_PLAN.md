# SINE discriminator — full project plan

Replaces `ROADMAP.md`, which was framed as if from a standing start and written
in software jargon. This is written from the accumulated record: `HANDOFF.md`
(59 sections), the three benchmark legs, Sergei's own eye reviews, and the
retractions. Nothing here is new invention — it is the project as it now stands.

---

## 1. What the project is

Sergei finds and curates SINE families in genomes. A search returns anywhere
from **zero to a few thousand** genomic hits, and that set is a mess: some hits
are real copies of a real family, some are junk, several families are often
mixed together, and the element boundaries are usually wrong.

Sorting that out is currently done by eye — read the alignment, recognise
groups, pull one out, build a consensus, search again with it, repeat until the
families are clean. The purpose of this project is to make that judgement
**measurable**, so a computer performs it and a person checks it.

It serves two users at once: **Sergei's own curation**, and **other people
running it on their own genomes**. Both, not one.

### The three starting situations

This is the branch point for everything the tool does, and getting it wrong is
what made the earlier plan read as alien.

| situation | what the tool must do |
|---|---|
| **A. No SINE library exists** for this genome | discover families from scratch: group the hits, decide which groups are real, fix their boundaries, emit consensuses |
| **B. A library exists but needs adjustment** | take each existing consensus, check whether the copies actually support it, correct the boundaries, split mixtures, prune contamination |
| **C. A library exists and is trusted** | use it as ground truth: calibrate and validate the criteria against known-correct answers |

Situation **C is what has been done so far** — human Dfam, Sergei's curated
`tim/`, and the Timema AnnoSINE benchmark were all used to calibrate. That work
is what makes A and B trustworthy, and it is not a detour from the project; it
is its foundation.

---

## 2. What "is this a SINE?" actually means

Taken from Sergei's own judgements across this project. These are the criteria
the tool must reproduce, and every one of them has already changed the code.

**Positive signals**

- **The flanks must not align.** *"its flanks are good — not alignable."* Copies
  sit in unrelated genomic contexts; that is the signature of an insertion. The
  flanks must be de-gapped and pushed against the element to be judged, never
  aligned.
- **A sharp edge on each side**, judged independently. *"right end is fine but
  left one is really bad"*; *"if its right edge is defined at gggagat… it
  becomes quite clear."*
- **Copies agree with each other** across the element, above genomic background.

**Things that do NOT disqualify**

- **High divergence.** *"no caveats, pure clear SINE with high divergence."*
- **Gaps and insertions in the tail**, if the right edge consolidates once gaps
  are removed. *"many individual level long insertions in the tail region but if
  gaps are removed it becomes consolidated at right edge… its clearly SINE."*
- **Truncated copies.** *"more like a SINE"*, never *"not an element."*
- **Two subfamilies in one set.** *"subfamily ambiguity doesn't matter, we answer
  sine/not-sine, not this sine or that sine."*

**Things that need a different answer than yes or no**

- **A mixture of two different things** — *"needs additional work before
  verdict."* Split it and re-run each part; do not score it.
- **Too few copies** — *"very weak signal over too few copies."* A verdict on
  n=6 is an impression, not a measurement. n=1 got *"wtf???"*.
- **Real but turbulent** — *"should be regarded as SINE but with caution and need
  manual reinspection."*
- **Badly presented rather than bad** — ancient families whose *"gappy flanks"*
  make them unreadable until de-gapped.

**What cannot be judged inside one alignment**

- **Flank uniqueness must be proven genome-wide.** *"enough uniq left flanks to
  lean towards SINE — requires post-processing with proving uniqness of at least
  some flanks on whole-genome level."*

**The boundary is part of the verdict, not a later step.** *"consensus needs
refinement and shortening"* was the correct call on a family the tool was
rejecting outright.

---

## 3. What exists and is proven

### The verdict engine — finished

Given an alignment with the consensus as the first row, it returns a 0–100
score, named sub-cases instead of a binary, and withholds judgement on mixtures.

Calibrated against four independent sources:

- **Synthetic sets with known truth.** Sets known *not* to be SINEs never score
  above **47**. Sets known to *be* SINEs never score below **90**. Nothing lands
  between — the scale separates cleanly with room to spare.
- **Sergei's curated `tim/` subfamilies** — all 42 alignments accepted, and the
  100-best-copy version of all 14 subfamilies scored exactly 100. Their
  boundaries had been independently confirmed on all 24 sides, so this tests the
  tool against boundaries someone else verified.
- **64 human Dfam families** — 61 accepted.
- **RepeatMasker's own human annotation** (1.9 million elements, never seen by
  this tool) — **95.6 %** of the 5,306 benchmark copies land on an annotated
  SINE; 74.2 % carry the same family name.

### Group discovery from SubFam — validated, loop not yet run

Sergei's manual "pick a group, peel it off, repeat" loop. Automated clustering of
the same 161 chunk-consensuses he clustered by hand, with the group consensus
built from **the original member sequences** rather than from the
consensi-of-50s (his correction), **recovered 7 of his curated subfamilies
blind** — `t2` at 99.5 %, `t3-1` at 100 %, `t1_1` at 99.2 %, `t8-1` at 98.5 %,
`t6-1` at 98.0 %, plus `t6-3` and `t6-4`.

### Boundary handling — half

- **Left edge: solved.** Detected within 1–3 bp of the true start on every
  known-good set.
- **Inward trim exists** but only reports; it does not rewrite the consensus.
- **Outward extension is detected and not acted on.**
- **Right edge: the obvious method is wrong** — see §5.

### Presentation — benchmark pages live

`benchmark.html` (three legs, filters, every row linking to the alignment at
three geometries) and `index.html` (36 alignments with per-position graphs).

---

## 4. What does not exist

1. **Nothing takes a bare consensus and a genome and produces a scored
   alignment.** Every result so far used an alignment somebody else built. This
   is the single biggest gap between "a method" and "a tool".
2. **Nothing writes a corrected consensus out.** The tool says the boundary is
   wrong; it never fixes it.
3. **No repeat-until-stable loop.** Correct the boundary, re-search, re-align,
   re-score, until it stops moving.
4. **No right-edge rule.**
5. **Genome-wide flank uniqueness is measured but useless so far.** It was built
   and calibrated (median 0.975 unique across 64 human families) but does not
   separate the case it was built for — `LmeSINE1c`, a false positive, scores
   0.917 while a real family scores 1.000. Its copies genuinely are in unique
   DNA; they are unrelated loci a weak consensus aligned. Needs rethinking, not
   more tuning.
6. **No reading of AnnoSINE2 / RepeatModeler / HiTE / RepeatMasker output.**
   Each writes a different file format; something must read each one.

---

## 5. The one genuine unknown

**Where does the element end on the 3′ side?**

The obvious method — move the edge to where copy-to-copy similarity stops —
fails on every SINE, and fails the same way. Copy identity collapses in the
poly-A tail because the A-run length differs between copies:

| | element body | last 90 bp |
|---|---|---|
| AluY | identity 0.93, A+T 0.26 | identity **0.27**, A+T **0.63** |
| curated t1_1 | identity 0.96 | identity **0.31**, A+T **0.70** |

So an identity-based edge finder amputates the poly-A tail from every SINE it
touches — deleting a defining feature. The rule must also watch base
composition: identity falling **while A+T rises** means still inside the tail;
identity falling **at ordinary genomic composition** means the element has
ended.

The pieces exist — terminator and A-run detection, and an A+T track per
position. The rule does not.

**This is the go/no-go for the whole refinement half of the project**, which is
why it is step 1.

---

## 6. Order of work

**Step 1 — the 3′ edge rule.**
Test: on all 64 human Dfam families, whose true ends are known, the corrected
edge must land within a few bases of the curated end **without cutting off the
poly-A tail**. If this cannot be made to work, boundary correction cannot be
built as conceived and the plan changes shape.

**Step 2 — run the peel loop, both ways.**
Written, never executed. Run it with mixtures treated as blocking (split and
re-cluster, converging toward the finer 14 subfamilies) and as acceptable
(converging toward the coarser 8). Compare both against `t1`–`t8` and v4; keep
whichever recovers more real families with fewer false ones.

**Step 3 — consensus + genome → scored alignment, as one command.**
Search, extract with flanks, align, score. Assembly of existing parts. This is
what makes it a tool someone else can run.

**Step 4 — the repeat-until-stable loop.**
Test: start from a deliberately wrong boundary and converge on the curated one.
Needs a fixed reference point and a hard iteration cap, or it can drift onto a
subfamily or onto nothing.

**Step 5 — read the upstream formats** (AnnoSINE2, RepeatModeler, HiTE,
RepeatMasker), so situation A and B inputs both work.

**Step 6 — the results page** (§7).

**Running in parallel:** AnnoSINE2 on three genomes with no SINE library —
*Pomacea canaliculata* (snail), *Acanthaster planci* (starfish), *Hydra
vulgaris* (cnidarian). First prospective test: three phyla, no known answer.

---

## 7. The results page

One page per run, published to the site, holding:

1. **A diagram of the run** — the stages with the counts that flowed through
   each: candidates in → loci found → groups formed → accepted → families out.
2. **The families found** — one row each: consensus length, copy count, score,
   flags, links opening the alignment at 50 bp / 400 bp / justified.
3. **A boundary picture per family** — identity and A+T along the sequence, with
   the supplied edge and the corrected edge both marked, so a boundary change
   can be judged by eye in one glance.
4. **What was left over** — groups that formed no acceptable family, so nothing
   disappears silently.
5. **Provenance** — genome, upstream tool, parameters, date.

---

## 8. Honest state

| part | built | proven |
|---|---|---|
| verdict engine | yes | yes, four independent ways |
| group discovery | clustering + consensus | 7 subfamilies recovered blind |
| left edge | yes | within 1–3 bp |
| right edge | no | obvious method disproven |
| consensus rewriting | no | — |
| repeat-until-stable | no | — |
| genome + consensus → alignment | parts only | — |
| reading upstream formats | no | — |
| genome-wide flank uniqueness | built | does not discriminate |
| results page | benchmark pages live | rendered, verified |

**The judgement is ~85 % done. The tool around it is ~40 %.** The hard,
original part — deciding whether a locus set is a family, with named reasons —
is finished and tested. The machinery that turns it into something runnable is
not.

**Two thresholds are fitted to very few examples** and must be revisited before
anything acts on them automatically: the mixture test (fitted to three cases)
and the over-long-consensus test (fitted to one).
