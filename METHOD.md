# How Sergei tells a SINE from a non-SINE

The method the tool has to reproduce. Steps 1–4 were dictated directly on
2026-09-01. Steps 5 onward are **reconstructed from the record** — his review
messages across this project, quoted verbatim where possible — and are marked as
such so he can correct rather than re-dictate.

---

## The goal

> "I was looking for a way to tell sine from non-sine."

And the shape of the solution, stated more than once and worth keeping at the
top because it is what everything else serves:

> "replace the expert in analyzing candidate sine by expertly looking at
> alignments"

The tool is not a search tool and not a library builder. **It replaces the
judgement a person makes when reading a 100-locus alignment.** Everything below
is a description of that judgement.

---

## Step 1 — Start from a preliminary candidate consensus *(dictated)*

> "I am using preliminary candidate consensus to search"

It is *preliminary*: it may be wrong, and whatever is wrong with it shapes
everything downstream. Too long → the hits carry flanking junk. A chimera → the
hits are two different things. A fragment → only part of the element is found.

Sources of that first consensus: a de novo scan (SINEderella), another tool's
output (AnnoSINE2, RepeatModeler, HiTE), an existing library entry that needs
adjusting, or a hand-built cluster consensus.

## Step 2 — Search the genome with it *(dictated)*

Hits come back. The count is not fixed — thousands, hundreds, tens, or zero.

## Step 3 — Inspect the hits *(dictated)*

> "after this search returns hits, i inspect the hits"

## Step 4 — Inspection means ~100-locus alignments with flanks *(dictated)*

> "i do some 100-loci alignments with flanks and look at them"

- **~100 loci** per alignment — what can actually be read by eye
- **Several** such alignments, not one
- **With flanks**, because the flanks are where the element's edge shows itself

What the several alignments are *(answered 2026-09-01)*:

> "random samples, best top hits, everything which works"

So: a random sample of the hits, and the top-scoring hits, and any other view
that helps — they are different windows on the same hit set, not different
candidate families. This is exactly the `top100` / `rand100` / `subfam` split
already present in the Tal data.

---

## Step 5 — What is looked for in that alignment *(reconstructed)*

### The decisive signal: the flanks must NOT align

> "its flanks are good — not alignable"

Copies of a real SINE sit in unrelated genomic contexts. Similarity is high
across the element and collapses to background immediately outside it. That
collapse — sharp, on both sides — is the signature of an insertion.

If the flanks *do* resemble each other, the copies are not independent
insertions: a satellite, a duplicated region, or copies sitting inside one
larger repeat.

**Consequence for display:** flanks must be **de-gapped and pushed against the
element** to be judged. An aligner will happily align unrelated flanks and make a
clean insertion look ragged.

> "not presented properly with aligned gappy flanks, need proper work"

### The edges, judged independently

Left and right are separate questions and one can be fine while the other is not.

> "right end is fine but left one is really bad, good example of grey zone which
> looks more like LINE"

> "right end is a bit wobbly unstable" → grey zone

> "if its right edge is defined at gggagat… it becomes quite clear"

That last one matters: **the boundary is part of the judgement, not a later
step.** Re-defining the edge can turn an unclear set into a clear one.

### Copies must agree with each other across the element

Above genomic background. Note this is copy-to-copy agreement, which is not the
same as agreement with the supplied consensus — the consensus may be wrong.

---

## Step 6 — What does NOT disqualify *(reconstructed)*

- **High divergence.**
  > "no caveats, pure clear SINE with high divergence"
- **Gaps and insertions in the tail**, if the right edge consolidates once gaps
  are removed.
  > "many individual level long insertions in the tail region but if gaps are
  > removed it becomes consolidated at right edge… its clearly SINE"
- **Truncated copies** — still copies of that family.
  > "indeed grayish on the left edge but still more like a SINE"
- **Two subfamilies in one set.**
  > "subfamily ambiguity doesn't matter, we answer sine/not-sine, not this sine
  > or that sine"
- **A minority of bad copies.** A set can be a clear SINE with 11 problem
  sequences worth noting.
  > "11 lower sequences need manual attention… but overall this is a SINE no
  > doubt in its majority"

---

## Step 7 — Answers other than yes or no *(reconstructed)*

- **A mixture of two different things** — not a bad set, an unfinished one.
  > "truly a mixture… top proper about half longer similarity sequences need
  > separate analysis as a good SINE candidate, then bottom more discordant
  > shorter ones need additional fresh re-run of whole analysis. This
  > non-homogenous selection problem is different from bad SINE set and needs
  > additional work before verdict."
- **Too little evidence** — a verdict on a handful of copies is an impression.
  > "very weak signal over too few copies, very low chance being a SINE"
  > "too few copies but looks like SINE with huge discordant insertions, needs
  > more evidence from genome"
- **Real but needs re-checking.**
  > "can be SINE, but very turbulent at many places, should be regarded as SINE
  > but with caution and need manual reinspection"
- **Badly prepared rather than bad.**
  > "unclear situation with left end, could be just not extended enough, more
  > like grey zone or technically badly prepared"
- **Consensus wrong rather than family absent.**
  > "weak old sine, consensus needs refinement and shortening, otherwise
  > difficult but very real SINE"

---

## Step 8 — What cannot be decided inside one alignment *(reconstructed)*

**Flank uniqueness has to be proven against the whole genome**, not within the
set:

> "nested copies, i agree, but there are enough uniq left flanks to lean towards
> SINE — requires post-processing with proving uniqness of at least some flanks
> on whole-genome level"

## Step 9 — Then: refine and repeat *(reconstructed)*

If the set looks good, build a consensus from it and search again with that —
each pass sharpening the consensus and the boundaries.

For the SubFam route specifically, the loop is a **peel**:

> "I do this by looking at alignment, pick up groups of similar sequences,
> copy-pasting them to separate alignment and if it looks good then after
> realignment I create their consensus and use it for future re-scan and
> refinement in next sinederella run. After selecting more obvious groups in
> subfam, I continue with removing them from subfam alignment and continue with
> remaining ones, until no groups remain (only noisy ungroupable sequence)."

With the correction that the group consensus must be built from **the original
member sequences**, not from the chunk-consensuses:

> "the algorithm should not create consensus from these 'consensi of 50s' used
> in subfam but look back into those original 50sequence chunks"

---

## Marked uncertain — for Sergei to correct

2. **Step 5: is there an order?** Does he look at flanks first, edges first, or
   take the whole picture at once?
3. **Step 9: when does he stop iterating?** When the consensus stops changing,
   when the alignment looks clean, or on some other signal?
4. **Zero hits, or ten hits** — what does he do then? Abandon the candidate, or
   is a small set still worth pursuing?
