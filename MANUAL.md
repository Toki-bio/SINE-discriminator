# How the SINE discriminator works

A plain description of what the tool actually computes today, which variables
combine into the verdict, over what ranges, what has been shown to work, and
what has not. Everything here is read off the code, not remembered.

Companion documents: `HANDOFF.md` is the full working history including every
retraction; this file is the summary.

---

## 1. The question it answers

Given a set of genomic loci proposed as one family — typically the output of a
`sear`/`blastn` search from a consensus — decide whether it is **a real SINE
family**, and if not, say which of several different things went wrong.

It deliberately does **not** answer "which subfamily is this". Sergei's rule:
*"subfamily ambiguity doesn't matter, we answer sine/not-sine, not this sine or
that sine."* Two subfamilies mixed in one set is still a set of SINEs, so the
subfamily split is reported as a note and costs no score.

The output is a **0–100 score, a list of named flags, and sometimes a withheld
verdict** — not a yes/no.

---

## 2. The input it needs

One aligned FASTA per candidate family, in which **the first sequence is the
known consensus**, named `CONSENSUS_*`.

That requirement carries the single most important design decision in the tool:

> **Boundaries come from the consensus row, not from a threshold.**
> The element is exactly the span between the consensus's first and last
> nucleotide. Every earlier attempt to *find* the boundary from the data needed
> a tuned cut-off and got it wrong; anchoring on a row that is already in the
> alignment has **zero free parameters**.

Each copy's element length is then simply the number of its own bases inside
that span.

---

## 3. The four evidence groups

| group | weight | question |
|---|---|---|
| **element** | 0.45 | is there an element at all, and do the copies support it? |
| **homogeneity** | 0.25 | is it *one* element, uniform across copies? |
| **uniqueness** | 0.20 | are the copies in independent genomic locations? |
| **insertion** | 0.10 | are there TSDs — positive evidence of an insertion event? |

Each is mapped to 0–1 by a saturating ramp `sat(x, lo, hi)`: below `lo` it is 0,
above `hi` it is 1, linear between.

### element (the gate)

```
g_elem = sat(id_all − flank_bg, 0.10, 0.30)^0.5 × sat(n_supported / n, 0.15, 0.98)^0.5
```

Two factors, both necessary: a **cliff** (copies match the consensus much better
than unrelated DNA does) and **breadth** (most copies participate).

`flank_bg` is the identity between copies *outside* the element, i.e. genomic
background — about 0.25 for unrelated DNA.

The 0.10–0.30 range is deliberately low. At 30 % divergence a real family sits
at identity 0.57 against a 0.29 background; that is already unambiguous.
Scoring on up to 0.55 made the statistic report **youth rather than existence**
and cost a clean diverged family a third of its element score.

### homogeneity

```
g_homog = 1 − sat(cv_core, 0.10, 0.50)
```

`cv_core` is the coefficient of variation of copy length. The range is
deliberately forgiving: truncated copies are still copies of that family.
Sergei on the truncated sets: *"more like a SINE"*, never *"not an element"*.
Scored on 0.05–0.25 it drove the truncated class to a mean of 9.7, which was
wrong.

### uniqueness

With 400 bp flank-decay data:

```
g_uniq = (1 − sat(decay_max, 75, 300)) × (1 − sat(edge_max, 0.45, 0.85))
```

**Both halves matter.** Distance alone rates a LINE fragment as isolated,
because its shared flank runs only ~50–75 bp before the copies truncate — yet
identity right at the edge is 0.89, which says the element does not end there.

Without decay data: **`g_uniq = 1.0`, and the absence is reported.** See §6.

### insertion

```
g_ins = sat(tsd_frac, 0.10, 0.60)
```

TSD presence only ever *adds*. Its absence is never held against a set, because
old copies lose their TSDs.

### assembly

```
rest  = g_homog^0.55 × g_uniq^0.45
score = 100 × min(1, g_elem × rest × (1 + 0.10 × g_ins))
```

**The element group gates the score rather than voting in it.** As a co-equal
term, a set with no element (g_elem 0.31) still scored 56 because its flanks
were unique and its lengths uniform — both true, and both irrelevant when
there is nothing there.

---

## 4. The flags — what the tool says instead of just a number

| flag | fires when | means |
|---|---|---|
| `NO_ELEMENT` | g_elem < 0.25 | no cliff: copies match each other no better than background |
| `CONTAMINATED` | bimodal identity, both modes real | a real family with junk mixed in — prunable |
| `SUBFAMILY_NOTE` | gap_excess > 0.03 | two subfamilies present; **does not affect the verdict** |
| `HETEROGENEOUS_SELECTION` | group median lengths differ ≥ 0.12 | two *different things* collected together — **verdict withheld** |
| `CONSENSUS_OVEREXTENDED` | a core window beats the whole span by ≥ 0.15 and the remainder sits at background | the consensus is longer than the element; trim to the reported window |
| `TRUNCATED_COPIES` | cv_core > 0.18, n_core ≥ 15 | many 5′-truncated copies |
| `FRAGMENT_OF_LONGER` | decay says similarity continues | the element runs past its annotated boundary |
| `SHARED_FLANKS` | flank identity − 0.25 > 0.15 | copies are not in independent contexts |
| `NOT_ISOLATED` | decay_max large | satellite, duplication, or host repeat |
| `INSUFFICIENT_COPIES` | n < 30 | score is an impression, not a measurement |
| `FLANKS_UNMEASURED` | no decay data, apparent flank sharing | reading unreliable at this flank width |
| `NO_FLANKS_PRESENT` | no measurable flank sequence | isolation/nesting/TSD unavailable |
| `RECOVERABLE_CORE` | n_core ≥ 20 but score < 55 | enough good copies to rebuild a clean family |

**A withheld verdict is a third answer.** `HETEROGENEOUS_SELECTION` sets a
`deferred` field: scoring a mixture answers neither question, so no number is
offered until it is split.

---

## 5. The two geometries — and the rule that follows

The same loci are scored on two different alignments:

| flanks | good for | bad for |
|---|---|---|
| **~50 bp** (as published) | how well copies match their consensus | anything about boundaries or context |
| **400 bp** (re-extracted) | where similarity *stops* | element identity in ancient families |

**Rule: element evidence from the short alignment, flank decay from the long
one. Never both from one geometry.**

This is measured, not assumed. At 400 bp the four ancient human MIR families
lose 19–37 points, entirely in the element group, while `flank_bg` stays flat —
the consensus smears from 1.59× to 2.41× its own length, because MAFFT cannot
absorb 400 bp of unalignable flank at 55 % identity. Young Alus at 88 % identity
are unaffected and actually *tighten*.

Conversely, without 400 bp flanks the isolation question cannot be asked at all:
every Timema candidate looked "nested" at 50 bp and every one of them proved to
be `ISOLATED_INSERTION` at 400 bp.

---

## 6. The one principle that caused three of the bugs

> **An absent measurement must never be scored as guilt.**

Three separate defects were the same mistake:

1. **Uniqueness fallback.** With no decay data, a 50–70 bp flank-sharing test
   stood in. Measured across the corpus it had *zero* discriminating power —
   all 30 genuine negatives score 0.0 with and without it — and its only effect
   was suppressing real SINEs. Removed; uniqueness is now neutral when
   unmeasured.
2. **Over-extended consensus.** `MIR1_Amn` scored 0.0 `NO_ELEMENT` because its
   consensus ran past the element into unrelated sequence, dragging identity
   down along the whole span. Support is now re-measured on the core window —
   judging support against a boundary just shown to be wrong is circular.
   0.0 → 100.0.
3. **Flank background from absent flanks.** The curated `subfam` alignments
   carry no genomic flanks; a handful of stub bases matched trivially and the
   background came out as **1.000**, which makes the cliff negative by
   definition. `t3-1`/`t3-2` scored 0.0 while their own `top100` scored 100.0.
   Both background computations now require ≥ 20 copies with ≥ 15 bp of flank.

---

## 7. What actually works — measured

| corpus | n | mean | extreme | crossing 50 |
|---|---|---|---|---|
| synthetic true negatives | 30 | 1.59 | max 47.0 | **0** |
| synthetic true positives | 432 | 99.96 | min 90.1 | **0** |
| **curated `tim/` ground truth** | **42** | **98.29** | min 71.3 | **0** |
| human Dfam | 64 | 94.49 | | 2 |
| Timema AnnoSINE candidates | 55 | 94.48 | | 3 |

A 43-point empty gap between the negatives and the positives, and a clean sweep
on the manually curated set whose boundaries were independently confirmed on all
24 sides.

Independent check: **95.6 % of the 5,306 human benchmark copies overlap a
RepeatMasker hg38 SINE annotation** (1,910,631 elements), 74.2 % by name.

### What works best

- **The cliff (element gate)** is the load-bearing statistic. Every genuine
  negative is rejected by it alone; nothing else is needed to reject them.
- **Consensus-row anchoring** — zero parameters, and it survived every corpus.
- **Contamination as bimodality**, not as a low tail. A diverged family has a
  unimodal identity distribution; a contaminated one has two modes. Both modes
  must be substantial: guarding only the lower one let a single copy at identity
  1.0 stand as an entire "clean family" and label the other 31 contamination.
- **400 bp decay** for isolation — the only measure that separates a real
  insertion from a satellite, and it is unanimous on both benchmark species.

### What did not work, and was removed

- **Turbulence / `NEEDS_REVIEW`.** Counted *columns* dipping below the family
  median. Easy to trigger at identity 0.88, impossible at 0.55 — so it was most
  sensitive where families are most homogeneous. It fired on 129 sets, 91 as the
  only flag, all scoring 90.1–100, and never moved a score. Every real problem
  it was supposed to catch was **row-wise** (subsets of copies), which
  `CONTAMINATED` and `SUBFAMILY_NOTE` already detect.
- **Within-set flank sharing** as a uniqueness proxy (see §6).
- **The gap rule** for contamination pruning — it manufactures families from
  noise. The MAD rule is used instead.

---

## 8. Flank adjustment — what the tool does and does NOT do

This deserves a direct answer because it is easy to assume more than exists.

**There is no adaptive search over flank widths.** The tool does not try several
flank lengths, score each, and keep the best. Flank width is a property of the
input alignment. Two fixed geometries are prepared *outside* the scorer
(~50 bp as published, 400 bp re-extracted from the genome by `bench_extend.py`),
and each is scored independently.

**One boundary adjustment does happen, and it is inward only.**
`core_window()` searches for a contiguous sub-span of the consensus that the
copies support much better than the whole:

- window length from 50 bp up to 80 % of the consensus, step 5
- start position step 5
- keeps the window maximising `inside − outside` identity
- fires only if `inside − whole ≥ 0.15` **and** `outside ≤ background + 0.18`

The second condition is what makes it a boundary test rather than a general
"some window is better" test — in any structured profile some window always
beats the mean, and a first version without it flagged clean `AluJr`, `AluSc`
and `AluY`. The `+0.18` admits weak residual similarity (`MIR1_Amn`'s tails sit
at 0.395 against 0.25 background) while staying far below any clean family's
~0.85.

When it fires, per-copy identity and support are recomputed **on that window**,
and the trim coordinates are reported so the consensus can be shortened.

**What is missing:** the symmetric outward case. If an element extends *past*
its annotated boundary, the tool detects it (`FRAGMENT_OF_LONGER`,
`ELEMENT_CONTINUES` from the decay curve) but does not attempt to extend the
consensus and re-score. That would be the natural next step and is not built.

---

## 9. Known limitations

1. **Whole-genome context is barely used.** A flank-uniqueness check exists
   (map 300 bp flanks back with `bwa`, read MAPQ) and is calibrated — median
   0.975 unique across 64 human sets — but it does **not** separate the case it
   was built for. `LmeSINE1c` (a false positive at 69.2, 7.1 % RepeatMasker
   overlap) scores 0.917, `AmnSINE1` scores 1.000. Its copies really do sit in
   unique DNA; they are unrelated loci a weak consensus aligned. No internal
   statistic catches this — external annotation does.
2. **Two thresholds are overfit** and should be revisited, not trusted:
   `HETEROGENEOUS_SELECTION` at 0.12 (fitted to three examples) and the
   `CONSENSUS_OVEREXTENDED` tail bound at `bg + 0.18` (fitted to one).
3. **Ancient families measure worse at 400 bp** (§5) — quantified, not
   compensated.
4. **The 40–75 band is not calibrated.** With a 43-point empty gap the exact
   threshold has not mattered; it will as grey-zone cases accumulate.
5. **Scope.** Validated on talpid SINEs, human Dfam, and Timema. A prospective
   test on three phyla with no SINE library (Mollusca, Echinodermata, Cnidaria)
   is running and has not reported yet.

---

## 10. Running it

```bash
python verdict.py <alignment_dir> [out.json]     # score a corpus
python flankdecay.py <ext_alignment_dir>          # 400bp decay -> flankdecay.json
python measure_c.py <alignment>                   # raw statistics for one set
python build_viewer_data.py                       # pack the report page payload
python rebuild_embed.py                           # regenerate _embed.js
```

`flankdecay.py` **overwrites** `flankdecay.json` rather than merging — a run on
one species will silently drop another's entries. Merge by hand, or the talpid
data disappears.
