# SINE discriminator — the plan

Structure: **problem → task → testing → propose the algorithm → build the tool.**

---

## 1. The problem

Measured on 72 of Sergei's own calls (`calls.tsv`, `calls_scored.tsv`):

**(a) The output shape is wrong.** He produces categories with an action
attached. The tool produces one number, and his categories collapse onto it:

| his call | n | mean score |
|---|---|---|
| SINE_MINOR_NOTE | 4 | **100.0** |
| SINE_CONSENSUS_WRONG | 1 | **100.0** |
| SINE_NEEDS_SUBGROUPING | 1 | **100.0** |
| SINE_MOST_PROBLEMATIC | 1 | **100.0** |
| SINE | 43 | 99.1 |

Four situations needing four different follow-up actions — shorten the
consensus, subgroup it, note the instability, do nothing — all score identically
to a clean SINE. Likewise MIXTURE_SPLIT_FIRST (79.6), GREY (89.5) and
ANCIENT_BADLY_PRESENTED (92.6) overlap completely while demanding opposite
responses.

**(b) Four outright failures.**

| his call | set | score |
|---|---|---|
| MOSAIC | ERI e1-4 | **0.2** |
| SINE_TOO_FEW | AluSx_short_ | **0.1** |
| NOT_SINE | LmeSINE1c | **69.2** |
| SINE (lowest) | — | 68.3 |

---

## 2. The task

Two things, in this order:

1. **Decide what the tool outputs.** Not a score — a **category plus the action
   it implies**. The categories are his, already observed in his own words.
2. **Find which measurable variables separate those categories**, and admit which
   ones no current variable can separate.

---

## 3. Testing

The calibration set is 72 calls across 16 categories. Everything below is
measured against those, never against synthetic truth.

### Test 1 — What separates his categories? *(all variables, existing calls)*
Compute **every** variable the project has on all 72 alignments — not just the
four groups that feed the score, but the raw ones: per-position identity and
coverage profiles, A+T, edge sharpness left and right, copy-length spread,
subfamily gap, contamination gap, flank decay where available, TSD fraction,
copy count, consensus-vs-copies agreement.
Then, per category pair, ask which variables differ.

**Output:** a table of variable × category with effect sizes.
**Decision:** which categories are already separable, and which are not.

### Test 2 — What is missing? *(the categories nothing separates)*
For every pair Test 1 cannot separate, look at those alignments directly and
name what a person sees that is not being measured. Sergei's own phrases point
at candidates already: *"11 lower sequences"* (a subset of copies, not a
column), *"island of instability"* in a specific region, *"right edge defined at
gggagat"* (a motif at the boundary), *"gappy flanks"* (a display property).

**Output:** a short list of proposed new variables, each with the case that
demands it.
**Decision:** which are worth implementing.

### Test 3 — Do the new variables work? *(implement and re-measure)*
Implement the shortlist, re-run Test 1.

**Output:** updated separation table.
**Gate:** a new variable earns its place only if it separates a pair that was
previously inseparable, without disturbing pairs that already worked.

### Test 4 — More calls, chosen where we are least sure
The 72 are unevenly spread: 43 are plain SINE, and nine categories have exactly
one example. One example cannot support a threshold.

Pick alignments where the proposed classifier is least certain, show them to
Sergei, record his calls. **Ask for judgements only where they change something**
— never as a general request to re-label.

**Output:** an expanded calibration set, weighted toward the thin categories.

### Test 5 — Held-out validation
Fit thresholds on part of the set, test on the rest. Report how often the
predicted category matches his, and where it fails.

**Gate:** if the classifier only works on the alignments it was fitted to, it is
not a tool.

---

## 4. Propose the algorithm

Written only after Tests 1–5, and stating for each category:

- the variables that identify it, with their ranges
- the confidence, given how many examples support it
- **the action it implies** — split, re-extract with longer flanks, shorten the
  consensus, gather more copies, check flanks genome-wide, accept, reject
- what it must NOT be confused with, and what prevents that

Plus, explicitly: the categories that remain unseparated, so they are known
limitations rather than silent errors.

---

## 5. Build the tool

Only after the algorithm is proposed and validated.

1. **One entry point**: candidate consensus + genome → hits → several ~100-locus
   alignments with flanks (random sample, top hits) → category + action per
   alignment. This is his steps 1–4, automated.
2. **The actions that can be automated**, each already half-present:
   split a mixture and re-run each part; shorten an over-long consensus;
   re-extract with longer flanks; prune contamination; de-gap flanks for display.
3. **The peel loop** for the SubFam route — written, never run; group consensus
   built from original members, not from consensi-of-50s.
4. **The iteration**: refine consensus → re-search → re-judge, until stable.
5. **The results page**: per run — a diagram of what flowed through, the
   families found with category and action, a boundary picture per family
   showing supplied vs corrected edges, the residue, and provenance.

---

## 6. Order and gates

| # | step | gate to pass before continuing |
|---|---|---|
| 1 | Test 1 — separation table | — |
| 2 | Test 2 — name what is missing | — |
| 3 | Test 3 — implement, re-measure | new variables must separate something previously inseparable |
| 4 | Test 4 — more calls where thin | at least 3 examples per category that carries an action |
| 5 | Test 5 — held-out validation | predicted category matches his on held-out alignments |
| 6 | propose the algorithm | — |
| 7 | build | algorithm validated |

**Running in parallel, no gate:** AnnoSINE2 on three genomes with no SINE
library — snail, starfish, hydra. First prospective test.

---

## 7. Rules

- His calls are ground truth. Curated `tim/` is ground truth. **Human Dfam is
  not** — it contains fragments, an ancestral gene, mixtures and single-copy
  entries.
- Record his words verbatim before doing anything with them.
- Compute on the servers.
- Every claim gets measured against `calls.tsv`, not against a synthetic corpus.
