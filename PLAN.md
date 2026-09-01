# SINE discriminator — the plan

Structure: **problem → task → testing → propose the algorithm → build the tool.**

---

## 1. The problem

What was done: I took the 72 alignments you have looked at and commented on in
this project, ran the scoring code on each one, and put its number next to what
you said about that same alignment. Your words are in `calls.tsv`, the numbers
in `calls_scored.tsv`.

**(a) The tool gives the wrong kind of answer.**

When you look at an alignment you name a situation and say what to do about it:
*this is a mixture, split it and re-run each part*; *this is real but the
consensus needs shortening*; *this is real but there are not enough copies to be
sure*.

The tool gives a number from 0 to 100. A number cannot carry that, and the proof
is that four different situations you described all come out as exactly 100.0:

| what you said | n | mean score |
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

**(b) Four alignments where the number itself is simply wrong.**

| what you said | set | score |
|---|---|---|
| MOSAIC | ERI e1-4 | **0.2** |
| SINE_TOO_FEW | AluSx_short_ | **0.1** |
| NOT_SINE | LmeSINE1c | **69.2** |
| SINE (lowest) | — | 68.3 |

---

## 2. The task

Two things, in this order:

1. **Decide what the tool should say.** Not a number — the name of the
   situation, plus what to do about it. The situations are yours, already
   recorded in your own words.
2. **Find which measurements tell those situations apart**, and be honest about
   which ones nothing currently measures.

---

## 3. Testing

There are 72 alignments you have judged, falling into 16 situations you named.
Everything below is checked against what you said, never against sequences I
generated.

### Test 1 — Which measurements tell your situations apart?
Compute **every** variable the project has on all 72 alignments — not just the
four groups that feed the score, but the raw ones: per-position identity and
coverage profiles, A+T, edge sharpness left and right, copy-length spread,
subfamily gap, contamination gap, flank decay where available, TSD fraction,
copy count, consensus-vs-copies agreement.
Then, for each pair of situations, ask which measurements actually differ
between them and by how much.

**Output:** a table of measurement against situation, showing how strongly each
one separates.
**Decision:** which situations can already be told apart, and which cannot.

### Test 2 — What is missing?
For every pair of situations Test 1 cannot tell apart, look at those alignments directly and
name what a person sees that is not being measured. Your own phrases point at
candidates already: *"11 lower sequences"* (a subset of copies, not a
column), *"island of instability"* in a specific region, *"right edge defined at
gggagat"* (a motif at the boundary), *"gappy flanks"* (a display property).

**Output:** a short list of new measurements to try, each tied to the alignment
that demands it.
**Decision:** which are worth building.

### Test 3 — Do the new measurements work?
Build the shortlist, re-run Test 1.

**Output:** an updated table.
**Gate:** a new measurement earns its place only if it tells apart a pair that
could not be told apart before, without spoiling a pair that already worked.

### Test 4 — More judgements, only where we are least sure
The 72 are unevenly spread: 43 are a plain SINE, and nine situations have
exactly one example. One example cannot fix a boundary between situations.

Pick the alignments where the rules are least certain, show those to you, record
what you say. **Ask only where your answer changes something** — never as a
general request to re-judge everything.

### Test 5 — Test it on judgements it has not seen
Set the rules using part of your judgements, then test them on the rest. Report
how often the tool names the same situation you did, and where it fails.

**Gate:** if it only works on the alignments it was tuned on, it is not a tool.

---

## 4. Propose the algorithm

Written only after Tests 1-5, and stating for each situation:

- the measurements that identify it, and in what range
- how much to trust it, given how many examples support it
- **the action it implies** — split, re-extract with longer flanks, shorten the
  consensus, gather more copies, check flanks genome-wide, accept, reject
- what it must NOT be confused with, and what prevents that

Plus, stated openly: the situations that still cannot be told apart, so they are
known limits rather than silent mistakes.

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
| 4 | Test 4 — more judgements where thin | at least 3 examples for every situation that carries an action |
| 5 | Test 5 — test on unseen judgements | the tool names the same situation you do, on alignments it was not tuned on |
| 6 | propose the algorithm | — |
| 7 | build | algorithm validated |

**Running in parallel, no gate:** AnnoSINE2 on three genomes with no SINE
library — snail, starfish, hydra. First prospective test.

---

## 7. Rules

- Your judgements are the standard. The curated `tim/` set is the standard.
  **Human Dfam is not** — it contains fragments, an ancestral gene, mixtures and single-copy
  entries.
- Record your words exactly as written before doing anything with them.
- Compute on the servers.
- Every claim gets checked against your judgements in `calls.tsv`, not against
  sequences I generated.
