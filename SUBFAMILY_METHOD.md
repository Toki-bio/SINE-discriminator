# Separating SINE subfamilies — current state

The task: take SINEderella's SubFam output and recover the subfamily partition a
human would make by eye, with a confidence figure attached. This is the record of
what works, what does not, and what his judgements have calibrated so far.

Written 2026-09-03. Supersedes the scattered notes in `ALGORITHM_NOTES.md`, which
remains the dated log.

---

## 1. What a subfamily is

MANUAL §6.1.6, his definition:

> A *subfamily* is a lineage within a family, defined by a **specific, shared,
> diagnostic pattern** of small indels/SNPs common to that lineage's copies (a
> synapomorphy) — **not** by the generic accumulation of private per-copy
> mutations from ordinary post-insertion decay, which is noise on top of the
> subfamily signal rather than the signal itself.

Two consequences that took a long time to act on:

- **Average pairwise identity is the wrong axis.** It measures the decay the
  definition calls noise.
- **Indels come first.** A rule that discards gaps cannot express the definition.

---

## 2. The pipeline that works

Input is SubFam's chunk consensuses (`input.clw`, degapped and realigned per
MANUAL §6.1.1) or, for a grey-zone call, individual copies from the 50-member
`.bnk` banks.

```
1. trim edges          sliding 20-col window, gap fraction > 0.6 → trim
                       (MSA-viewer getTrimBoundaries, his "optimal" preset)
2. variable sites      keep columns where the minority count >= 15 % of nSeq,
                       counting only sequences whose own span covers the column,
                       and counting internal gaps as a character
                       (MSA-viewer script.js:5555)
3. features            every (position, character) INCLUDING gap, present in
                       >= MIN_SET and <= 50 % of the pool; skip columns where one
                       BASE dominates > 80 % of non-gap
4. blocks              seeded growth over features by Jaccard >= 0.45 — take the
                       widest unused feature as seed, attach what agrees with the
                       seed, emit, repeat. NOT transitive linkage
5. groups              sequences carrying >= 60 % of a block's features
6. peel                remove the group, re-align the remainder when >= 100 have
                       gone (>= 40 in later rounds), repeat
7. score               ordering-adjacency + separator distribution (§4)
```

Every step 1–2 and the linkage choice in 4 came from reading his tools, not from
invention.

## 3. What each fix was worth, measured

Timema v4, 595 chunk consensuses, 13 curated subfamilies. Weighted mean purity of
the peeled groups against his partition:

| change | placed | purity | note |
|---|---|---|---|
| identity-based clustering | 111 / 596 | — | t3 never recovered |
| features, gaps skipped | 493 / 595 | 0.826 | t3-1/t3-2 collapse into one cluster at 58 % |
| **+ gap as a feature state** | **584 / 595** | **0.926** | t3-1 63 @ 100 %, t3-2 53 @ 94 % |
| **+ edge trim before clustering** | 583 / 595 | **0.937** | t1-4 24 @ 67 % → 34 @ 91 % |
| best swept setting (`MIN_SET 8`) | 591 / 595 | 0.941 | residue 4 |

The ceiling is 0.882 — the purity of the chunk labelling itself, since each chunk
is labelled by majority vote of its 50 members.

---

## 4. Scoring a proposed split

Two independent axes. Neither alone is enough.

### 4a. Ordering-adjacency — *how consistently do the groups sort apart*

Build the average-linkage tree on MAFFT's own 6-mer distance, take the **leaf
order**, and measure the fraction of adjacent pairs sharing a label. Normalise:

```
score = (observed − chance) / (max − chance)      max = (n−2)/(n−1)
chance = sum over groups of (size/n)^2
```

Computed on **variable sites at 15 %**, which sharpens it considerably.

Two earlier attempts that FAILED, recorded so they are not retried:

| test | t3 (he named it instantly) | verdict |
|---|---|---|
| mean 6-mer distance, between − within | **−0.105** | says unseparable; fails because t3-2's internal spread (0.994) exceeds its distance to t3-1 (0.888) |
| top-level guide-tree cut at k=2 | splits 119 vs 2 | peels outliers; the highest split is not where the structure is |
| **leaf-order adjacency** | **0.90** | works |

A group mean asks "are these clouds apart" and breaks on internal heterogeneity.
A cut asks "is the top split here". Adjacency only looks at neighbours, and
survives both.

### 4b. Separator distribution — *is the evidence independent*

Collapse consecutive variable sites into **events**, and count them **in a
reference sequence's coordinates, not in alignment columns** — gap columns split
one event into several and inflate the count (t1-4 reads 78 in column space, 31
in reference space).

Contiguous variable columns are **one mutational event**, not many. t3's "14 gap
columns" at 207–218 are one indel.

| at 15 %, reference coords | reference | variable cols | events | sites/event | median separator | largest event |
|---|---|---|---|---|---|---|
| **t3-1/t3-2** | 282 bp | 25 (**9 %**) | 13 | 1.9 | **19 bp** | 1 bp |
| **t1-4/t1-2** | 263 bp | 213 (**81 %**) | 31 | 6.9 | **1 bp** | **51 bp** |

**Events are not alignment noise** — realigned five ways (localpair ±`--ep
0.123`, `--auto`, `--retree 2`, `--globalpair`), 100 % of events in both pairs
survive every aligner. My hypothesis that packed events were artefacts was tested
and rejected.

The real reading:

- **few, strong, dispersed, consistent** → a real subfamily split (t3: 9 % of the
  element variable, columns up to |diff| 0.95, median 19 bp apart)
- **many, weak, packed, inconsistent** → one heterogeneous group (t1-4/t1-2:
  **81 %** of the element variable and still nothing above |diff| 0.66 — two pools
  each internally diverse and mutually overlapping)

**Fraction-variable was tested across pairs and does NOT work** — it measures
evolutionary distance, not separability. The clearly-distinct t2/t8-2 pair has the
highest value (0.70) of the three judged pairs. The discriminator is **`best`**,
the strongest single consistent column: 0.96 and 1.00 for splits, **0.66** for the
pair he merges. His instruction to reject anything below the t1-4 case sets the
bar at **>= 0.85 split, <= 0.70 merge**. See §4c.


### 4c. `best` — *is any of the difference consistent* (the primary test)

The strongest single column, max over |gap-fraction difference| and
occupancy-guarded |majority-base difference|. Measured on **copies**, never chunk
consensuses.

| pair | his call | frac_var | **best** |
|---|---|---|---|
| t3-1 / t3-2 | split, named at a glance | 0.08 | **0.96** |
| t1-4 / t1-2 | merge | 0.39 | **0.66** |
| t2 / t8-2 | distant, distinct | 0.70 | **1.00** |

**Threshold: >= 0.85 split, <= 0.70 merge**, from his instruction to reject
anything below the t1-4 case.

**Chunk consensuses give the wrong answer here.** At chunk level t1-2/t1-4 reads
`frac_var 0.13, best 0.86` — a clean split — because `cons -plurality 18` over 50
copies smooths each chunk. Across all 45 chunk-level pairs `best` is 0.97-1.00
almost everywhere: at chunk level nearly everything looks separable. Chunk level
proposes candidates; copy level judges them.

---

## 5. Calibration against his judgements

| pair | his words | adjacency (15 % var sites) | median separator |
|---|---|---|---|
| **t3-1 / t3-2** | named the insertion at a glance | **0.90** | 19 bp, 9 % variable |
| **t1-4 / t1-2** | "faint difference in 3' part by small deletion and few other places. still very borderline" | **0.55** | 1 bp, 81 % variable |

Provisional bands, two points only:

- **≥ 0.85 with dispersed evidence** — a separation he names without hesitating
- **~0.55, evidence packed** — grey zone; report as such, do not split
- **≈ 0** — no structure in the ordering

**Missing: an anchor at the bottom.** No case he has rejected outright.

### His curation calls recorded

- **t3 / t4 / t5** are one unit in his own hands — `tim/subfam/t345.al`. My
  repeated "cannot separate t3/t4/t5" was measuring something he had already
  decided.
- **t1-2 + t1-4** — "more plausible to merge". Supported: nothing separates them
  above |diff| 0.66, strays leak into t1-2, t1-4 fragments at all 26 swept
  settings. Merged they are 141,884 of 339,090 copies (42 %).
- The **TACAT indel** is at 3′ alignment columns 441–447: deleted in **59 %** of
  t1-4 copies against 16 % of t1-2. Real, partially penetrant, not a
  synapomorphy.

---

## 6. Consensus quality is a confounder

The curated `t1-4` consensus is defective in two measurable ways:

- **12 bp short** — 242 against a rebuilt 253 and t1-2's 254
- it asserts **TACAT at body columns 315–320 where only 24 % of its own copies
  agree** (21 read TACAT; 7 TACAC, 7 TATAT, 6 TAAAT)

Rebuilt with `conse` at 35 % over the top 100 by bitscore, once: identity to t1-2
falls from 0.930 to 0.897. Part of that may be the longer consensus reaching into
more variable sequence rather than real divergence — not offered as settled.

**A grey-zone call cannot be trusted while a consensus is wrong.** Rebuild first.

---

## 7. Method rules learned from him

- **Use the general subfam alignment as the source.** Extracting one group and
  re-aligning it alone destroys the comparison that makes a separation visible.
- **Drop to the 50-member banks for a grey-zone call.** `cons -plurality 18` is a
  majority rule over 50 copies, so a diagnostic carried by a minority within each
  chunk cannot survive it. Any weak or partially penetrant difference is
  invisible at chunk level regardless of parameters.
- **Edges are used only when they carry firm evidence well above the noise level
  characteristic of simple repeats** — otherwise trim or downweight. There is a
  mechanical reason: k-mers spanning a ragged edge are mostly unique, so they
  enter a distance's denominator but never its shared count, inflating distance
  systematically.
- **Never report a single-setting result as a property of the data.** His
  `analyzeClusterability` exists for this: "a single run cannot distinguish
  settings that are too strict from data with no structure to find". Every
  subfamily I called unseparable came from one setting; t3 and t6-1 both
  dissolved when swept.
- **The nearest competitor falls out of the peel for free** — the group its
  strays leak into.


### 7a. Never compare SubFam chunk names across runs

`run_subfam_per_sf.sh` runs SubFam separately per subfamily, and each run numbers
its chunks `input_001…` from scratch. Nine `.al` files all containing
`input_001.bnk` hold nine **different** chunk sets sharing names. I compared them
by name, got 100 % overlap, and reported that the files were all the same 200
chunks — measuring a naming collision.

Compare by sequence, or by the genomic coordinates of the members in the `.bnk`.

---

## 8. Open

- **t1-4 is the one genuine remaining failure**: 88 chunks, fragments at every one
  of 26 settings, best 41 @ 88 %. His call is to merge it into t1-2.
- ~~Are t1-4's packed events alignment artefacts?~~ **Tested and rejected.**
  Realigned five ways, 100 % of events survive in both pairs. In reference
  coordinates t1-4 has 31 events (not the 78 I reported from alignment-column
  space) covering **81 % of the element**, against t3's 25 columns over **9 %**.
  The rule that falls out: *few, strong, dispersed, consistent* = a real split;
  *many, weak, packed, inconsistent* = one heterogeneous group. This is the
  strongest support yet for merging t1-2 and t1-4.
- **`EXCL_MIN` is inert** — four values give byte-identical output. Remove it.
- **The peel reports a partition, not a confidence.** It should emit the
  adjacency score and separator statistics per group, so a grey-zone group is
  labelled rather than silently shattered into eight clusters that look like
  eight findings.
- **saq s1–s9** — recreated from scratch: 9 groups, 598 of 600 placed, three of
  his consensuses recovered at identity 1.000. Over-splits s7g (3 pieces) and s8
  (2), loses the three smallest (5, 7, 20 chunks). His original chunk→group
  labels could not be found; the KIT shell history shows the curation was done
  visually in MSA-viewer, not on the server, so it would be in browser
  localStorage or a local download. **`asSINEment` cannot reconstruct it** —
  it assigns only 65 of 600, because 10/10 unanimity almost never holds between
  consensuses that are 96–99 % identical to one another.
- Bottom-of-scale calibration anchor still missing.

---

## 9. Tools

| file | what |
|---|---|
| `peel_features.py` | the peel; env-overridable parameters; trim + gap-state features |
| `sweep_peel.py` | parameter sweep, after `analyzeClusterability` |
| `mafft_dist.py` | MAFFT's 6-mer distance, ported from `disttbfast.c` / `mltaln9.c` and verified |
| `tree_cut.py` | UPGMA cut with cut-height and next-merge height |
| `sweep_table.txt` | the 26-run sweep |

Do not add to these without first checking `REINVENTED.md`.
