# Algorithm notes — general findings

What holds across cases, as opposed to what happened on a particular day. The
day-by-day record is `FINDINGS.md`; his method is `METHOD.md`; this file is what
the algorithm should be built on.

---

## 1. A rejection has several independent reasons, and they must be reported separately

`aca_SINE_1` (starfish) is the clearest case so far. Sergei's reading: *"bad and
it includes mosaicism ... not only mosaic, but also weak identity, too long
identity with bad flanks"*. All confirmed by measurement:

| reason | measured | comparison |
|---|---|---|
| weak identity | id_all **0.552** | 0.72–0.83 in the accepted ones |
| **flanks not independent** | flank background **0.450** | 0.27–0.31 in the accepted ones |
| almost no boundary | cliff **0.102** | 0.45–0.53 in the accepted ones |
| mosaic | patch2d **0.143** | 0.061–0.097 in the accepted ones |

The tool scored it **26.3 with no flags at all** — the right number for the
wrong reason, and no explanation. **A verdict without its reasons is not the
answer he needs**, and this is the same defect as four of his properties all
scoring exactly 100.

**Rule: every reason that fires should be reported, not just the one that
happened to move the score most.**

## 2. Flank background is a first-class signal, not a normalisation constant

`aca_SINE_1`'s flank background of 0.450 against a genomic 0.25 says its copies
sit in related contexts — they are not independent insertions. That is one of
his decisive criteria (*"its flanks are good — not alignable"*), and it happened
to be visible only because the cliff computation needs it.

It should be reported in its own right whenever it exceeds background by a clear
margin, with the interpretation attached: satellite, segmental duplication, or
copies inside one larger repeat.

## 3. The mosaic is row-wise WITHIN column, which is why every earlier measure missed it

Neither per-column nor per-copy statistics detect it:

- column jumpiness: mosaic 0.261, legitimate candidate 0.239 — no separation
- block structure (runs vs chance): the *legitimate* set is blockier — backwards
- mean identity: says "less conserved", which is a different statement

What responds is **per-copy variation across regions** (`patch2d`): for each
copy, its concordance with the consensus in each 20 bp window, then how much
that varies *within* that copy. A clean family is uniformly concordant; a mosaic
has copies that are concordant in some regions and discordant in others, with
the regions differing between copies.

Ordered all four new candidates exactly as he judged them (0.061 / 0.097 /
0.143 / 0.150). Four points, so not yet a threshold.

## 4. AnnoSINE copy counts are a floor, not an estimate

On the two snail candidates, AnnoSINE reported 537 and 339 genomic copies. A
full SINEderella run on the same consensuses found **6,016 and 2,848** — 11x and
8x more.

So a low AnnoSINE copy count is not evidence against a candidate, and any
threshold on copy number must be applied after a proper search, never to the
discovery tool's own count.

## 5. Assignment leak transfers across species as a real-family signal

In Timema, the candidates Sergei judged real had leak 0.00–0.18 % while the ones
he judged noisy ran 65–98 %. Both snail candidates come back at **0.00 %**.

Leak does not answer sine/not-sine on its own — he was explicit that subfamily
ambiguity is a different question — but zero leak across thousands of copies is
strong evidence the locus set is coherent.

## 6. Being unrelated to known elements is normal, not suspicious

All-vs-all blastn over SINEBase (230 sequences): **23 % isolated at e ≤ 1e-5,
55 % at e ≤ 1e-20**, median 2 neighbours, maximum 15. SINEBase is small clusters
plus a large isolated tail, not one connected family.

So "no similarity to known SINEs" is a property a real SINE family routinely
has. It must never count against a candidate.

## 7. A SINE–LINE partnership shows at the LINE's 3′ terminus or not at all

Searching a candidate against 1,784 L1 and 7,117 other LINE entries produced
hits of 14–28 bp at e ≥ 0.14, all sitting **1,200–4,600 bp inside** the LINEs.
A genuine partnership puts the SINE's tail on the LINE's last ~50 bp.

**Test the position of the hit, not its existence.** Short internal hits to long
LINEs are background and will always be found.

## 8. Structural evidence is per-copy, and its strength varies

For the two snail candidates:

| | A-box | B-box | TSD |
|---|---|---|---|
| `pom_SINE_0` | 2 mismatches | **0 mismatches** | **86 % of copies** |
| `pom_SINE_1` | 3 mismatches | **0 mismatches** | **35 % of copies** |

A perfect B-box is the hardest part to acquire by chance, so both have an intact
pol-III promoter. But TSD fraction separates them sharply, and it agrees with
every other measure of which is the stronger candidate. **TSD fraction is worth
more than TSD presence.**

## 9. Display faults are not data faults, and the two must not be confused

Flanks that look bad have had three different causes so far, only one of which
was real:

1. **Column width set by one outlier copy** — one copy with a 125 bp flank pads
   seventy others with gaps. Fixed by capping the width to what the copies
   actually use.
2. **A percentile cap is not enough on skewed lengths**, and a constant floor
   cannot help when copies carry ~10 bases. Fixed by choosing the width that
   keeps the panel under 25 % gaps, with the floor following the copies' median.
3. **The copies genuinely have no flanks** — no presentation can fix this, and
   the tool should say so instead of showing a wall of dashes.

A hypothesis that flanks in the aligner distort the element alignment was
**tested and disproved**: re-aligning the element alone changes element gaps by
at most 0.04, in both directions.

## 10. Flank similarity ISLANDS — a signal the score is blind to

Sergei on `aca_SINE_0`: *"aligned left flank has very faint but non-random
islands of similarity"*, and then *"these islands can be far off on the flanks
(both directions)"*.

He is right, and my first reading was wrong: I dismissed the far-upstream ones
as low-coverage artifacts. With a null that controls for how many copies are
present at each column, they are strongly significant — on `aca_SINE_0`,
**49 patches covering 426 columns, up to z = 69**, on both sides and hundreds
of bases out from the element. The coverage control is the whole difference:
without it a column where 12 copies of 100 happen to agree looks like noise
alongside one where 98 do.

Scanned uniformly over all 50 measurable alignments of the prospective corpus
(z > 8 sustained over at least 6 columns, against a null that controls for
coverage):

| alignment | island columns | flank_bg | score |
|---|---|---|---|
| `hyd_SINE_9` | 739 | 0.914 | 0.0 |
| **`aca_SINE_0`** | **426** | **0.271** | **100.0** |
| `aca_SINE_1` (rejected) | 424 | 0.450 | 26.3 |
| `hyd_SINE_17` | 303 | 0.543 | 100.0 |
| `hyd_SINE_0` | 139 | 0.614 | 0.0 |
| everything else | <= 75 | ~0.30 | |

**`aca_SINE_0` scores a perfect 100 and has as much localised flank similarity
as the set he rejected.** The score cannot see it.

### Islands must be measured on the LONG flank, never the display alignment

Same shape as the flank-decay rule, and it bites hard here:

| `aca_SINE_0` measured on | island columns |
|---|---|
| the full 400 bp flank | **426** |
| the trimmed display alignment | **6** |

Trimming cuts the flank panel to the width the copies actually fill. That is
right for judging a boundary and it is what fixed the walls of dashes, but it
throws away exactly the far-out columns Sergei said the islands sit in. Two
flank geometries, two purposes, and they must never be taken from one file.

### Only one of the two high-island sets is invisible to the score

I first wrote that `hyd_SINE_9` was a second case. It is not, and the
distinction is the whole point:

| | island columns | flank_bg | score | seen? |
|---|---|---|---|---|
| `hyd_SINE_9` | 739 | **0.914** | 0.0, SHARED_FLANKS | **yes** |
| `aca_SINE_0` | 426 | **0.271** | 100.0, no flag | **no** |

Taking "invisible" to mean over 150 island columns with a flank average below
0.40, `aca_SINE_0` is the **only** set in the whole 50-alignment corpus.

`hyd_SINE_9`'s flanks are similar *everywhere*, so an average catches it and
the existing rule already rejects it. `aca_SINE_0`'s flanks are ordinary on
average and similar only in patches. **The islands measure earns its place
precisely because it separates these two, which the average cannot.**

### Why the existing flank measure misses it

`flank_bg` is a single average over the whole flank. `aca_SINE_0`'s is 0.271 -
indistinguishable from background - because a few hundred island columns are
diluted by a thousand ordinary ones. **An average cannot detect a localised
signal.** The islands have to be found as islands.

### Validated against his own calls, and it is a CAUTION signal

Run over the whole labelled corpus (673 sets), measuring the **fraction** of
flank inside an island rather than the raw count. Of the 55 sets he has judged:

| his call | n | island fraction |
|---|---|---|
| plain `SINE` | 31 | **<= 0.067** |
| the three flank-caution calls | 3 | **0.19 - 0.37** |
| `NOT_SINE` / `UNUSABLE` / `MOSAIC` | 7 | **>= 0.50** |

Three bands with wide empty gaps. The middle band is exactly the three sets he
asked for flank checks on and nothing else:

- `ERI__eri__e2-4` - *"requires post-processing with proving uniqness of at
  least some flanks on whole-genome level"* (0.195, flank_bg 0.281)
- `ERI__eri__e2-3` - *"with caution and need manual reinspection"* (0.194, 0.371)
- `NEGSEGDUP__eri__r00` - *"unclear situation with left end"* (0.368, 0.326)

**All three have an ordinary flank average.** `NEGSEGDUP__eri__r00` is a
labelled segmental duplication whose flank_bg is 0.326 - the average cannot
tell it from a clean family, and only the islands can.

By class the measure lands on segmental duplications (0.871), LINE ORFs (0.819)
and satellites (0.523), against 0.025 for POS. Every one of the 11 sets with
flank_bg >= 0.40 also has island fraction >= 0.15, so nothing is lost by using
it. It fires on **14 of 273** sets, all in those classes plus hedgehog ERI, and
on **none** of POS, MIXED10, MIXED30, NEGJITTER, NEGSPLICE, NEGTRUNC5, NEGRAND
or NEGCHIM.

**It must be reported, never subtracted.** He still leans SINE on two of the
three calibration sets; penalising would turn `e2-4` (97.4, which he accepts
with a caveat) into a rejection. `verdict.py` emits `FLANK_ISLANDS` at fraction
>= 0.10 and changes no score - verified, 0 of 273 scores moved.

### Where the prospective candidates land

`aca_SINE_0` reads **0.171** at flank_bg 0.271: between his plain-SINE ceiling
and his caution floor, 2.5x the highest plain SINE he judged. It is the only
one of the 52 new-species alignments above 0.10 with an ordinary flank average.

And the signal is in the **top hits only** - `top100` 0.171, `rand100` under
0.03, while `hyd_SINE_9` is 0.331 on both views. So it is not a family sitting
in a shared context; it is a family whose most similar copies do. That points
at a subset inside a larger duplication, and the check is to search those
particular flanks against the starfish genome.

### What it might mean, unresolved

**Still to do:** pull `aca_SINE_0`'s top-hit flanks and search them against the
starfish genome, which is the only thing that will say which of the four
explanations above is right.

## 11. A background measured next to the element is not a background

The single most damaging class of bug in this tool, hit three separate times:
comparing the element against something that contains the element.

`hyd_SINE_0` is the clearest case. Its consensus covers 79 bp of a 208 bp
element, so ~150 bp of SINE sat in what the code called flank. Flank identity
came out **0.614 against an element identity of 0.560** — the background above
the signal — and the tool concluded there was no element, on a family Sergei
reads as a good SINE and RepBase calls SINE2-2B_HM.

**Rule: take the background from far enough out that the element cannot reach
it.** The flank-decay profile already walks outward in 25 bp steps, so its far
end is unrelated DNA wherever the consensus happens to stop. Both the cliff and
the support threshold now use that value, and nothing else in the corpus moved.

Same shape as the earlier failures: the cliff measured over core copies only
(circular), and the flank background measured on element-only alignments
(1.000). Whenever a denominator can contain the numerator, it will.

## 12. Every measurement has its own geometry, and they must not share a file

Three measurements now need three different views of the same loci:

| measurement | needs | breaks if given the other |
|---|---|---|
| boundary by eye | flank trimmed to what copies fill | a wall of dashes |
| flank islands | full 400 bp | 426 island columns become 6 |
| flank decay, microsatellite | full 400 bp | far field never reached |

The trimmed alignment is a **display artefact**. Nothing should ever be measured
on it. This has now caused two separate wrong answers, so it is a rule rather
than an observation.

## 13. Microsatellite content is an axis nothing else in the corpus tests

Sergei asked for it by name on `hyd_SINE_7` (*"combination of microsatellites to
me, not sine - not caught by your filters!"*) and `hyd_SINE_17` (*"surrounded by
long 2-nt microsatellites"*).

Measured as the fraction of ungapped sequence inside a tandem repeat of period
1–6, median over copies, separately for element and flank:

| | element | flank |
|---|---|---|
| `hyd_SINE_7` | **0.517** | 0.000 |
| `hyd_SINE_17` | **0.229** | **0.253** |
| hydra SINEs he accepted | <= 0.092 | 0.000 |
| **every class of the 673-set corpus** | **0.000** | 0.000 |

The last row is the point: POS, NEGSAT, NEGSEGDUP, NEGLINEORF and all the
synthetic classes sit at zero. **No labelled class exercises this property**, so
no threshold could have been fitted from the corpus — the bands exist only
because he judged the hydra sets. A curated negative set only covers the failure
modes someone thought to build.

Two numbers, not one: repeat *in the element* means the consensus collects loci
by repeat content rather than descent; repeat *in the flank* means the copies
sit in tracts where flanks cannot show independence and mapping is unreliable.

## 14. A reason that cannot change the verdict is decoration

Note 1 said every reason that fires should be reported. The inverse turned out to
be just as wrong: three reasons fired and the score ignored them.

- `hyd_SINE_7`: 52 % simple repeat, scored **77**
- `hyd_SINE_6`: 25 of 100 copies in the core, scored **93**
- `hyd_SINE_2`: 8 copies, scored **95**, its own flag text reading *"this score
  is an impression rather than a measurement"*

He rejected all three. A number that reads as acceptance while the text beneath
it says the evidence is not there is worse than no number at all.

**Rule: some reasons cap the score rather than nudging it.** Not a subtraction —
there is no principled amount to subtract — a cap at 45, just under the
acceptance line, saying the question cannot be answered yes on this evidence.
Currently `MICROSATELLITE_ELEMENT`, `SMALL_CORE`, `INSUFFICIENT_COPIES`.

## 15. When the tool can fix the fault, it should fix it

`CONSENSUS_UNDEREXTENDED` told the user to extend the consensus and re-run. It
now does it:

| stage | length | score |
|---|---|---|
| AnnoSINE seed | 110 bp | 0.0 NO_ELEMENT |
| extended by the decay distance | 310 bp | 0.2, flagged **too long** |
| trimmed to the supported window | **211 bp** | **100.0 clean** |

RepBase's SINE2-2B_HM is 208 bp. The loop converges because the two consensus
rules are mirrors: extension is bounded by where similarity reaches background,
and over-extension is caught by the window the copies actually support. Neither
rule needed changing to make the loop work.

This is Sergei's second starting situation — "a library needing adjustment" —
and it is the one place where the tool can do the adjustment rather than hand
back advice.

## 16. Evidence against and evidence missing are different answers

The score answers "how much evidence that this is a SINE family". Three
situations look alike in a single number and must not:

| | example | right answer |
|---|---|---|
| evidence against | 52 % of the element is simple repeat | reject |
| evidence missing | 8 copies, no flanks | **cannot say** |
| evidence for | 93 of 100 copies in a clean core | accept |

Collapsing the middle into a low score put `HUM__hum__AluYb9` - a real Alu
subfamily - at 45.0 because it has few copies. The scorer now carries
`assessable: false` for `INSUFFICIENT_COPIES` and `NO_FLANKS_PRESENT`, and the
page shows it as a third state rather than a rejection. Human rejections went
from 1 to 0.

This is the same rule as "an absent measurement must never be scored as guilt",
which has now been broken four separate times in four different places. It is
worth checking every new reason against it before adding it.

## 17. Unmeasured is not zero, and it is not the maximum either

Three bugs in flank decay, found together when the profile was regenerated over
the whole corpus for the first time, all one mistake:

- a satellite was called on **distance alone**, so a set at 0.255 flank identity
  - background - was called a satellite because its noisy profile never
  confidently dropped
- an unmeasurable distance **defaulted to the 400 bp maximum**, which collapsed
  the uniqueness term to zero
- the distance was **claimed beyond where the profile was readable**: copies
  that carry only 50 bp of flank were credited with 400 bp of similarity

Together these rejected 10 sets including a positive and a set he calls a good
sine, dropping the minimum POS score from 90.1 to 0.2.

**Rule: a measurement may only claim the range it actually covered.** Where the
data runs out, the answer is "not measured", which feeds nothing. And a
classification that depends on a distance must also require the thing being
measured to be present - distance means nothing until there is elevated
similarity to decay from.

## 18. A consensus threshold is a fraction of the sequences, never a count

Building family consensuses from his de novo scan, I reused SubFam's own call,
`cons -plurality 18`. Six of the eleven consensuses then got **zero blastn hits
against the very genome their loci came from**, and none scored above 50.

`-plurality 18` asks for eighteen sequences to agree. Inside SubFam that runs on
chunks of about fifty, where eighteen is a little over a third - a sensible
majority. Applied to a family of twenty it demands near-unanimity, so most
columns come out N; stripping the N's then splices non-adjacent columns into a
chimera that matches nothing.

Rebuilt as a plain column majority - the most common base wherever at least half
the copies have one - the same loci give 47 usable consensuses, all SINE-sized.

**Any parameter borrowed from another stage has to be re-read as a fraction of
whatever it is now being applied to.** The same trap as measuring a background
next to the element: a number that was right in one context, carried unchanged
into another.

## 19. Grouping a de novo locus set: use the provenance, not the clustering

SubFam on the 3,408 loci from the scorpion de novo scan gives 69 chunk
consensuses whose median pairwise identity is **0.205** - near random. That is
not SubFam failing. It chunks by input order after a rough reorder, so it works
when the input is already one family, which is exactly how SINEderella uses it.
A de novo sweep is the opposite case: those loci came from 792 different query
families, so every 50-locus chunk spans many of them.

For a de novo set the grouping is already in the data - **the query that found
each locus**. Grouping by query name gives 47 families of 20+ loci, all
SINE-sized once their consensuses are built properly.

## 20. A discovery tool's candidate list is redundant, and it must be collapsed first

AnnoSINE proposed **160** candidates for zebrafish. All-vs-all blastn gives 424
redundant pairs, and single-link clustering at 80 % identity over 60 % of the
shorter sequence collapses them to **87 families** — one cluster alone holding
**55 near-identical seeds**, with others of 12 and 8, and 83 singletons.

So **46 % of that candidate list is the same families proposed over and over.**
Sergei suspected it on sight: *"maybe zebrafish 160 candidates can be
deduplicated? too much redundancy suspected."*

This matters beyond compute. Every per-class count, every "n of N accepted", and
any threshold fitted on a candidate list is distorted by whichever family
happened to be proposed 55 times. Deduplicate before scoring, and pick the
representative by genomic copy count rather than by order.

It also cuts real time: each candidate costs a blastn against the whole genome,
so on a 1.7 Gb assembly the redundant 73 were about an hour of wasted search.

Related to note 4 — a discovery tool's copy counts are a floor, not an estimate —
and the same lesson from the other side: its candidate *list* is not a family
list either, in both directions.

## 21. A process-matching pattern that can match itself will kill the wrong thing

`pgrep -f "script.py"` matches its own command line, so it reports a dead job as
alive. Already recorded, already fixed with `pgrep -f "script[.]py"`.

The same trap in the other direction is worse: **`pkill -f auto_pipeline.sh`
matched the ssh command that contained that string and killed the session before
the rest of the command ran.** It happened three times tonight - once killing a
zebrafish run mid-flight, twice silently dropping the `sed` that was meant to
follow.

Two rules:

- bracket a character in every `pgrep`/`pkill` pattern: `auto[_]pipeline.sh`
- never put `pkill` and the work it is clearing the way for in the same remote
  command; the kill can take the shell with it


---

## 2026-09-02 — peel on diagnostic columns, not identity

MANUAL §6.1.6 defines a subfamily as a **shared diagnostic pattern** — a
synapomorphy — and says the accumulation of private per-copy mutations is
*noise on top of the subfamily signal*. Average pairwise identity measures
exactly that noise, so the identity-based peel loop was answering the wrong
question. Rebuilt on features (`peel_features.py`).

**A group is a block of co-occurring (position, character) features**, and a
sequence joins if it carries ≥60 % of the block. That tolerance is the point:
one noisy column no longer expels a member, which is what makes the exact
co-occurrence sets in `SINEClusterer` fragment on 600 noisy chunks.

**Seeded growth, not single linkage, over the features.** Linking features
transitively chains them the same way single linkage chained t3/t4/t5 into one
143-chunk cluster — measured: union-find gave 5 blocks and placed 111 of 596
sequences. Taking the widest unused feature as a seed and attaching only what
agrees *with that seed* gave 26 blocks and placed 503.

### Result on Timema, against his t1–t8

| | identity-based | feature-based |
|---|---|---|
| placed | 111 / 596 | **503 / 596** |
| weighted purity | — | **0.958** |
| t3 | never recovered | 105 @ 96 %, 42 @ 98 % |
| groups | 5 of 8 | all of t1, t2, t3, t6, t7, t8 |
| re-align | never fired | fired twice |

Per group: t1 81 @100 %, t2 29 @100 %, t3 105 @96 %, t7 32 @100 %, t8 56 @100 %,
t6 7 @100 %. Residue 93, mostly t1.

**t1 comes out as four separate pure clusters** — 81, 72, 20, 16, all 100 % —
which is his own v4 curation splitting it into `t1_1 t1-2 t1-3 t1-4`. The
substructure is real and the method finds it without being told.

**His re-alignment rule is load-bearing, and this is the first run that showed
it.** Round 2, after 330 sequences left and the remainder was re-aligned, found
a t3 block of **169 co-occurring features** at 98 % purity. In round 1 that block
did not exist. Everything peeled in round 1 under the identity method, so the
effect could never be observed.

### Correction to the ground truth itself

The chunk labels were being built by matching member coordinates against
`assignment_full.tsv` with only `(+)` and `(-)` suffixes — silently dropping the
23.9 % of loci whose strand is `(+,-)`. Adding that suffix raised mean chunk
purity **0.953 → 0.975** and moved chunks between t3 and t4. Earlier t3/t4
conclusions were drawn on labels distorted by that gap.

---

## 2026-09-02 — re-running the page, and two silent-failure bugs

Re-scored the whole corpus with the current code to update the site. Result
first, because it is the boring one: **not a single set changed score.** 273
sets, published vs fresh run, zero differences. The lowest true positive is
81.5 of 132, exactly what the page already claims. The graph *data* was not
obsolete — the prose around it was.

### Bug 1 — the scorer silently degrades to "no evidence"

`verdict.py` loads four side inputs, each in a bare `try/except → {}`:

```python
try:    _DECAY   = json.load(open("flankdecay.json"))
except: _DECAY   = {}
```
…and the same for `flank_islands.json`, `microsat.json`, `region_part.json`.

A missing or *wrong* file is therefore indistinguishable from "this set has no
flank signal", and no-evidence scores **high**. Running the corpus without them
put every negative at 100 and gave a separation of −55 — with no error, no
warning, nothing in the log.

Worse, the files exist under two naming schemes. The corpus versions are
`flankdecay_corpus.json`, `islands_corpus.json`, `microsat_corpus.json`,
`corpus_regionpart.json` (673/673 coverage of `aln_c`). The generic names that
`verdict.py` actually opens hold the **new-species** data — 119 keys, **0 of
which match any corpus set**. Anyone who runs the scorer in that directory
scores the corpus with the flank evidence switched off and is told nothing.

This is the same failure shape as the zero-byte `genome.sizes` that made all 150
zebrafish alignments contain only the consensus row: an input that is absent or
wrong reads as an input that is legitimately empty. **Load these with an
explicit key-overlap check against the corpus and refuse to score at 0 %
overlap.**

### Bug 2 — step7's background is deflated by soft-masking

`step7_boundary_refine.sh` samples 300 random 70 bp genomic regions and compares
them **positionally and case-sensitively**:

```python
ident = (1 - sum(1 for i in range(min_len) if a[i] != b[i]) / min_len) * 100
```

The Timema genome is **48.4 % soft-masked**, so lowercase-vs-uppercase counts as
a mismatch. Measured on `genome.clean.fa`:

| | mean identity | elevated_frac(>45 %) |
|---|---|---|
| case-sensitive (as written) | **13.640 %** | 0.000 % |
| case-insensitive | **27.176 %** | 0.100 % |

13.640 reproduces the `background_identity_pct` in his Timema
`boundary_refinement.tsv` exactly, and scorpion's 14.586 is the same effect.
The true background is ~27 %, near the 0.25 the discriminator assumes.

**Scope of the damage: the reported number, not the calls.** The decision
variable is `elevated_frac`, which is 0.000 case-sensitive against 0.100
case-insensitive, both far under the 5 % threshold — so boundaries are
unaffected. But `background_identity_pct` is wrong in every
`boundary_refinement.tsv` written so far. One `.upper()` fixes it.

### What was actually obsolete on the page

- **"Subfamily structure — unreliable — confounded by age and copy number."**
  Superseded: 503 of 596 Timema chunks at 0.958 purity on diagnostic columns.
- **"Natural negatives — missing — no LINEs, satellites or processed
  pseudogenes yet."** Contradicted by the page's own data, which already
  contains **44** of them: 14 LINE, 16 chimera, 8 mosaic, 3 satellite,
  3 segmental duplication.
- **No subfamily section at all**, despite that being the strongest result.

Added a Subfamilies section, fixed both rows, added the nav entry. Verified by
rendering the page in headless Chrome: no console errors, table fits, nav has
the new entry.

---

## 2026-09-02 — the alignments were never built with `sear`

`grep -rn "sear " *.py *.sh` over the whole repo returns **nothing**. Every
alignment in this project was built by `candidate_to_aln.py`, which is my
reimplementation of his search:

```python
FLANK = 400
blastn -query q -subject genome -evalue 1e-10 -word_size 11 -dust no
       -outfmt '6 sseqid sstart send pident length bitscore sstrand'
```

then a hand-written single-pass overlap dedupe and `bedtools slop`.

**`sear -b100` already does this, better, and was sitting there the whole time:**
top N by bitscore → ±50 bp flanks → `getfasta -s` → query prepended as row 1 →
`mafft --localpair --maxiterate 1000 --ep 0.123`. Correctly oriented by
construction. What differs:

| | his `sear` | my `candidate_to_aln.py` |
|---|---|---|
| engine | `ssearch36`, Smith-Waterman | `blastn -word_size 11` |
| hit filter | span ≥80 % of query length, ≥65 % identity | `evalue 1e-10`, my own `hi-lo >= 40` |
| merge | `bedtools merge -s` | hand-written single-pass overlap loop |
| flank | 50 | **400** |
| top-100 alignment | built in, `-b100` | rebuilt by hand |
| orientation | stranded extraction, consensus row 1 | `--adjustdirection`, lost |

That single `FLANK = 400` is upstream of everything in this conversation:
`conse` on 400 bp-flanked copies calls **541 bases for a 253 bp element**, and
the over-extension then propagates through
`measure_c.py:91` / `justify_all.py:24`, which take the element window as the
consensus's first-to-last non-gap column with no validation.

**Sixth reimplementation.** SubFam, conse, sine_consensus.sh,
step7_boundary_refine.sh, step8a, and now sear.

### His calibration call on NEGTRUNC5__saq__s5_5seqs

He read the alignment by eye and said the consensus should start at **AGTTCGA**.
Measured, copy occupancy through the consensus:

| consensus base | occupancy |
|---|---|
| 1 | 0.34 |
| 25 | 0.51 |
| 43 | 0.56 |
| 61 | 0.70 |
| **64 = AGTTCGA** | **0.77** |
| 70–88 | 0.82 / 0.74 / 0.77 / 0.79 |

Occupancy rises monotonically and **plateaus at AGTTCGA**. My `OCC >= 0.50`
threshold fires at base 25 — 39 bases too early. **The boundary he sees is where
the occupancy curve flattens, not where it crosses a fixed fraction.**
`GTTCGA` is the tRNA B-box core, so these copies are 5′-truncated back to the
B box, which is exactly what the NEGTRUNC5 class is constructed to be.

His second observation on the same file — "flanks are incompletely sampled in
many copies and the left end of the flank is hardly identifiable" — is the same
fact from the other side: below 0.5 occupancy half the copies have no base
there, so there is no flank edge to see.

### Flank sizes

**50 L / 70 R**, his instruction. (`extract_alignments.sh` uses 30/70;
`step8a_extract_alignments.sh` uses `BASE_UP=50`, `BASE_DOWN=70`. The step8a
values are the ones to use.)

### What this means for the rebuild

Not a patch to the alignments. The product has to be rebuilt from the search up:

1. `sear -b100 <consensus> <genome>` for the top-100 view, `sear` + sampling for
   rand100 — his tool, his filters, his orientation handling.
2. Flanks 50 L / 70 R.
3. Keep a separate 400 bp extraction **only** for decay/island/microsatellite
   measurement, which is what `run_species.py` always said it was for — and
   never publish those files to the viewer.
4. Boundary from the occupancy plateau, calibrated against his AGTTCGA call.
5. Block on consensus length vs median copy length, so a 541-base consensus for
   a 253 bp element cannot ship.
