# SINE / non-SINE discriminator — design spec

Carry-over from a claude.ai planning session. Everything below is design, not code.
Nothing has been implemented yet.

---

## 0. Context

- Author: Sergei Kosushkin, NAU / bioinformatics, transposable element genomics.
- Existing tooling to build on: `sear2k` (SINE genome search, sliding 2k windows),
  `SINEderella` (discovery pipeline), `SubFam` (subfamily consensus), `FaSort10`,
  `SINE_orth_loc`, `LINE_walker`. Standard: MAFFT, ssearch36, bedtools, seqkit, vsearch, EMBOSS cons.
- Existing curated material usable as a labeled corpus: sq2 subfamilies across
  111 squamate species (`sq2_Eryx_all_subfam.fa`, 3,989 seqs), Tal repo
  subfamilies d1–d5 / g1–g7 / t1–t6 across talpids, Squam2 / Anolis consensus work.
- Servers: KIT, Monsoon. KIT runs fail2ban — never issue rapid repeated SSH attempts.

---

## 1. The problem

**Input at the decision point:** a set of genomic loci proposed as ONE family —
output of a `sear2k` search from a pre-existing consensus, or from a tool such as
AnnoSINE. True element boundaries are NOT yet established.

**Required output:** true SINE family vs. not, plus refined boundaries, plus a
cleaned copy list, plus machine-readable reason codes.

Two distinct failure modes need separate verdicts:
1. the set is not a SINE family at all;
2. the set is a mixture — real copies plus contaminants.

## 1.1 The manual criterion being formalized

A SINE family is defined operationally from the shape of a flanked multiple
alignment of ~100 randomly sampled copies. Three tiers, all currently judged by eye
in about a second:

**Tier 1 — alignability.** The copies align to each other at all; a coherent
consensus exists.

**Tier 2 — boundary structure.** The alignment is a three-part step function:
robustly aligned core, flanked left and right by regions where pairwise similarity
collapses to genomic background (the level of arbitrary same-length loci). The
collapse must be *sharp*, not gradual. Two systematic exceptions sit slightly above
background and are expected, not disqualifying:
- an A-rich zone ~10 bp upstream of the core (insertion-site signature);
- residual tail features bleeding into the right flank.

**Tier 3 — internal homogeneity.** Within the core, conservation varies smoothly —
some segments better preserved, others decayed — but is NOT mosaic. No kaleidoscope
where different copy subsets carry different blocks. Copy lengths near-uniform,
within ~10%.

## 1.2 Constraints established during discussion (these are binding)

- **TSDs are one-sided evidence only.** Many SINEs lack them. Presence carries a
  strong positive likelihood ratio (negatives essentially never show them); absence
  carries LR ≈ 1. TSDs may raise a borderline score, never reject a candidate.
  Same treatment for the A-rich zone and tail motifs. Formalize as monotone
  one-sided constraints in the model.
- **Do not flatten the alignment to a mean profile.** A mean identity profile is
  invariant to permuting which copies carry which variants — exactly the
  information Tier 3 depends on. All statistics are computed per copy first, then
  examined as a distribution or a matrix. Aggregation happens only at the end.
- **Nested SINEs must not sink a family.** When a subset of copies sits inside
  another repeat, flank uniqueness drops for that subset only. Diagnose and set
  aside; do not discard the whole set.
- **Consensus position is normally sufficient to set boundaries.** Column
  information content collapses past the element edge because flanks have nothing
  real to agree on. Edge cases exist but are rare and can be neglected.
  Per-copy breakpoints are therefore demoted from boundary-*setting* to
  boundary-*testing*: anchor on the consensus edge, then measure each copy's
  deviation from that anchor. The *spread* of deviations is the discriminative
  quantity — no consensus-level statistic can see it.
  - boundary position ← consensus IC crossing (one robust number)
  - boundary quality ← distribution of per-copy residuals around it

---

## 2. Key methodological ideas

**Per-family permutation null.** Instead of a fixed background-identity threshold
that breaks across genomes and GC contents: permute which left flank is paired with
which copy. Composition preserved exactly, homology destroyed. Score the observed
core-vs-flank contrast against that distribution. Self-calibrating per genome and
per family age.

**Widen flanks to 300 bp.** 50 bp is too short to establish that similarity reaches
background *and stays there*. Costs nothing; stabilizes both the sharpness estimate
and the permutation null.

**Homogeneity as matrix rank (replaces the flattened profile).**
Build M = copy × window identity-to-consensus matrix (rows = individual loci,
columns = sliding windows along the alignment). A true family should be
*separable*: `M[i,j] ≈ a_i · b_j`, where `a_i` = that copy's age and `b_j` =
positional conservation. One factor for how old the copy is, one for which parts of
the element are conserved, nothing copy-specific in between.
Take the SVD; the fraction of variance in the first component is the homogeneity
statistic. Rank-1 → clean family. Rank ≥ 2 → copy subsets have their own positional
structure, i.e. mosaic. Single interpretable number, cheap, preserves per-locus
continuity by construction.

**Mosaicism = block-wise phylogenetic incongruence.** Cross-check the SVD with
off-the-shelf recombination detection (PHI test, GARD, 3SEQ, Neighbor-Net delta) and
with Mantel congruence between sliding-window distance matrices. Note this permits
legitimate subfamily structure — which is *global* and congruent across windows —
while flagging block-swapped composites.

**LINE discrimination falls out of breakpoint spread.** 5′-truncated LINE fragments
give a fixed right edge and a ragged, staircase left edge. Tight-one-side /
broad-other-side is the signature. Bimodal on either side = composite or two merged
families.

**Nesting diagnosis via the copy × copy flank identity matrix** (computed separately
for left and right flanks; never aggregated to a single number):
- true SINE, unique flanks → background everywhere;
- subset inside a host repeat → an elevated **block** among just those copies;
- satellite / segdup → elevated across the whole matrix;
- LINE 3′ end → elevated one side only, graded rather than block-structured.

Practically: identify the elevated block, exclude those copies from the Tier-2
statistic but keep them in the family, evaluate the cliff on the unique-flank
component. Require a sufficient *absolute number* of unique-flank copies (20–30 is
plenty for the permutation null), never a fraction of the total. Trimmed/median
statistics give much of this robustness even before block detection.
Cheap independent pre-screen: count genome-wide hits of each copy's own flank —
high copy number flags nesting without needing a repeat annotation.

**Avoid circularity in sampling.** If copies are drawn from the seed search's own
coordinates, boundaries are inherited from the seed and Tier 2 is partly
self-fulfilling. Re-search with nhmmer at a permissive threshold, extract
coordinates independently, then sample. Sample stratified by bitscore decile —
NOT the top 100 hits, which are the least informative.

---

## 3. Side quest (do this FIRST — it calibrates everything else)

Subsampling is a label-free stability test.

Take 3–5 families known cold from sq2 and Tal with thousands of copies. Draw K = 20
independent subsets at each of n = 25, 50, 100, 200. Run the measurement routine on
all of them. Yields four things at once:

1. **Boundary reproducibility** — do the K consensus-edge calls land in the same
   place? Spread = the boundary rule's real precision, comparable against the
   manual call.
2. **n_min** — whether 100 is right, 50 suffices, or old degenerate families need
   200. n_min almost certainly scales with family age; measure it rather than guess.
3. **Error bars on every feature** — variance of `cliff_z`, `rank1_frac`, `len_cv`
   across subsets of the same family. Thresholds cannot be set without knowing how
   much each statistic wobbles from sampling alone.
4. **A QC runnable on new families with no ground truth** — if two random subsets of
   a candidate disagree on boundaries or consensus, that is evidence of a mixed set,
   independent of any trained model.

Two extra variants:
- **age-stratified subsets** (all young copies vs all old copies) — does the
  boundary rule hold at high divergence?
- **deliberately contaminated subsets** (90% real + 10% random loci) — finds the
  contamination level at which the consensus edge blurs. This directly sets the
  pruning threshold in the refinement loop.

Manual check: dump all K boundary calls as a single overlay against the alignment,
one row per subset, marks at each call — disagreement visible at a glance instead of
requiring K files opened. The existing MSA viewer
(`toki-bio.github.io/MSA-viewer`) already does most of that rendering.

**Caution:** subsets drawn from a family curated by boundary inspection will agree
with the boundary rule partly by construction. Draw from the raw `sear2k` hit set
for the family, not from the curated final copy list.

---

## 4. Pseudocode

### 4.1 Core measurement — identical in training and application

```
function MEASURE(loci, genome, flank=300):

    seqs  = extract(genome, loci, flank_left=flank, flank_right=flank)
    aln   = mafft(seqs)                        # --adjustdirection
    cons  = consensus(aln, threshold=0.5)
    N, L  = shape(aln)

    ##### A. per-copy identity landscape ########################
    # rows = copies, cols = sliding windows. NEVER collapse rows.
    W = 25; STEP = 5
    M = matrix(N, n_windows)
    for i in copies:
        for j in windows:
            M[i,j] = identity(aln[i, window_j], cons[window_j])

    ##### B. boundaries: consensus anchor + per-copy residuals ##
    ic          = information_content_per_column(aln)
    anchorL     = leftmost  column where ic crosses background threshold
    anchorR     = rightmost column where ic crosses background threshold

    for i in copies:
        prof_i      = smooth(M[i, :])
        bg_i        = median(prof_i over outer 100bp of both flanks)
        edgeL[i]    = local crossing of prof_i near anchorL
        edgeR[i]    = local crossing of prof_i near anchorR
        resL[i]     = edgeL[i] - anchorL       # residual, the real variable
        resR[i]     = edgeR[i] - anchorR
        sharpL[i]   = logistic_slope(fit_sigmoid(prof_i near edgeL[i]))
        sharpR[i]   = logistic_slope(fit_sigmoid(prof_i near edgeR[i]))
        implied_len[i] = edgeR[i] - edgeL[i]

    v.resL_spread   = IQR(resL)                # tight = real edge
    v.resR_spread   = IQR(resR)
    v.res_asymmetry = log(resR_spread / resL_spread)   # LINE 5' truncation
    v.resL_modality = dip_test(resL)           # bimodal = composite
    v.resR_modality = dip_test(resR)
    v.sharp_L       = median(sharpL)
    v.sharp_R       = median(sharpR)

    ##### C. per-family permutation null ########################
    for r in 1..500:
        permuted = aln with flanks reassigned to random copies
        null_r   = core_vs_flank_contrast(permuted)
    observed     = core_vs_flank_contrast(aln)
    v.cliff_z    = (observed - mean(null)) / sd(null)
    v.cliff_p    = empirical_p(observed, null)

    ##### D. homogeneity via matrix rank ########################
    # M ~ a_i * b_j  (copy age x positional conservation) => rank 1
    # keep raw values; separability is multiplicative, do NOT mean-centre rows
    Mc          = M restricted to core columns
    S           = svd(Mc).singular_values
    v.rank1_frac = S[1]^2 / sum(S^2)           # high = clean family
    v.rank2_frac = S[2]^2 / sum(S^2)           # high = mosaic

    v.phi_p      = PHI_test(aln[:, core])      # independent cross-check
    v.mantel_min = min over window pairs of
                     mantel(dist_matrix(win_a), dist_matrix(win_b))

    ##### E. flank uniqueness matrix — nesting detection #########
    FL = pairwise_identity_matrix(aln[:, left_flank])
    FR = pairwise_identity_matrix(aln[:, right_flank])
    for F in (FL, FR):
        blocks     = detect_elevated_blocks(F)      # spectral / HAC
        block_frac = fraction of copies inside any block
        global_elev= median(F off-block) - genomic_background
    v.nested_frac_L / v.nested_frac_R       = block_frac
    v.global_flank_elev_L / _R              = global_elev   # segdup/satellite
    v.flank_graded_R = is_graded_not_blocky(FR)             # LINE 3' end

    unique_set = copies in NO elevated block
    if |unique_set| >= 20:
        recompute B and C on unique_set only   # nesting must not sink a family
    v.n_unique = |unique_set|

    ##### F. tandem / satellite check ###########################
    v.neighbor_hit  = frac of copies whose flank matches the BODY of an
                      adjacent-coordinate copy
    v.coord_cluster = clustering index of loci coordinates

    ##### G. length ############################################
    v.len_cv         = sd(ungapped_len) / mean(ungapped_len)
    v.implied_len_cv = sd(implied_len) / mean(implied_len)
    v.big_indel_frac = frac of indels >20bp shared by >10% of copies

    ##### H. one-sided bonuses — may only ADD evidence ##########
    v.tsd_frac    = frac of copies with 4-20bp direct repeat flanking their
                    OWN edgeL/edgeR
    v.arich_score = A+T excess in [anchorL-20, anchorL-5], vs null
    v.tail_score  = poly-A / (TG)n / simple-repeat signal right of anchorR
    # each enters the model with a monotone one-sided (LR >= 1) constraint

    ##### I. sanity ############################################
    v.n_copies     = N
    v.frac_aligned = frac of input seqs retained by aligner
    v.core_identity= median(M over core)       # family age proxy

    return v, aln, anchorL, anchorR, resL, resR, unique_set, M
```

### 4.2 Step 1 — learning

```
TRAIN:
    corpus = []
    for fam in curated_positives:              # sq2, Tal d/g/t, Squam2
        corpus.append( MEASURE(fam.loci), label=SINE )

    for fam in natural_negatives:              # LINE, MITE, satellite, segdup,
        corpus.append( MEASURE(fam.loci), label=fam.class )   # random loci,
                                                              # processed pseudogene
    # synthetic hard negatives — the part that teaches precision
    for fam in curated_positives:
        corpus.append(MEASURE(shift_boundaries(fam, +-20..100)), NOT_SINE)
        corpus.append(MEASURE(splice_halves(fam_a, fam_b)),      NOT_SINE)
        corpus.append(MEASURE(dilute(fam, random_loci, 30%)),    MIXED)
        corpus.append(MEASURE(truncate_5prime(fam, geometric)),  NOT_SINE)

    # simulated families — unlimited labels, known ground truth
    for p in grid(age, copy_number, truncation_rate, subfamily_depth,
                  nested_fraction, host_genome_GC):
        corpus.append( MEASURE(simulate_family(p)), label=p.truth )

    # variable selection
    drop v where |corr(v, v')| > 0.9           # keep the interpretable one
    rank v by mutual_information(v, label)
    report per-variable separation, POSITIVE vs EACH negative class separately
           # must be able to SEE which variable kills which negative

    model = logistic_regression(selected_v)    # or shallow GBM
            with monotone constraints on section-H variables
    calibrate(model, isotonic)
    thresholds = choose(model, target_precision=0.95 on real-only holdout)
    # holdout = real curated data ONLY. never validate on simulations.
```

Notes on the corpus:
- **Processed pseudogenes** (sharp boundaries, poly-A, TSDs) genuinely look
  SINE-like; probably deserve their own class rather than a plain negative label.
- **Segmental duplications / low-copy paralogs** have similarity extending *through*
  the flanks — no boundary at all. Easy but essential negative.
- **Satellites** are cheaply excluded by `neighbor_hit` + coordinate clustering.
- Simulation is credible here because the generative process is well understood:
  plant a consensus into real genomic background with TSDs, vary age,
  substitution/indel model, truncation rate, copy number, subfamily depth.
- Curated positives should be labeled with the manual boundary calls, so
  `resL_spread` / `resR_spread` can be validated against what would be drawn by hand.

### 4.3 Step 2 — application

```
CLASSIFY(consensus, genome):

    loci = sear2k(genome, consensus)           # or AnnoSINE output
    loci = dedupe(loci); loci = drop_tandem_adjacent(loci)

    if |loci| < 20: return INSUFFICIENT

    sample = stratified_sample(loci, n=100, by=bitscore_decile)
             # NOT top-100: top hits are the least informative

    for pass in 1..3:
        v, aln, anchorL, anchorR, resL, resR, unique_set, M =
                MEASURE(sample, genome)
        p = model.predict(v)

        outliers = copies with (residual far from mode)
                            OR (row of M poorly explained by rank-1 fit)
        if |outliers| == 0 or pass == 3: break
        sample = sample - outliers             # prune, realign, rescore

    verdict = SINE      if p > t_hi
              NOT_SINE  if p < t_lo
              REVIEW    otherwise

    reasons = []
    if v.res_asymmetry     > c1: reasons += "ragged 5' edge - LINE-like"
    if v.cliff_z           < c2: reasons += "no boundary above background"
    if v.rank2_frac        > c3: reasons += "block-incongruent - mosaic"
    if v.global_flank_elev > c4: reasons += "flanks non-unique - segdup/satellite"
    if v.nested_frac_L     > c5: reasons += "N copies nested in host repeat"
    if v.len_cv            > c6: reasons += "length heterogeneous"
    if v.tsd_frac          > c7: reasons += "(+) TSDs present"

    return { verdict, p, reasons,
             refined_bounds = (anchorL, anchorR),
             clean_copies   = sample,
             nested_copies  = loci - unique_set,
             diagnostics    = v }
```

---

## 5. Where the real engineering sits

Everything else is standard. These two are not:

1. `detect_elevated_blocks(F)` — block structure in the copy × copy flank identity
   matrix. Needs to distinguish a genuine block from noise at low copy numbers, and
   must not fire on global elevation (segdup/satellite), which is a different verdict.
2. Robust core-run detection in section B when a copy is fragmentary — the
   consensus anchor helps a lot here but does not fully solve it.

## 6. Immediate next steps

1. Run the section-3 subsampling side quest on 3–5 known sq2 / Tal families.
   Verify the consensus-IC boundary rule against manual calls before building
   anything else.
2. Instrument `SINEderella` / `sear2k` to dump `MEASURE()` inputs for every family
   already curated across the squamate and talpid sets — a labeled corpus likely
   already exists in hand.
3. Implement `MEASURE()` alone, with no model, and eyeball its outputs on known
   positives and negatives. Thresholds and model come after the statistics are
   trusted.

## 7. Open items flagged but not resolved

- Points from the original discussion the author marked as unclear and intends to
  return to for clarification.
- Whether `MIXED` should be a third verdict class or handled entirely by the
  pruning loop.
- Genomic background estimation for `global_flank_elev` — per genome, once, or
  per candidate?
