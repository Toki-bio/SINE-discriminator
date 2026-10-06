# SubFam: what to do with the loci left over after the full chunks of 50

**The problem.** SubFam orders the loci (`mafft --retree 0 --reorder`), cuts chunks of 50, and since the repo version of 2026-09-05 DELETES the last chunk if it has
fewer than 50 loci. Before that it kept it, but `cons -plurality 18` needs 18 agreeing loci, so a short chunk gave a consensus of almost only gaps (junk row).

**Test** (`last_chunk_test.py`, run on KIT, 26 sets x 5-6 replicates): real loci of ccr with the group of every locus known from the 7 peeled groups
(17,550 loci), sampled in 5 compositions (all 7 groups, g4+g6+g7, g2+g3, g5+new13, and ~100-locus runs); each sample ordered as SubFam does, cut by 50,
leftover r = N mod 50 (1-49). Options: A drop (current), B keep short (old, plurality 18), B2 keep short with the plurality scaled to its size,
C last 50 loci (overlapping the previous chunk), D merge into the previous chunk, E balanced split (all chunks 49-50).

**1. Dropping the leftover is biased, not random.** The dropped loci are the tail of the guide-tree order, always the same group:

| composition | dropped loci that are... | share of that group in all loci |
|---|---|---|
| all 7 groups (35 replicates, 870 dropped) | g7 59 %, g2 36 %, g3 0 %, g4 0 % | g7 7 %, g2 20 %, g3 39 % |
| g4+g6+g7 | g7 91 % | 20 % |
| g2+g3 | g2 97 % | 33 % |
| g5+new13 | new13 100 % | 72 % |
| ~100-locus runs (42 replicates) | a whole group disappeared in 0.81 of the replicates | |

**2. Overlapping (C) dilutes the tail.** The last 50 loci are the tail group plus ~25 loci of its neighbour, and ~25 loci are counted twice. Mean purity of the end chunk
(share of its commonest group), leftover r = 11-35, mixtures: keep short with scaled plurality 0.650; E balanced 0.609; C overlap 0.596; D merge 0.530.
Informativeness of the consensus (share of non-gap columns): B2 0.69, E 0.59, C 0.58, D 0.52, B (unscaled plurality 18) 0.35.
For r = 36-49 all options but D are alike (0.57-0.58). For r = 1-10 the short chunk has ~4 loci (purity 0.83 trivially); E and C mix it into a neighbour (0.56-0.59).

**3. Recommendation.** Keep the short chunk and scale the consensus threshold to its size (`plurality = 36 % of its loci`, 18 for 50, never below 2): nothing is lost, nothing
double-counted, the row is a clean consensus of the tail group. If equal chunk sizes are required, use the balanced split instead of the overlap: same purity as the overlap, no
double counting. `SubFam.last_chunk.sh` implements both (`SUBFAM_LAST=keep|balanced|drop`, default `keep`); not installed anywhere: it is a proposal for the SubFam repo.

**Toy test of the script** (KIT, 112 loci): keep -> chunks 50/50/12; balanced -> 38/37/37; drop -> 50/50 (as the current SubFam). Single-locus leftovers are skipped (mafft needs two).
With nothing left over the consensus equals the original's up to MAFFT's own run-to-run variation: **the original SubFam.sh run twice on the same input gives different `.cons`
files** (`--threadit $(nproc)`; the rest of the pipeline uses `--threadit 0` for this reason, SINEderella c34d052).
