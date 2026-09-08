#!/usr/bin/env bash
# Fast publish: step7 already done; skip slow border loop (use MAFFT_VAL if re-run later).
set -euo pipefail
RUN=/staging/tmp/scorpions/oma/run_oma
export PEEL_FLAGS=/staging/tmp/SINEderella/oma_peel_border_flags.tsv
export PEEL_ALN_DIR=/staging/tmp/scorpions/oma/sd/rebuild
export SKIP_STEP7=1
export SKIP_BORDER_LOOP=1
export PATH="/staging/conda/envs/bioinfo/bin:/staging/miniconda3/bin:/usr/bin:$PATH"
rm -rf "$RUN/results/alignments"
exec bash /staging/tmp/SINEderella/run_publish_alignments.sh "$RUN" oma "$RUN/alignments"
