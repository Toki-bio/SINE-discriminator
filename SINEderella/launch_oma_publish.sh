#!/usr/bin/env bash
set -euo pipefail
RUN=/staging/tmp/scorpions/oma/run_oma
cd "$RUN"
rm -rf results/alignments border_loop rebuilt_consensus
export PEEL_FLAGS=/staging/tmp/SINEderella/oma_peel_border_flags.tsv
export PEEL_ALN_DIR=/staging/tmp/scorpions/oma/sd/rebuild
export PATH="/staging/conda/envs/bioinfo/bin:/staging/miniconda3/bin:/usr/bin:$PATH"
exec bash /staging/tmp/SINEderella/run_publish_alignments.sh "$RUN" oma "$RUN/alignments"
