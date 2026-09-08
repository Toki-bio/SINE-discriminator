#!/usr/bin/env bash
# Full publish: step7 already done; border loop + step8a + boundary_justify.
set -euo pipefail
RUN=/staging/tmp/scorpions/oma/run_oma
cd "$RUN"
export PEEL_FLAGS=/staging/tmp/SINEderella/oma_peel_border_flags.tsv
export PEEL_ALN_DIR=/staging/tmp/scorpions/oma/sd/rebuild
export SKIP_STEP7=1
export SKIP_BORDER_SCAN=1
export PATH="/staging/conda/envs/bioinfo/bin:/staging/miniconda3/bin:/usr/bin:$PATH"
export DISC=/staging/tmp/sinedisc
rm -rf results/alignments border_loop rebuilt_consensus
exec bash /staging/tmp/SINEderella/run_publish_alignments.sh "$RUN" oma "$RUN/alignments"
