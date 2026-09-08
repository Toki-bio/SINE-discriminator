#!/usr/bin/env bash
# Rebuild row-0 consensus from copy columns + re-justify existing step8 alignments.
set -euo pipefail
RUN="${RUN:-/staging/tmp/scorpions/oma/run_oma}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DISC="${DISC:-$(dirname "$SCRIPT_DIR")}"
export PATH="/staging/conda/envs/bioinfo/bin:/staging/miniconda3/bin:/usr/bin:$PATH"
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
shopt -s nullglob
for dir in "$RUN/results/alignments" "$RUN/alignments"; do
  [[ -d "$dir" ]] || continue
  for f in "$dir"/*.aln.fa; do
    log "rebuild $f"
    python3 "$DISC/rebuild_consensus_row.py" "$f"
    python3 "$DISC/boundary_justify.py" "$f"
    python3 "$DISC/trim_display_flanks.py" "$f" --mode "${TRIM_DISPLAY_MODE:-occupancy}" || true
  done
done
log "Done rebuild_consensus_row + boundary_justify + trim_display_flanks"
