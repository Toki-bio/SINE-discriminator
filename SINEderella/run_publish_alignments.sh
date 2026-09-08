#!/usr/bin/env bash
# Default SINEderella publish path for Tal-style report alignments.
#
# Replaces bare extract_top100_rand100_subfam.sh for all runs where flanks
# may not be unique at 50+70 bp (i.e. all production publish paths).
#
# Chain (existing tools, wired):
#   1. step7_boundary_refine.sh  — widen flanks to background or 1000 bp cap
#   2. needs_border_loop.py      — list subfamilies needing element border loop
#   3. border_loop_subfam.py     — outward BED moves, rebuild consensus (5′ then 3′)
#   4. step8a_extract_alignments.sh — extract with confirmed boundary extensions
#   5. rebuild_consensus_row.py  — row 0 = copy majority at CONS_EDGE (0.50)
#   6. boundary_justify.py       — copy-supported element window + display justify
#   7. trim_display_flanks.py      — drop abandoned end columns (occupancy or hybrid)
#
# Usage:
#   run_publish_alignments.sh <RUN_ROOT> <SPECIES_CODE> [OUT_DIR]
#
# Environment:
#   SINEDERELLA_BIN  — directory with step7/step8a (default /staging/tmp/SINEderella)
#   PEEL_FLAGS       — optional TSV of subfamilies flagged at peel step4
#   PEEL_ALN_DIR     — optional dir of peel loci alignments ({sf}.aln)
#   DISC             — Python tools dir (default: repo root, parent of SINEderella/)
#   SKIP_BORDER_SCAN — set to 1 to skip needs_border_loop --scan
#   TRIM_DISPLAY_MODE — occupancy (oma publish) or hybrid (general corpus)
set -euo pipefail

RUN_ROOT="${1:?usage: $0 <RUN_ROOT> <SPECIES_CODE> [OUT_DIR]}"
SPECIES="${2:?usage: $0 <RUN_ROOT> <SPECIES_CODE> [OUT_DIR]}"
OUT_DIR="${3:-$RUN_ROOT/alignments}"

SINEDERELLA_BIN="${SINEDERELLA_BIN:-/staging/tmp/SINEderella}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DISC="${DISC:-$(dirname "$SCRIPT_DIR")}"
PEEL_FLAGS="${PEEL_FLAGS:-}"
PEEL_ALN_DIR="${PEEL_ALN_DIR:-}"

export PATH="/staging/conda/envs/bioinfo/bin:/staging/miniconda3/bin:/usr/bin:$PATH"

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }

if [[ "${SKIP_STEP7:-0}" != "1" ]]; then
log "step7: boundary refinement (50 bp steps, 1000 bp cap)"
"$SINEDERELLA_BIN/step7_boundary_refine.sh" "$RUN_ROOT" 50 50 1000
else
log "step7: skipped (SKIP_STEP7=1)"
fi

if [[ "${SKIP_BORDER_LOOP:-0}" == "1" ]]; then
  log "border loop: skipped (SKIP_BORDER_LOOP=1)"
  cp -f "$RUN_ROOT/step2/step2_output/assigned.fasta" \
    "$RUN_ROOT/step2/step2_output/assigned.publish.fasta"
  cp -f "$RUN_ROOT/consensuses.clean.fa" "$RUN_ROOT/consensuses.publish.fa"
else
log "detecting subfamilies needing border loop"
NEED_ARGS=(python3 "$SCRIPT_DIR/needs_border_loop.py" "$RUN_ROOT")
[[ -n "$PEEL_FLAGS" && -f "$PEEL_FLAGS" ]] && NEED_ARGS+=(--peel-flags "$PEEL_FLAGS")
[[ "${SKIP_BORDER_SCAN:-0}" != "1" ]] && NEED_ARGS+=(--scan)
mapfile -t BORDER_SFS < <("${NEED_ARGS[@]}" || true)

if ((${#BORDER_SFS[@]})); then
  log "border loop: ${#BORDER_SFS[@]} subfamilies: ${BORDER_SFS[*]}"
  for sf in "${BORDER_SFS[@]}"; do
    PEEL_ALN=""
    if [[ -n "$PEEL_ALN_DIR" ]]; then
      if [[ -f "$PEEL_ALN_DIR/${sf}.aln" ]]; then
        PEEL_ALN="$PEEL_ALN_DIR/${sf}.aln"
      elif [[ -f "$PEEL_ALN_DIR/${sf#oma_}.aln" ]]; then
        PEEL_ALN="$PEEL_ALN_DIR/${sf#oma_}.aln"
      fi
    fi
    python3 "$SCRIPT_DIR/border_loop_subfam.py" "$RUN_ROOT" "$sf" \
      ${PEEL_ALN:+--peel-aln "$PEEL_ALN"} \
      --sine-script "$SINEDERELLA_BIN/sine_consensus.sh" || {
        log "WARN border loop failed for $sf — continuing with seed consensus"; }
  done
  python3 "$SCRIPT_DIR/apply_border_to_assigned.py" "$RUN_ROOT"
  python3 "$SCRIPT_DIR/merge_rebuilt_consensus.py" "$RUN_ROOT"
else
  log "border loop: none flagged"
  cp -f "$RUN_ROOT/step2/step2_output/assigned.fasta" \
    "$RUN_ROOT/step2/step2_output/assigned.publish.fasta"
  cp -f "$RUN_ROOT/consensuses.clean.fa" "$RUN_ROOT/consensuses.publish.fa"
fi
fi

for need in step7_boundary_refine.sh step8a_extract_alignments.sh; do
  [[ -x "$SINEDERELLA_BIN/$need" ]] || {
    echo "ERROR: missing $SINEDERELLA_BIN/$need" >&2; exit 1; }
done
[[ -f "$DISC/boundary_justify.py" ]] || {
  echo "ERROR: missing $DISC/boundary_justify.py" >&2; exit 1; }
[[ -f "$DISC/rebuild_consensus_row.py" ]] || {
  echo "ERROR: missing $DISC/rebuild_consensus_row.py" >&2; exit 1; }

log "step8a: extract alignments (boundary_refinement.tsv + publish loci/consensus)"
STEP8_OUT="$RUN_ROOT/results/alignments"
mkdir -p "$STEP8_OUT"
# step8a reads consensuses.clean.fa and assigned.fasta — swap in publish versions
cp -f "$RUN_ROOT/consensuses.clean.fa" "$RUN_ROOT/consensuses.clean.fa.pre_publish.bak"
cp -f "$RUN_ROOT/step2/step2_output/assigned.fasta" \
  "$RUN_ROOT/step2/step2_output/assigned.fasta.pre_publish.bak"
cp -f "$RUN_ROOT/consensuses.publish.fa" "$RUN_ROOT/consensuses.clean.fa"
cp -f "$RUN_ROOT/step2/step2_output/assigned.publish.fasta" \
  "$RUN_ROOT/step2/step2_output/assigned.fasta"
"$SINEDERELLA_BIN/step8a_extract_alignments.sh" "$RUN_ROOT" "$SPECIES"
cp -f "$RUN_ROOT/consensuses.clean.fa.pre_publish.bak" "$RUN_ROOT/consensuses.clean.fa"
cp -f "$RUN_ROOT/step2/step2_output/assigned.fasta.pre_publish.bak" \
  "$RUN_ROOT/step2/step2_output/assigned.fasta"

# Prefer step8a output; fall back to OUT_DIR symlink/copy
if [[ "$OUT_DIR" != "$STEP8_OUT" ]]; then
  mkdir -p "$OUT_DIR"
  cp -a "$STEP8_OUT"/*.aln.fa "$OUT_DIR/" 2>/dev/null || true
fi

log "rebuild_consensus_row: row 0 from copy columns in each MSA"
shopt -s nullglob
for f in "$STEP8_OUT"/*.aln.fa; do
  python3 "$DISC/rebuild_consensus_row.py" "$f"
done
if [[ "$OUT_DIR" != "$STEP8_OUT" ]]; then
  for f in "$OUT_DIR"/*.aln.fa; do
    [[ -f "$f" ]] || continue
    python3 "$DISC/rebuild_consensus_row.py" "$f"
  done
fi

log "boundary_justify on published alignments"
shopt -s nullglob
for f in "$STEP8_OUT"/*.aln.fa; do
  python3 "$DISC/boundary_justify.py" "$f"
done
if [[ "$OUT_DIR" != "$STEP8_OUT" ]]; then
  for f in "$OUT_DIR"/*.aln.fa; do
    [[ -f "$f" ]] || continue
    python3 "$DISC/boundary_justify.py" "$f"
  done
fi

TRIM_MODE="${TRIM_DISPLAY_MODE:-occupancy}"
[[ -f "$DISC/trim_display_flanks.py" ]] || {
  echo "ERROR: missing $DISC/trim_display_flanks.py" >&2; exit 1; }
log "trim_display_flanks (mode=$TRIM_MODE)"
for f in "$STEP8_OUT"/*.aln.fa; do
  python3 "$DISC/trim_display_flanks.py" "$f" --mode "$TRIM_MODE" || true
done
if [[ "$OUT_DIR" != "$STEP8_OUT" ]]; then
  for f in "$OUT_DIR"/*.aln.fa; do
    [[ -f "$f" ]] || continue
    python3 "$DISC/trim_display_flanks.py" "$f" --mode "$TRIM_MODE" || true
  done
fi

log "Done. Alignments in $STEP8_OUT (and $OUT_DIR if different)"
