#!/bin/bash
# Publish oma report assets: all-consensi alignment + step1 SubFam input (30k sample).
set -euo pipefail
RUN_ROOT="${1:-/staging/tmp/scorpions/oma/run_oma}"
OUT_DIR="${2:-$RUN_ROOT/alignments}"
SITE_ALIGN="${3:-}"   # optional local site alignments dir to copy into

CONS="$RUN_ROOT/consensuses.clean.fa"
SUBIN="$RUN_ROOT/genome.clean_step1/subfam_input/input.clw.al"
mkdir -p "$OUT_DIR"

echo "MAFFT all consensuses -> oma_consensuses.aln.fa"
mafft --auto --quiet "$CONS" > "$OUT_DIR/oma_consensuses.aln.fa"

echo "Step1 SubFam input (30k loci sample) -> oma_subfam_input.aln.fa"
cp -f "$SUBIN" "$OUT_DIR/oma_subfam_input.aln.fa"

n_cons=$(grep -c '^>' "$CONS" || true)
n_sub=$(grep -c '^>' "$SUBIN" || true)
echo "done: $n_cons consensuses, $n_sub subfam_input rows"
ls -la "$OUT_DIR/oma_consensuses.aln.fa" "$OUT_DIR/oma_subfam_input.aln.fa"

if [[ -n "$SITE_ALIGN" && -d "$SITE_ALIGN" ]]; then
  cp -f "$OUT_DIR/oma_consensuses.aln.fa" "$SITE_ALIGN/"
  cp -f "$OUT_DIR/oma_subfam_input.aln.fa" "$SITE_ALIGN/"
  cp -f "$RUN_ROOT/consensuses.clean.fa" "$SITE_ALIGN/oma_consensuses_v2.fa" 2>/dev/null || true
fi
