#!/bin/bash
# Re-run ssearch36 per subfamily; save pctid TSV to plots/ for KDE comparison.
set -euo pipefail
RUN_ROOT="${1:-/staging/tmp/scorpions/oma/run_oma}"
S2="$RUN_ROOT/step2/step2_output"
PLOTS="$S2/plots"
CONS="$RUN_ROOT/consensuses.clean.fa"
SUBFAM_DIR="$S2/subfamilies"
THREADS="${THREADS:-8}"
export PATH="/staging/conda/envs/bioinfo/bin:/usr/bin:$PATH"

mkdir -p "$PLOTS"
tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

for sf_fasta in "$SUBFAM_DIR"/*.fasta; do
  [[ -f "$sf_fasta" ]] || continue
  sf=$(basename "$sf_fasta" .fasta)
  out="$PLOTS/${sf}_pctid.tsv"
  [[ -s "$out" ]] && continue
  awk -v sf="$sf" '/^>/{name=$1; sub(/^>/,"",name); found=(name==sf)}
    found{print}' "$CONS" > "$tmpdir/cons.fa"
  [[ -s "$tmpdir/cons.fa" ]] || continue
  awk '/^>/{h=$0; sub(/^>/,"",h); split(h,p,"|"); print ">"p[1]; next} {print}' \
    "$sf_fasta" > "$tmpdir/copies.fa"
  [[ -s "$tmpdir/copies.fa" ]] || continue
  ssearch36 -Q -n -z 11 -E 100 -T "$THREADS" -m 8 \
    "$tmpdir/cons.fa" "$tmpdir/copies.fa" > "$tmpdir/sim.m8" 2>/dev/null || true
  awk -F'\t' '{
    seq=$2; pctid=$3+0; bs=$12+0
    if(!(seq in best) || bs>best[seq]){ best[seq]=bs; pid[seq]=pctid }
  } END { for(s in pid) printf "%s\t%.2f\n", s, pid[s] }' \
    "$tmpdir/sim.m8" > "$out"
  echo "pctid $sf $(wc -l < "$out") lines"
done
echo "Done pctid export to $PLOTS"
