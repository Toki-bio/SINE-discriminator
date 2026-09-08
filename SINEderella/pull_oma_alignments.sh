#!/usr/bin/env bash
# Pack oma publish alignments for pull to Windows site repo.
set -euo pipefail
SRC=/staging/tmp/scorpions/oma/run_oma/results/alignments
OUT=/staging/tmp/scorpions/oma/oma_publish_alignments.tgz
tar -czf "$OUT" -C "$SRC" $(ls "$SRC"/oma_*_top100.aln.fa "$SRC"/oma_*_rand100.aln.fa "$SRC"/oma_*_subfam.aln.fa 2>/dev/null | xargs -n1 basename)
ls -lh "$OUT"
base64 -w0 "$OUT" > "${OUT}.b64"
wc -c "${OUT}.b64"
