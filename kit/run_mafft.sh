#!/bin/bash
# Align every corpus set with MAFFT, in parallel. Idempotent: skips finished sets.
set -u
OUT=/data/W/toki/SINE_disc
mkdir -p $OUT/aln_c $OUT/logs

align_one() {
    f="$1"
    b=$(basename "$f" .fa)
    o=/data/W/toki/SINE_disc/aln_c/$b.aln.fa
    if [ -s "$o" ]; then return 0; fi
    mafft --retree 2 --maxiterate 0 --adjustdirection --quiet --thread 1 "$f" > "$o.part" \
        2> /data/W/toki/SINE_disc/logs/$b.err
    if [ -s "$o.part" ]; then mv "$o.part" "$o"; else rm -f "$o.part"; echo "FAIL $b"; fi
}
export -f align_one

ls $OUT/sets_c/*.fa | parallel -j 24 --bar align_one {} 2>&1 | tail -20
echo "ALIGNED: $(ls $OUT/aln_c/*.aln.fa 2>/dev/null | wc -l) / $(ls $OUT/sets_c/*.fa | wc -l)"
echo MAFFT_STAGE_DONE
