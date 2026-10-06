#/bin/sh
# original idea by me, impproved by Nikita, added xargs parallel alignment of chunks
# Variant of SubFam.sh with a choice for the loci left over after the full chunks (SUBFAM_LAST=keep|balanced|drop), see SubFam_last_chunk/README.md.

T0=$(date +%s); T1=$(date +%s)
Secs_HMS() { echo "  done in $(( ${1} / 3600 ))h $(( (${1} / 60) % 60 ))m $(( ${1} % 60 ))s"; }
#
if [[ $# -eq 0 ]] ; then
    echo 'No arguments supplied!'
    exit 1
fi
Files=($(which "mafft") $(which "cons") $(which "seqkit") $(which "seqret") "$1" )
for f in "${Files[@]}" ; do
if [ ! -f "$f" ]; then
echo "File $f: not found"
exit 1
fi
done
if [ -z "$2" ]
        then BnkSz=50
        else BnkSz=$2
fi
# What happens to the loci left over after the full chunks:
#   keep      (default) the short last chunk is kept; its consensus threshold is scaled to its size
#   balanced  all chunks get (nearly) the same size, differing by at most 1: nothing lost, nothing counted twice
#   drop      the short last chunk is deleted (the SubFam.sh of 2026-09-05; it silently removes the tail of the guide-tree order)
LAST=${SUBFAM_LAST:-keep}
NSEQ=$(grep -c '^>' "$1")
#
echo Reordering bank and making consensus sequences
if [ "$LAST" = "balanced" ] && [ "$NSEQ" -gt "$BnkSz" ]; then
        SPLIT="seqkit split2 -p $(( (NSEQ + BnkSz - 1) / BnkSz )) -O ./"
else
        SPLIT="seqkit split2 -s $BnkSz -O ./"
fi
mafft --thread $(nproc) --threadtb $(nproc) --threadit $(nproc) --nuc --quiet --retree 0 --reorder $1 \
        | $SPLIT > /dev/null  2>&1 \
                && for f in stdin.part*.fasta; do mv -f -- "$f" "${f/stdin.part/${1%.*}}"; done \
                && for f in *.fasta; do mv -f -- "$f" "${f%fasta}bnk"; done \
        || exit 1

for b in $(find . -maxdepth 1 -name "${1%.*}_*.bnk" -type f | sort); do
        Count=$(grep -c '^>' "$b")
        if [ "$Count" -lt 2 ] || { [ "$LAST" = "drop" ] && [ "$Count" -lt "$BnkSz" ]; }; then
                echo "Skipping short batch $b ($Count sequences < $BnkSz)"
                rm -f "$b"
        elif [ "$Count" -lt "$BnkSz" ] && [ "$LAST" != "balanced" ]; then
                echo "Keeping short batch $b ($Count sequences); consensus threshold scaled to its size"
        fi
done
#
echo Aligning consensus sequences
find -type f -name "${1%.*}_*.bnk" -print0 |
 xargs -0 -t -I % -P $(nproc) sh -c "mafft --thread $(nproc) --threadtb $(nproc) --threadit $(nproc) --nuc --reorder --quiet '%' > '%'.al"

for b in $(find . -name "${1%.*}_*.al"); do
    N_IN=$(grep -c '^>' "$b")
    PLUR=$(( (N_IN * 36 + 50) / 100 ))          # 36 % of the loci: 18 for 50 loci (as before), 4 for 12, never below 2
    [ "$PLUR" -lt 2 ] && PLUR=2
    cat $b \
    | cons -plurality $PLUR -name $(basename "${b%.*}") -filter \
    | awk '!/^>/{gsub(/[Nn]/, "-")}1' > $(basename "${b%.*}").cons
    rm $b
done;

Secs_HMS "$(($(date +%s) - ${T1}))"; T1=$(date +%s)
#
echo Final alignment
cat ${1%.*}_*.cons > ${1%.*}.clw
mafft --thread $(nproc) --threadtb $(nproc) --threadit $(nproc) --localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --quiet ${1%.*}.clw \
        | seqret -filter -osformat2 msf > ${1%.*}.msf
Secs_HMS "$(($(date +%s) - ${T1}))";
#
echo Job completed
Secs_HMS "$(($(date +%s) - ${T0}))"
