#!/usr/bin/env bash
# Align the 4 lanes of paired FASTQ to GRCh38 (no_alt + hs38d1) with minimap2 -x sr,
# then merge, mark duplicates and index. Output: $WORK/bam/$SAMPLE.md.bam
#
# minimap2 (rather than bwa-mem2) because building a bwa-mem2 index needs more
# RAM than this 31 GB machine has; minimap2 sr is fast and adequate for SV/CNV,
# depth and read-backed phasing, which is all the BAM is used for here.
. "$(dirname "$0")/env.sh"

BAM_DIR="$WORK/bam"
mkdir -p "$BAM_DIR"

# Reference: decompress, index, minimap2 index
if [[ ! -s "$REF.fai" ]]; then
  gunzip -c "$DATA/ref/GRCh38_no_alt_hs38d1.fna.gz" > "$REF"
  bio samtools faidx "$REF"
fi
[[ -s "$REF.mmi" ]] || bio minimap2 -x sr -t "$THREADS" -d "$REF.mmi" "$REF"

for lane in L001 L002 L003 L004; do
  out="$BAM_DIR/$lane.bam"
  [[ -s "$out" ]] && { echo "skip $lane"; continue; }
  r1="$DATA/${SAMPLE}_HGWCNDSX7_S16_${lane}_R1_001.fastq.gz"
  r2="${r1/_R1_/_R2_}"
  rg="@RG\tID:HGWCNDSX7.$lane\tSM:$SAMPLE\tLB:$SAMPLE\tPL:ILLUMINA\tPU:HGWCNDSX7.$lane"
  echo "$(date +%T) align $lane"
  bio minimap2 -ax sr -t "$THREADS" -R "$rg" "$REF.mmi" "$r1" "$r2" \
    | bio samtools fixmate -m -u - - \
    | bio samtools sort -@4 -m 1500M -T "$BAM_DIR/tmp_$lane" -o "$out.tmp" -
  mv "$out.tmp" "$out"
done

final="$BAM_DIR/$SAMPLE.md.bam"
if [[ ! -s "$final.bai" ]]; then
  echo "$(date +%T) merge + markdup"
  bio samtools merge -@8 -u -o - "$BAM_DIR"/L00?.bam \
    | bio samtools markdup -@8 -s -f "$BAM_DIR/markdup.stats" - "$final"
  bio samtools index -@8 "$final"
fi
echo "$(date +%T) ALIGN_DONE $final"
