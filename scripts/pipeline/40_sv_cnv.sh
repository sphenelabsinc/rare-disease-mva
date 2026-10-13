#!/usr/bin/env bash
# After alignment: structural variants (Manta), read depth (mosdepth), and the
# inputs for the mosaic-aneuploidy check (per-chromosome depth + B-allele fraction).
. "$(dirname "$0")/env.sh"

BAM="$WORK/bam/$SAMPLE.md.bam"
SV="$WORK/sv"
DEPTH="$WORK/depth"
mkdir -p "$SV" "$DEPTH"

# Primary contigs only, for Manta's call regions
bio bash -c "awk '\$1~/^chr([0-9]+|X|Y)\$/ {print \$1\"\t0\t\"\$2}' '$REF.fai' | bgzip -c > '$SV/primary.bed.gz' && tabix -f -p bed '$SV/primary.bed.gz'"

# Manta (separate env: python2-based)
if [[ ! -s "$SV/manta/results/variants/diploidSV.vcf.gz" ]]; then
  rm -rf "$SV/manta"
  "$MM" run -n manta configManta.py --bam "$BAM" --referenceFasta "$REF" \
    --callRegions "$SV/primary.bed.gz" --runDir "$SV/manta"
  "$MM" run -n manta "$SV/manta/runWorkflow.py" -j "$THREADS" -g 24
fi

# Depth: 1 kb bins genome-wide, plus per-base over the MVA genes for exon-level CNV
bio mosdepth -t 4 -n --fast-mode --by 1000 "$DEPTH/$SAMPLE.1kb" "$BAM"
printf "chr15\t40159022\t40223123\tBUB1B\nchr11\t95787964\t95834693\tCEP57\nchr5\t890883\t921357\tTRIP13\n" > "$DEPTH/mva_genes.bed"
bio mosdepth -t 4 --by "$DEPTH/mva_genes.bed" "$DEPTH/$SAMPLE.mva_genes" "$BAM"
echo SV_CNV_DONE
