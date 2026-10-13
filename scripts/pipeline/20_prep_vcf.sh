#!/usr/bin/env bash
# Prepare the provided VCF: add "chr" prefix (submission format uses chr15 etc.),
# keep primary contigs, split multiallelics, left-align against the reference,
# then QC (bcftools stats, sex check, runs of homozygosity).
# Output: $WORK/vcf/$SAMPLE.norm.vcf.gz (all FILTER values kept; filter downstream)
. "$(dirname "$0")/env.sh"

VCF_DIR="$WORK/vcf"
QC_DIR="$WORK/qc"
mkdir -p "$VCF_DIR" "$QC_DIR"

chrmap="$VCF_DIR/chr_rename.txt"
{ for c in $(seq 1 22) X Y; do echo -e "$c\tchr$c"; done; echo -e "MT\tchrM"; } > "$chrmap"
regions=$(cut -f1 "$chrmap" | paste -sd,)

norm="$VCF_DIR/$SAMPLE.norm.vcf.gz"
bio bcftools view -r "$regions" "$RAW_VCF" -Ou \
  | bio bcftools annotate --rename-chrs "$chrmap" -Ou \
  | bio bcftools norm -m-any -f "$REF" --check-ref w -Ou \
  | bio bcftools view -e 'ALT="*"' -Oz -o "$norm" --write-index=tbi

# QC
bio bcftools stats -f PASS "$norm" > "$QC_DIR/bcftools_stats.txt"
# Sex: heterozygosity on chrX outside PARs and count of chrY calls
bio bcftools view -H -f PASS -r chrX:2781480-155701382 -i 'GT="het"' "$norm" | wc -l > "$QC_DIR/chrX_het.txt"
bio bcftools view -H -f PASS -r chrX:2781480-155701382 -i 'GT="AA"' "$norm" | wc -l > "$QC_DIR/chrX_homalt.txt"
bio bcftools view -H -f PASS -r chrY:2781480-56887902 "$norm" | wc -l > "$QC_DIR/chrY_calls.txt"
# ROH (consanguinity / UPD signal); default 0.4 AF where no population AF available
bio bcftools view -f PASS -Ou "$norm" | bio bcftools roh --AF-dflt 0.4 -G30 -O r -o "$QC_DIR/roh.txt"
echo PREP_DONE "$norm"
