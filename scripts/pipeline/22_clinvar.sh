#!/usr/bin/env bash
# Annotate the normalised VCF with ClinVar (GRCh38 weekly release 2026-10-04).
# ClinVar uses 1..22/X/Y/MT contig names, so it is renamed to the chr style of
# 20_prep_vcf.sh before annotation. Adds INFO/CLNSIG, CLNREVSTAT, CLNDN and puts the ClinVar
# variation ID in the ID column; records without a ClinVar match are left unchanged.
# Input:  $WORK/vcf/$SAMPLE.norm.vcf.gz
# Output: $WORK/vcf/$SAMPLE.clinvar.vcf.gz (read by 25_vep.sh)
. "$(dirname "$0")/env.sh"

CLINVAR_DATE="20261004"
CLINVAR_DIR="$RES/clinvar"
VCF_DIR="$WORK/vcf"
mkdir -p "$CLINVAR_DIR"

src="$CLINVAR_DIR/clinvar_$CLINVAR_DATE.vcf.gz"
if [[ ! -f "$src" ]]; then
  curl -fL -o "$src" "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2026/clinvar_$CLINVAR_DATE.vcf.gz"
  curl -fL -o "$src.tbi" "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2026/clinvar_$CLINVAR_DATE.vcf.gz.tbi"
fi

chr_clinvar="$CLINVAR_DIR/clinvar_$CLINVAR_DATE.chr.vcf.gz"
if [[ ! -f "$chr_clinvar" ]]; then
  bio bcftools annotate --rename-chrs "$VCF_DIR/chr_rename.txt" "$src" -Oz -o "$chr_clinvar" --write-index=tbi
fi

bio bcftools annotate -a "$chr_clinvar" \
  -c ID,INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN \
  "$VCF_DIR/$SAMPLE.norm.vcf.gz" -Oz -o "$VCF_DIR/$SAMPLE.clinvar.vcf.gz" --write-index=tbi
echo CLINVAR_DONE "$VCF_DIR/$SAMPLE.clinvar.vcf.gz"
