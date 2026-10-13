#!/usr/bin/env bash
# Annotate the normalised, ClinVar-tagged VCF with Ensembl VEP 113 (offline cache):
# consequences on the MANE/canonical transcript, HGVS, gnomAD v4 exome/genome AFs,
# SIFT/PolyPhen. Output: $WORK/vep/$SAMPLE.vep.vcf.gz
. "$(dirname "$0")/env.sh"

CACHE="$RES/vep"
OUT="$WORK/vep"
mkdir -p "$OUT"

# Check a chromosome dir, not the cache root: a partial (per-chromosome) extract may exist.
if [[ ! -d "$CACHE/homo_sapiens/113_GRCh38/1" ]]; then
  tar -xzf "$CACHE/vep113.tar.gz" -C "$CACHE"
fi

# vep has no chrM/"chr" problem: the cache uses 1..22 names, VEP maps chr-prefixed input itself.
vep --offline --cache --dir_cache "$CACHE" --cache_version 113 --assembly GRCh38 \
  --fasta "$REF" --fork "${VEP_FORKS:-6}" --buffer_size 20000 \
  --input_file "$WORK/vcf/$SAMPLE.clinvar.vcf.gz" --format vcf \
  --vcf --compress_output bgzip --output_file "$OUT/$SAMPLE.vep.vcf.gz" --force_overwrite \
  --symbol --hgvs --canonical --mane --biotype --numbers --variant_class \
  --sift b --polyphen b --af_gnomade --af_gnomadg --max_af \
  --pick --pick_order mane_select,canonical,appris,tsl,biotype,ccds,rank,length \
  --no_stats
bio tabix -f -p vcf "$OUT/$SAMPLE.vep.vcf.gz"
echo VEP_DONE
