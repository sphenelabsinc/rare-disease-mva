#!/usr/bin/env bash
# Phenotype-driven prioritisation with Exomiser 14 (hiPHIVE + OMIM), 8 HPO terms.
. "$(dirname "$0")/env.sh"
EXO="$RES/exomiser"
OUT="$WORK/exomiser"
mkdir -p "$OUT"
sed -e "s|__VCF__|$WORK/vcf/$SAMPLE.norm.vcf.gz|" -e "s|__OUT__|$OUT|" \
  "$ROOT/config/exomiser/genome-analysis.yml" > "$OUT/analysis.yml"
cd "$EXO/exomiser-cli-14.0.0"
# Exomiser reads ./application.properties; generate it so the data path stays out of git
{ cat "$ROOT/config/exomiser/application.properties"; echo "exomiser.data-directory=$EXO/data"; } > application.properties
java -Xmx${EXOMISER_HEAP:-8g} -jar exomiser-cli-14.0.0.jar --analysis "$OUT/analysis.yml"
echo EXOMISER_DONE
