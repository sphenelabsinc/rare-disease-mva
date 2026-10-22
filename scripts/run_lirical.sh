#!/usr/bin/env bash
# Phenotype-only ranking with LIRICAL v2.4.1, run locally (complements Exomiser's combined score).
# HPO release 2026-09-01 (phenotype.hpoa #version 2026-09-02). Pre-set reading: BUB1B in the
# top 20 counts as agreement with Exomiser.
#
# One-time install (into data/resources, gitignored; no admin rights needed):
#   curl -L -o jdk21.tar.gz https://api.adoptium.net/v3/binary/latest/21/ga/mac/x64/jdk/hotspot/normal/eclipse
#   tar -xzf jdk21.tar.gz -C data/resources/jdk
#   curl -L -o lirical.zip https://github.com/TheJacksonLaboratory/LIRICAL/releases/download/v2.4.1/lirical-cli-2.4.1-distribution.zip
#   unzip lirical.zip -d data/resources/lirical
#   java -jar data/resources/lirical/lirical-cli-2.4.1/lirical-cli-2.4.1.jar download -d data/resources/lirical/data
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RES="$ROOT/data/resources"
JAVA="$(ls -d "$RES"/jdk/*/Contents/Home/bin/java | head -1)"
JAR="$RES/lirical/lirical-cli-2.4.1/lirical-cli-2.4.1.jar"
LDATA="$RES/lirical/data"
OUT="$ROOT/data/work/lirical"
HPO="HP:0002859,HP:0000121,HP:0004322,HP:0001508,HP:0003202,HP:0001622,HP:0001518,HP:0200067"

# The current HGNC release maps a few lncRNAs/pseudogenes to a shared Ensembl or Entrez ID
# (e.g. LINC00595/LINC00856 -> ENSG00000230417), which makes LIRICAL abort on a duplicate key.
# Keep the first row for each ID and blank the ID in later duplicates; no MVA-related gene is affected.
HGNC="$LDATA/hgnc_complete_set.txt"
[[ -f "$HGNC.orig" ]] || cp "$HGNC" "$HGNC.orig"
awk -F'\t' -v OFS='\t' '
  NR == 1 { for (i = 1; i <= NF; i++) { if ($i == "entrez_id") e = i; if ($i == "ensembl_gene_id") g = i }; print; next }
  { if ($e != "" && seen_e[$e]++) $e = ""; if ($g != "" && seen_g[$g]++) $g = ""; print }
' "$HGNC.orig" > "$HGNC"

mkdir -p "$OUT"
"$JAVA" -jar "$JAR" prioritize -d "$LDATA" -p "$HPO" -o "$OUT" -x lirical -f tsv -f html 2>&1 | grep -v -E 'SLF4J|DEBUG' || true

# Phenotype-only output ranks diseases; map them to genes with the same HPO release and keep
# each gene's best disease rank -> lirical_genes.tsv (gene_rank is the rank among genes).
G2D="$LDATA/genes_to_disease.txt"
[[ -f "$G2D" ]] || curl -sfL -o "$G2D" \
  https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-09-01/genes_to_disease.txt
"$ROOT/.venv/bin/python" - "$OUT" "$G2D" <<'EOF'
import sys
import pandas as pd
out, g2d = sys.argv[1], sys.argv[2]
d = pd.read_csv(f"{out}/lirical.tsv", sep="\t", comment="!", thousands=",")
m = d.merge(pd.read_csv(g2d, sep="\t")[["gene_symbol", "disease_id"]],
            left_on="diseaseCurie", right_on="disease_id")
genes = m.sort_values("rank").drop_duplicates("gene_symbol").reset_index(drop=True)
genes.insert(0, "gene_rank", genes.index + 1)
cols = ["gene_rank", "gene_symbol", "rank", "diseaseCurie", "diseaseName", "posttestprob", "compositeLR"]
genes[cols].to_csv(f"{out}/lirical_genes.tsv", sep="\t", index=False)
print(genes[genes.gene_symbol == "BUB1B"][cols].to_string(index=False))
EOF
echo "LIRICAL_DONE $OUT"
