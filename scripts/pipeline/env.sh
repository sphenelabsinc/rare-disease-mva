# Shared paths/tools for the pipeline. Source from any step: . scripts/pipeline/env.sh
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DATA="$ROOT/data"
REF="$DATA/ref/GRCh38_no_alt_hs38d1.fa"
WORK="$DATA/work"
RES="$DATA/resources"
LOGS="$DATA/logs"
SAMPLE="WGS_EX2312012"
RAW_VCF="$DATA/WGS_EX2312012_HGWCNDSX7.vcf.gz"
THREADS="${THREADS:-16}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
MM="${MICROMAMBA:-$HOME/.local/bin/micromamba}"
# Put the bioconda env first on PATH (so vep finds the env perl, not /usr/bin/perl)
export PATH="$MAMBA_ROOT_PREFIX/envs/bio/bin:$PATH"
bio() { "$@"; }
mkdir -p "$WORK" "$LOGS"
