"""Shared paths and QC thresholds for the Track 1 pipeline.

Every QC threshold lives here so the report can cite one place, and changing a
cutoff never requires editing a script. Ported from the Plan's pipeline/config.py.

QC thresholds only count and flag variants; nothing is removed from the input to
Exomiser (which uses FILTER=PASS via its own PASS_ONLY mode).
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
WORK = DATA / "work"

# --- Inputs -------------------------------------------------------------------
SAMPLE = "WGS_EX2312012"
PROBAND_ID = "PROBAND01"  # as used in track1_submission_template.csv
RAW_VCF = DATA / "WGS_EX2312012_HGWCNDSX7.vcf.gz"
NORM_VCF = WORK / "vcf" / f"{SAMPLE}.norm.vcf.gz"  # from scripts/pipeline/20_prep_vcf.sh
PHENOTYPE_DOCX = DATA / "Challenge_Clinical_Phenotype_1.docx"

# --- Phase 0: reference genome --------------------------------------------------
GRCH38_CHR1_LENGTH = 248_956_422
GRCH37_CHR1_LENGTH = 249_250_621
PHASE0_DIR = WORK / "phase0"

# --- QC funnel (file-internal quality only, no population data) ------------------
PRIMARY_CONTIGS = [str(c) for c in range(1, 23)] + ["X", "Y", "MT"]
REQUIRE_PASS = True          # GATK hard filters: QD2, MQ40, FS60, ...
MIN_GQ = 20                  # genotype quality (phred): ~99% confidence in the call
MIN_DP = 10                  # per-sample read depth
HET_AB_RANGE = (0.20, 0.80)  # alt allele balance for het calls; outside = possible artefact or mosaic
HOM_ALT_MIN_AB = 0.80        # alt allele balance for hom-alt calls

QC_DIR = WORK / "qc"
QC_FUNNEL = QC_DIR / "qc_funnel.tsv"
CANDIDATE_QC = QC_DIR / "candidate_qc.tsv"
