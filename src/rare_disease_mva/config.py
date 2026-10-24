"""Shared paths and QC thresholds for the Track 1 pipeline.

Every QC threshold lives here so the report can cite one place, and changing a
cutoff never requires editing a script.

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

# --- Candidate rules (see rank.py) -------------------------------------------------
# Population rarity, gnomAD popmax as reported by Exomiser (MAX_FREQ, percent units).
MAX_FREQ_PCT_RECESSIVE = 0.1    # AR / XR: < 0.1%
MAX_FREQ_PCT_DOMINANT = 0.01    # AD / XD / MT: < 0.01%
# Clustered low-allele-balance calls are likely reads mismapped from a paralog. Cluster size and
# window follow GATK VariantFiltration (--cluster-size 3 --cluster-window-size 10).
CLUSTER_SIZE = 3
CLUSTER_WINDOW_BP = 10
CLUSTER_MAX_AB = 0.20           # a clustered call is flagged only if its alt allele balance is below this
# A heterozygous call on non-PAR chrX cannot occur in an XY sample (one X copy); such calls are
# reads from a paralog or pseudogene placed on chrX. Sex from the chrX hom-alt/het ratio.
SAMPLE_SEX = "XY"
CHRX_PAR = [(10_001, 2_781_479), (155_701_383, 156_030_895)]  # GRCh38 PAR1, PAR2

# --- EPCR scoring (see epcr.py) ----------------------------------------------------
# EPCR = w_pheno * phenotype_score + w_acmg * ACMG tier (+ w_panel * panel, control only)
# phenotype_score is the Exomiser gene phenotype score, not the combined score: the combined
# score already contains the variant score, which the ACMG tier also scores.
EPCR_WEIGHTS = (0.6, 0.4, 0.0)
EPCR_SENSITIVITY = [(0.6, 0.4, 0.0), (0.5, 0.5, 0.0), (0.75, 0.25, 0.0), (0.5, 0.3, 0.2)]
EPCR_CONTROL = (0.5, 0.3, 0.2)  # panel bonus added back; reported, never submitted
# Exomiser phenotype score of a gene with no phenotype match (median over all genes in this run);
# the baseline-corrected sensitivity arm uses max(0, (score - baseline) / (1 - baseline)).
PHENO_BASELINE = 0.5
ACMG_TIER = {"P": 1.0, "LP": 0.75, "VUS": 0.4, "LB": 0.1, "B": 0.0}
CIS_DISTANCE_BP = 100   # two hets closer than this are treated as one allele in cis
CIS_PENALTY = 0.3
HEDGE_EPCR = 0.08       # fixed band 0.05-0.1, no formula
SECONDARY_EPCR = 0.10   # fixed band 0.05-0.15, no formula
MAX_MANUAL_ADJUST = 0.10
