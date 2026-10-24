# rare-disease-mva

Bioinformatics workspace for the
[Rare Disease, Real Kid: MVA Hackathon 2026](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026).
Submission window: 24 Aug – 24 Oct 2026.

## Prerequisites

- Python 3.11+ (tested on 3.13)
- `make`
- ~150 GB free disk (85 GB dataset plus caches and intermediates)
- A HuggingFace account with access to
  [`SageBio/mva-hackathon-2026-data`](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data)

On Debian/Ubuntu, `venv` and `make` ship as separate packages:

```sh
sudo apt install -y python3-venv make
```

## Getting started

```sh
# 1. Create the virtualenv and install dependencies
make env
source .venv/bin/activate

# 2. Add your HuggingFace token (https://huggingface.co/settings/tokens)
cp .env.example .env
$EDITOR .env           # set HF_TOKEN=hf_...

# 3. Download the dataset (~85 GB, 11 files) into ./data
make data
```

`make data` runs `scripts/download_data.py`, which downloads one file at a
time, retries failures, and skips files that are already complete. You can
safely interrupt it and run it again to resume.

To fetch a single file instead:

```sh
hf download SageBio/mva-hackathon-2026-data WGS_EX2312012_HGWCNDSX7.vcf.gz \
  --repo-type dataset --local-dir ./data
```

## Analysis

All outputs go to `data/work/` and tools/databases to `data/resources/` (both gitignored).
Paths and thresholds live in `src/rare_disease_mva/config.py`; the shell steps share
`scripts/pipeline/env.sh`.

### 1. Variant pipeline (workstation, bioconda tools)

Run in numeric order. Needs the reference (`data/ref/`), the Exomiser 14.0.0 data release 2406
and the VEP 113 cache under `data/resources/`.

| Script | Step | Main output |
|---|---|---|
| `scripts/pipeline/10_align.sh` | FASTQ → BAM (minimap2), duplicate marking | `data/work/bam/` |
| `scripts/pipeline/20_prep_vcf.sh` | chr names, primary contigs, split multiallelics, left-align; QC stats, sex check, ROH | `data/work/vcf/*.norm.vcf.gz` |
| `scripts/pipeline/22_clinvar.sh` | ClinVar 2026-10-04 annotation | `data/work/vcf/*.clinvar.vcf.gz` |
| `scripts/pipeline/25_vep.sh` | Ensembl VEP 113 annotation | `data/work/vep/` |
| `scripts/pipeline/30_exomiser.sh` | Exomiser genotype + phenotype ranking (8 HPO terms) | `data/work/exomiser/` |
| `scripts/pipeline/40_sv_cnv.sh` | Manta SVs and mosdepth depth | `data/work/sv/`, `data/work/depth/` |

### 2. Genotype QC (laptop)

Counts and flags only; nothing is removed from the Exomiser input. Reads the normalised VCF
if present, otherwise the raw VCF.

```sh
python -m rare_disease_mva.qc funnel                         # genome-wide funnel -> data/work/qc/qc_funnel.tsv
python -m rare_disease_mva.qc candidates <submission.csv>    # GQ/DP/AB per submitted variant -> data/work/qc/candidate_qc.tsv
```

### 3. Phenotype-only ranking (LIRICAL)

`scripts/run_lirical.sh` runs LIRICAL v2.4.1 with the HPO 2026-09-01 release, once with all
8 HPO terms and once without the parental term HP:0200067, and maps diseases to genes.
One-time install steps (portable JDK 21, LIRICAL, its data) are in the script header.
Output: `data/work/lirical/lirical_genes.tsv`, `lirical_no_parental_genes.tsv`.

### 4. Mechanism gene panel

```sh
python -m rare_disease_mva.panel     # -> config/gene_panels_merged.tsv
```

Merges `config/gene_panels.tsv` with PanelApp green genes (R457, Childhood solid tumours),
the genes of HPO-annotated MVA diseases and literature genes, recording each source. Used to
check whether top genes are known mechanism genes; it does not weight the ranking or EPCR.

### 5. Candidates and EPCR

```sh
python -m rare_disease_mva.rank <exomiser.variants.tsv> config/gene_panels_merged.tsv --secondary LZTR1
                                     # -> data/work/rank/raw_ranking.tsv, candidates.tsv, exclusions.tsv
python -m rare_disease_mva.epcr data/work/rank/candidates.tsv --submission <out.csv>
                                     # -> data/work/epcr/epcr_scored.tsv, epcr_sensitivity.tsv
python -m rare_disease_mva.convergence <exomiser.genes.tsv> <exomiser.variants.tsv> <submission.csv>
                                     # -> data/work/rank/convergence.tsv (agreement check, reported only)
```

`rank.py` builds one hypothesis per gene from Exomiser's contributing variants, ranked by the
unweighted Exomiser gene combined score (recessive genes take the two best hets or a
homozygote; artefact-prone genes are flagged), and adds a hedge row pairing the top recessive
gene's strongest allele with its rarest other het. Variants at or above the rarity threshold
(gnomAD popmax < 0.1% recessive, < 0.01% dominant) are dropped first, and candidate variants in a
cluster of ≥ 3 calls within 10 bp with allele balance < 0.2 are flagged as mapping artefacts
(read from the normalised VCF, or the raw VCF via `--vcf`). In an XY sample, heterozygous calls on
non-PAR chrX are also flagged as mapping artefacts, and a variant marked `fail` in
`config/read_review.tsv` (with a written reason) is excluded as a failed read review. Panel membership is recorded but does
not change the order. `epcr.py` scores EPCR = 0.6 × Exomiser gene phenotype score + 0.4 × ACMG tier
(weaker allele for compound hets; the phenotype score, not the combined score, so the variant is
not counted twice), applies the exclusion rules, enforces that the top row is strictly highest,
and reruns the ranking under alternative weights and phenotype-term definitions. The candidate-table format is in the module
docstring. ACMG classes for manually reviewed variants come from `config/acmg_review.tsv`.

### Configuration files

| File | Contents |
|---|---|
| `src/rare_disease_mva/config.py` | paths, QC thresholds, rarity, cluster and chrX ploidy rules, EPCR weights |
| `config/exomiser/` | Exomiser analysis settings and HPO terms |
| `config/gene_panels.tsv` | phenotype-driven differential gene list |
| `config/gene_panels_merged.tsv` | merged mechanism panel with sources |
| `config/acmg_review.tsv` | manual ACMG/AMP review, one row per evidence code with reasons |
| `config/read_review.tsv` | read-level review decisions (pass/fail) with reasons |

### Tests

```sh
pytest
```

## Layout

```
config/      analysis settings, gene panels, manual ACMG review
data/        raw HF dataset, resources and outputs (gitignored, never commit)
docs/        specs and instructions
notebooks/   exploratory analysis (committed without outputs)
scripts/     variant pipeline (scripts/pipeline/), LIRICAL run, data download
src/         reusable pipeline code (rare_disease_mva package)
tests/       unit tests
```

## Conventions

- Read inputs from `data/`; don't hardcode absolute paths.
- Load secrets from `.env` via `python-dotenv`.
- Parse VCFs with `cyvcf2`.
- Explore in `notebooks/`; move reusable logic into `src/`.
