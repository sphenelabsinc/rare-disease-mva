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
python -m rare_disease_mva.convergence <exomiser.genes.tsv> <exomiser.variants.tsv> <out.csv>
                                     # -> data/work/rank/convergence.tsv (agreement check, reported only)
```

`rank.py` builds one hypothesis per gene from Exomiser's contributing variants, ranked by the
unweighted Exomiser gene combined score (recessive genes take the two best hets or a
homozygote; artefact-prone genes are flagged), and adds a hedge row pairing the top recessive
gene's strongest allele with its rarest other het. Variants at or above the rarity threshold
(gnomAD popmax < 0.1% recessive, < 0.01% dominant) are dropped first, and candidate variants in a
cluster of ≥ 3 calls within 10 bp with allele balance < 0.2 are flagged as mapping artefacts
(read from the normalised VCF, or the raw VCF via `--vcf`). Panel membership is recorded but does
not change the order. `epcr.py` scores EPCR = 0.6 × Exomiser gene phenotype score + 0.4 × ACMG tier
(weaker allele for compound hets; the phenotype score, not the combined score, so the variant is
not counted twice), applies the exclusion rules, enforces that the top row is strictly highest,
and reruns the ranking under alternative weights and phenotype-term definitions. The candidate-table format is in the module
docstring. ACMG classes for manually reviewed variants come from `config/acmg_review.tsv`.

### Configuration files

| File | Contents |
|---|---|
| `src/rare_disease_mva/config.py` | paths, QC thresholds, rarity and cluster rules, EPCR weights |
| `config/exomiser/` | Exomiser analysis settings and HPO terms |
| `config/gene_panels.tsv` | phenotype-driven differential gene list |
| `config/gene_panels_merged.tsv` | merged mechanism panel with sources |
| `config/acmg_review.tsv` | manual ACMG/AMP review, one row per evidence code with reasons |

### Tests

```sh
pytest
```

## Layout

```
config/      analysis settings, gene panels, manual ACMG review
config/track2/   curated Track 2 inputs (causal edges, targets, rubric)
data/        raw HF dataset, resources and outputs (gitignored, never commit)
docs/        specs and instructions
notebooks/   exploratory analysis (committed without outputs)
reports/     derived summary tables and general write-ups (patient-specific
             reports are kept local; see "Data policy" below)
reports/track2/  Track 2 stage write-ups and tables (see its own README)
scripts/     variant pipeline (scripts/pipeline/), LIRICAL run, data download
scripts/track2/  Track 2 stage scripts (s0–s4, public APIs only)
src/         reusable pipeline code (rare_disease_mva package)
tests/       unit tests
```

## Data policy and what is not in this repository

This repository is public. **No patient data, genotype data or genome-wide
table is committed, in any form.** `data/` is gitignored in full, and the
patterns at the end of `.gitignore` catch the same file types anywhere in the
tree in case a script writes outside `data/`. Git history has been checked: no
file under `data/` has ever been committed.

The committed files under `reports/` are tables derived from public resources
only (PubMed, AlphaFold DB, Ensembl/Compara, STRING, Reactome, ChEMBL, Open
Targets) and write-ups that discuss BUB1B, its pathway and candidate drugs in
general terms.

**Patient-specific reports are left out of git.** They contain the proband's
individual variant calls and read counts, secondary and carrier findings, and
sample QC metrics, none of which has been published by the family or the
organisers. They were submitted to the hackathon directly and are gitignored
here:

| Kept local | Contents |
|---|---|
| `reports/track1/` | Track 1 methods write-ups, submission CSV and report PDF |
| `reports/track2/*_track2_report.{md,pdf}` | Track 2 submitted report |
| `reports/track2/01_variant_dossier.md` | Stage 1 dossier (variant calls, read support) |

Submission files are named `sphenelabs_*` (report PDFs and Markdown, the Track
1 CSV), and that prefix is gitignored anywhere in the tree, so they stay out of
git even if moved.

The methods are described in the code, configs and the committed stage
write-ups; anything that refers to these files can be regenerated locally from
the challenge dataset.

Genotype-scale data is held locally only and is deleted within 30 days of the
submission window closing.

### Expected `data/` layout

File names and purposes only — none of these files are in the repository. Get
them with `make data` (challenge dataset, requires `HF_TOKEN`) and the
`scripts/pipeline/` steps (references and databases).

```
data/
  Challenge_Clinical_Phenotype_1.docx        case description (not committed)
  WGS_EX2312012_HGWCNDSX7.vcf.gz{,.tbi}      challenge VCF
  WGS_EX2312012_..._L00{1..4}_R{1,2}_001.fastq.gz   raw reads, 4 lanes
  ref/
    GRCh38_no_alt_hs38d1.fa{,.fai,.mmi}      reference, index, minimap2 index
  resources/
    clinvar/      clinvar.vcf.gz, clinvar.chr.vcf.gz (+ .tbi)
    vep/          VEP 113 cache (homo_sapiens/)
    exomiser/     Exomiser 14.0.0 CLI, 2406_hg38, 2406_phenotype
    gencode/      gencode.v46.basic.gtf.gz, genes.bed
    1kg/          chr15 BUB1B 4 Mb panel subset (+ .tbi)
  work/           all pipeline outputs: bam/ vcf/ vep/ exomiser/ sv/ depth/
                  qc/ rank/ epcr/ lirical/ track2/cache/
  logs/           per-step logs
```

### Which committed files depend on which ignored inputs

| Committed | Needs (not committed) |
|---|---|
| `src/rare_disease_mva/qc.py` | `data/work/vcf/*.norm.vcf.gz`, or the raw VCF |
| `src/rare_disease_mva/rank.py`, `epcr.py`, `convergence.py` | Exomiser TSVs under `data/work/exomiser/` |
| `config/acmg_review.tsv` | manual review of the challenge VCF calls |
| `scripts/track2/s0`–`s4` | nothing patient-derived; public APIs, cached in `data/work/track2/cache/` |
| `reports/track2/tables/*.tsv` | regenerated from those public APIs alone |

So Track 2 is fully reproducible from a clean checkout with network access;
Track 1 needs the challenge dataset and the pipeline run first.

### Verifying before a push

```sh
git ls-files | git check-ignore --stdin -v     # expect no output
git log --all --pretty=format: --name-only | sort -u | grep -E '^data/'   # expect no output
```

## Conventions

- Read inputs from `data/`; don't hardcode absolute paths.
- Load secrets from `.env` via `python-dotenv`.
- Parse VCFs with `cyvcf2`.
- Explore in `notebooks/`; move reusable logic into `src/`.
