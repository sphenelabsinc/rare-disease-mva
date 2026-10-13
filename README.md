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

## Layout

```
data/        raw HF dataset (gitignored, never commit)
docs/        specs and instructions
notebooks/   exploratory analysis
scripts/     one-off utilities (data download)
src/         reusable pipeline code (rare_disease_mva package)
```

## Conventions

- Read inputs from `data/`; don't hardcode absolute paths.
- Load secrets from `.env` via `python-dotenv`.
- Parse VCFs with `cyvcf2`.
- Explore in `notebooks/`; move reusable logic into `src/`.
