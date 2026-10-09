# rare-disease-mva

Bioinformatics workspace for the **Rare Disease Real Kid: MVA Hackathon 2026**
(https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026).
Submission window: 24 Aug – 24 Oct 2026.

## Environment

Python 3.11+ with a local venv. Bootstrap once:

```sh
make env          # creates .venv and installs deps
source .venv/bin/activate
```

To download the dataset (~85 GB):

```sh
cp .env.example .env   # add your HF_TOKEN
make data
```

Data lives in `./data/` and is git-ignored. Never commit it.

## Project layout

```
data/          # raw HF dataset (gitignored)
docs/          # specs and instructions
notebooks/     # exploratory analysis
src/           # reusable pipeline code
```

## Key dataset files

- `WGS_EX*.vcf.gz` – whole-genome sequencing variant call files
- Dataset repo: `SageBio/mva-hackathon-2026-data` (HuggingFace, type=dataset)

## Conventions

- All analysis starts from files in `data/`; never hardcode absolute paths
- Load secrets via `python-dotenv` from `.env`
- Use `cyvcf2` for VCF parsing (faster than PyVCF)
- Notebooks go in `notebooks/`, reusable logic moves to `src/`
