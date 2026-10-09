You are a bioinformatician specializing in rare disease participating in https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026

Step 1 - download data locally. Do NOT commit (make sure to gitignore) the actual data.

Rare Disease, Real Kid: MVA Hackathon 2026 - Dataset
Challenge Space: SageBio/rare-disease-real-kid-mva-hackathon-2026
Submission Period: 24 August 2026 – 24 October 2026
Dataset Size: ~85 GB
Quickstart
Python
from huggingface_hub import hf_hub_download

# Download a specific file into a local './data' folder
file_path = hf_hub_download(
    repo_id="SageBio/mva-hackathon-2026-data",
    filename="WGS_EX2312012_HGWCNDSX7.vcf.gz",  # Replace with target file
    repo_type="dataset",
    local_dir="./data",
    token="YOUR_HF_TOKEN"
)

print(f"File downloaded to: {file_path}")

CLI
# Download a specific file into a local './data' folder
hf download SageBio/mva-hackathon-2026-data \
  WGS_EX2312012_HGWCNDSX7.vcf.gz \
  --repo-type dataset \
  --local-dir ./data \
  --token YOUR_HF_TOKEN

Note: YOUR_HF_TOKEN can be omitted from either approach if you've previously logged in (either via hf auth login or by setting the HF_TOKEN environment variable).

Full Data Download
The compressed dataset is ~85 GB across 11 files. We strongly recommend having at least 100 - 150 GB of storage to account for file caching, indexing, and pipeline intermediate files.

Python
from huggingface_hub import snapshot_download

# Download challenge files into a local './data' folder
dataset_dir = snapshot_download(
    repo_id="SageBio/mva-hackathon-2026-data",
    repo_type="dataset",
    local_dir="./data",
    ignore_patterns=["README.md", ".gitattributes"],
    token="YOUR_HF_TOKEN"
)

print(f"Dataset downloaded to: {dataset_dir}")

CLI
# Download challenge files into a local './data' folder
hf download \
  SageBio/mva-hackathon-2026-data \
  --repo-type dataset \
  --local-dir ./data \
  --exclude "README.md" ".gitattributes" \
  --token YOUR_HF_TOKEN