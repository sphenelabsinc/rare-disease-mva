"""Download dataset files one at a time with retry logic.

Usage:
    python scripts/download_data.py

Skips files that already exist at the correct size.
"""
import os
import time
from pathlib import Path

from dotenv import load_dotenv
from huggingface_hub import HfApi, hf_hub_download

load_dotenv()

REPO_ID = "SageBio/mva-hackathon-2026-data"
LOCAL_DIR = Path("./data")
SKIP_PATTERNS = {"README.md", ".gitattributes"}
MAX_RETRIES = 5
RETRY_DELAY = 30  # seconds

token = os.environ.get("HF_TOKEN")
if not token:
    raise SystemExit("ERROR: HF_TOKEN not set in .env")

LOCAL_DIR.mkdir(exist_ok=True)

api = HfApi()
files = [
    f for f in api.list_repo_files(REPO_ID, repo_type="dataset", token=token)
    if f not in SKIP_PATTERNS
]
print(f"Found {len(files)} files to download")

for i, filename in enumerate(files, 1):
    dest = LOCAL_DIR / filename
    # Check if already complete by comparing size
    if dest.exists():
        info = api.get_paths_info(REPO_ID, filename, repo_type="dataset", token=token)
        if info and info[0].size and dest.stat().st_size == info[0].size:
            print(f"[{i}/{len(files)}] SKIP {filename} (complete)")
            continue

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            print(f"[{i}/{len(files)}] Downloading {filename} (attempt {attempt})")
            hf_hub_download(
                repo_id=REPO_ID,
                filename=filename,
                repo_type="dataset",
                local_dir=str(LOCAL_DIR),
                token=token,
            )
            print(f"[{i}/{len(files)}] Done: {filename}")
            break
        except Exception as e:
            print(f"  ERROR: {e}")
            if attempt < MAX_RETRIES:
                print(f"  Retrying in {RETRY_DELAY}s...")
                time.sleep(RETRY_DELAY)
            else:
                print(f"  FAILED after {MAX_RETRIES} attempts: {filename}")

print("Download complete.")
