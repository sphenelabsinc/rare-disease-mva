.PHONY: env data clean

VENV := .venv
PYTHON := $(VENV)/bin/python
PIP := $(VENV)/bin/pip

## Bootstrap: create venv and install deps
env: $(VENV)/bin/activate

$(VENV)/bin/activate: pyproject.toml
	python3 -m venv $(VENV)
	$(PIP) install --upgrade pip
	$(PIP) install -e ".[dev]"
	@touch $(VENV)/bin/activate
	@echo "Activate with: source $(VENV)/bin/activate"

## Download full dataset from HuggingFace (~85 GB, needs HF_TOKEN in .env)
# Downloads sequentially with per-file retry; skips already-complete files.
# HF_HUB_DISABLE_XET avoids xet protocol crashes on Python 3.14 / macOS.
data: env
	@test -f .env || (echo "ERROR: copy .env.example to .env and set HF_TOKEN" && exit 1)
	HF_HUB_DISABLE_XET=1 $(PYTHON) scripts/download_data.py

clean:
	rm -rf $(VENV) __pycache__ .pytest_cache
