"""Merge the mechanism gene panel from the in-repo list, PanelApp and the literature.

The merged panel is used only to check whether top-ranked genes fall in known mechanism genes, and as the control
arm of the EPCR sensitivity analysis. It does not weight the ranking or enter EPCR.

Sources, each recorded in the `sources` column:
  repo:<category>        config/gene_panels.tsv (59 genes chosen from the phenotype)
  panelapp:<id>@<ver>    PanelApp green genes for R457 (1637) and Childhood solid tumours (243)
  literature             genes added from a literature review

Only panel IDs are sent to PanelApp; no patient data leaves the machine.

Usage:
    python -m rare_disease_mva.panel [out.tsv]   (default config/gene_panels_merged.tsv)
"""
import csv
import sys
from pathlib import Path

import requests

from rare_disease_mva import config as C

PANELAPP = "https://panelapp.genomicsengland.co.uk/api/v1/panels/{}/"
PANELAPP_IDS = {1637: "R457 Sarcoma of possible germline origin", 243: "Childhood solid tumours"}
LITERATURE = {
    "CENATAC": ("AR", "MVA (CENATAC); added from a literature review"),
    "ZWINT": ("", "kinetochore / SAC component; added from a literature review"),
}
REPO_PANEL = C.ROOT / "config" / "gene_panels.tsv"
DEFAULT_OUT = C.ROOT / "config" / "gene_panels_merged.tsv"


def panelapp_green(panel_id):
    d = requests.get(PANELAPP.format(panel_id), timeout=60).json()
    genes = sorted(g["gene_data"]["gene_symbol"] for g in d["genes"] if g["confidence_level"] == "3")
    return d["version"], genes


def merge():
    rows = {}
    with open(REPO_PANEL) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows[r["gene"]] = {"gene": r["gene"], "sources": [f"repo:{r['panel']}"],
                               "inheritance": r["inheritance"], "rationale": r["rationale"]}
    for pid, label in PANELAPP_IDS.items():
        version, genes = panelapp_green(pid)
        for g in genes:
            row = rows.setdefault(g, {"gene": g, "sources": [], "inheritance": "",
                                      "rationale": f"PanelApp green: {label}"})
            row["sources"].append(f"panelapp:{pid}@{version}")
    for g, (inh, why) in LITERATURE.items():
        row = rows.setdefault(g, {"gene": g, "sources": [], "inheritance": inh, "rationale": why})
        row["sources"].append("literature")
    return sorted(rows.values(), key=lambda r: r["gene"])


def main(out=DEFAULT_OUT):
    rows = merge()
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene", "sources", "inheritance", "rationale"])
        for r in rows:
            w.writerow([r["gene"], ";".join(r["sources"]), r["inheritance"], r["rationale"]])
    n_multi = sum(len(r["sources"]) > 1 for r in rows)
    print(f"{len(rows)} genes written to {out} ({n_multi} from more than one source)")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_OUT)
