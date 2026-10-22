"""Merge the mechanism gene panel from the in-repo list, PanelApp, HPO and the literature.

The merged panel is used only to check whether top-ranked genes fall in known mechanism
genes, and as the control arm of the EPCR sensitivity analysis. It does not weight the
ranking or enter EPCR.

Sources, each recorded in the `sources` column:
  repo:<category>        config/gene_panels.tsv (59 genes chosen from the phenotype)
  panelapp:<id>@<ver>    PanelApp green genes for R457 (1637) and Childhood solid tumours (243)
  hpo:<release>          genes of every disease named "Mosaic variegated aneuploidy" in the
                         HPO annotations (phenotype.hpoa + genes_to_disease.txt)
  literature             genes added from a literature review

Only panel IDs are sent to PanelApp and the HPO files are public downloads; no patient
data leaves the machine.

Usage:
    python -m rare_disease_mva.panel [out.tsv]   (default config/gene_panels_merged.tsv)
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import requests

from rare_disease_mva import config as C

PANELAPP = "https://panelapp.genomicsengland.co.uk/api/v1/panels/{}/"
PANELAPP_IDS = {1637: "R457 Sarcoma of possible germline origin", 243: "Childhood solid tumours"}
HPO_RELEASE = "v2026-09-01"
HPO_URL = "https://github.com/obophenotype/human-phenotype-ontology/releases/download/{}/{}"
HPO_DIR = C.DATA / "resources" / "hpo" / HPO_RELEASE
MVA_NAME = "mosaic variegated aneuploidy"
INHERITANCE = {"HP:0000007": "AR", "HP:0000006": "AD", "HP:0001417": "XL"}
LITERATURE = {
    "CENATAC": ("AR", "MVA (CENATAC); added from a literature review"),
    "ZWINT": ("", "kinetochore / SAC component; added from a literature review"),
    "SLF2": ("AR", "Atelis syndrome 1 (OMIM:620184); MVA-like chromosomal instability, literature review"),
    "SMC5": ("AR", "Atelis syndrome 2 (OMIM:620185); MVA-like chromosomal instability, literature review"),
}
REPO_PANEL = C.ROOT / "config" / "gene_panels.tsv"
DEFAULT_OUT = C.ROOT / "config" / "gene_panels_merged.tsv"


def panelapp_green(panel_id):
    d = requests.get(PANELAPP.format(panel_id), timeout=60).json()
    genes = sorted(g["gene_data"]["gene_symbol"] for g in d["genes"] if g["confidence_level"] == "3")
    return d["version"], genes


def hpo_file(name):
    path = HPO_DIR / name
    if not path.exists():
        HPO_DIR.mkdir(parents=True, exist_ok=True)
        r = requests.get(HPO_URL.format(HPO_RELEASE, name), timeout=300)
        r.raise_for_status()
        path.write_bytes(r.content)
    return path


def hpo_mva_genes():
    """gene -> (sorted disease IDs, inheritance), for diseases whose HPO annotation name says MVA."""
    names, modes = {}, defaultdict(set)
    with open(hpo_file("phenotype.hpoa")) as fh:
        rows = csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t")
        for r in rows:
            if MVA_NAME in r["disease_name"].lower():
                names[r["database_id"]] = r["disease_name"]
                if r["aspect"] == "I" and r["hpo_id"] in INHERITANCE:
                    modes[r["database_id"]].add(INHERITANCE[r["hpo_id"]])
    genes = defaultdict(set)
    with open(hpo_file("genes_to_disease.txt")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["disease_id"] in names:
                genes[r["gene_symbol"]].add(r["disease_id"])
    return {g: (sorted(ids), "/".join(sorted(set().union(*(modes[i] for i in ids)))))
            for g, ids in genes.items()}


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
    for g, (ids, inh) in sorted(hpo_mva_genes().items()):
        row = rows.setdefault(g, {"gene": g, "sources": [], "inheritance": inh,
                                  "rationale": f"HPO MVA annotation: {', '.join(ids)}"})
        row["sources"].append(f"hpo:{HPO_RELEASE}")
    for g, (inh, why) in LITERATURE.items():
        row = rows.setdefault(g, {"gene": g, "sources": [], "inheritance": inh, "rationale": why})
        row["inheritance"] = row["inheritance"] or inh
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
