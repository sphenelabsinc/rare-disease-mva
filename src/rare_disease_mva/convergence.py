"""Agreement check between the genotype+phenotype ranking, the phenotype-only ranking and the
mechanism panel. Reported only; nothing here changes the ranking or EPCR.

For each gene: Exomiser gene rank and combined score (best MOI), LIRICAL gene rank with all
8 HPO terms and without HP:0200067, mechanism-panel membership with sources, and the
best rare variant (MAX_FREQ <= 0.1%, highest Exomiser variant score).

Usage:
    python -m rare_disease_mva.convergence <exomiser.genes.tsv> <exomiser.variants.tsv>
        <submission.csv> [--genes G1,G2,...] [--out TSV]
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from rare_disease_mva import config as C
from rare_disease_mva.rank import ALT_PAIR_MAX_FREQ_PCT, _hgvs

MECHANISM_GENES = ("BUB1B", "CEP57", "TRIP13", "MAD1L1")
LIRICAL_TOP = ("KCNJ1", "SLC12A1", "HRAS")
COLS = ["gene", "exomiser_rank", "exomiser_moi", "exomiser_score", "lirical_rank",
        "lirical_rank_no_parental", "in_panel", "panel_sources", "best_rare_variant"]


def _tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def exomiser_genes(path) -> dict[str, dict]:
    """Best (lowest) rank per gene."""
    out = {}
    for r in _tsv(path):
        out.setdefault(r["GENE_SYMBOL"], r)
    return out


def lirical_ranks(path) -> dict[str, str]:
    out = {}
    for r in _tsv(path):
        out.setdefault(r["gene_symbol"], r["gene_rank"])
    return out


def best_rare_variants(path, genes) -> dict[str, str]:
    best: dict[str, tuple] = {}
    for r in _tsv(path):
        g = r["GENE_SYMBOL"]
        freq = float(r["MAX_FREQ"] or 0)
        if g not in genes or freq > ALT_PAIR_MAX_FREQ_PCT:
            continue
        key = (-float(r["EXOMISER_VARIANT_SCORE"] or 0), freq)
        if g not in best or key < best[g][0]:
            desc = (f"{r['CONTIG']}:{r['START']} {r['REF']}>{r['ALT']} {_hgvs(r['HGVS'])} "
                    f"{r['FUNCTIONAL_CLASS']} {r['GENOTYPE']} AF%={freq:g}")
            best[g] = (key, desc)
    return {g: d for g, (_, d) in best.items()}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("genes_tsv", type=Path)
    p.add_argument("variants_tsv", type=Path)
    p.add_argument("submission", type=Path)
    p.add_argument("--genes", default=",".join(MECHANISM_GENES + LIRICAL_TOP))
    p.add_argument("--panel", type=Path, default=C.ROOT / "config" / "gene_panels_merged.tsv")
    p.add_argument("--lirical", type=Path, default=C.WORK / "lirical")
    p.add_argument("--out", type=Path, default=C.WORK / "rank" / "convergence.tsv")
    a = p.parse_args(argv)

    with open(a.submission) as fh:
        submitted = [r["notes"].split()[0] for r in csv.DictReader(fh)]
    genes = list(dict.fromkeys(submitted + [g for g in a.genes.split(",") if g]))
    exo = exomiser_genes(a.genes_tsv)
    lir = lirical_ranks(a.lirical / "lirical_genes.tsv")
    lir_np = lirical_ranks(a.lirical / "lirical_no_parental_genes.tsv")
    panel = {r["gene"]: r for r in _tsv(a.panel)}
    rare = best_rare_variants(a.variants_tsv, set(genes))

    rows = []
    for g in genes:
        e = exo.get(g, {})
        rows.append({
            "gene": g, "exomiser_rank": e.get("#RANK", "not ranked"), "exomiser_moi": e.get("MOI", ""),
            "exomiser_score": e.get("EXOMISER_GENE_COMBINED_SCORE", ""),
            "lirical_rank": lir.get(g, "not ranked"), "lirical_rank_no_parental": lir_np.get(g, "not ranked"),
            "in_panel": int(g in panel), "panel_sources": panel[g]["sources"] if g in panel else "",
            "best_rare_variant": rare.get(g, "none"),
        })
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(*r.values(), sep="\t")


if __name__ == "__main__":
    main()
