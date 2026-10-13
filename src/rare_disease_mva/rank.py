"""Build the Track 1 submission (<=10 rows) from Exomiser output + the curated gene panel.

Each candidate is a gene-level hypothesis:
  * AR / XR genes -> a pair of contributing heterozygous variants (compound het) or
    one homozygous / hemizygous variant.
  * AD / XD genes -> a single variant.

Score = Exomiser gene combined score, boosted for genes on the phenotype-driven panel
(config/gene_panels.tsv). Variants failing the read-level review (manual_overrides)
are dropped. EPCR is assigned by rank with an explicit, documented ladder so the
top hypothesis carries most of the probability mass (F-max is computed at the
variant level, so a sharp top row is rewarded).

Usage:
    python -m rare_disease_mva.rank <exomiser.variants.tsv> <gene_panels.tsv> <out.csv>
"""
from __future__ import annotations

import csv
import sys
from collections import defaultdict
from dataclasses import dataclass, field

PROBAND_ID = "PROBAND01"
MAX_ROWS = 10
# Genes whose calls are dominated by reference/paralog artefacts in short-read WGS.
ARTEFACT_PRONE = ("HLA-", "KIR", "PRSS1", "MUC", "OR")
# Panel weight multipliers on the Exomiser combined score.
PANEL_BOOST = {"mva": 1.5, "mva_related": 1.2, "chromosome_instability": 1.2,
               "sarcoma_predisposition": 1.1, "growth_cancer": 1.1, "nephrocalcinosis": 1.0}
# Primary rows: the top hypothesis is given most of the mass; the rest decay.
EPCR_LADDER = [0.85, 0.06, 0.03, 0.02, 0.015, 0.01, 0.008, 0.006, 0.004, 0.002]
# Two hets closer than this are most likely one complex allele in cis, not a compound het.
CIS_DISTANCE = 100
CIS_PENALTY = 0.3
# Alternative second alleles for the top recessive gene must be at least this rare (percent, Exomiser units).
ALT_PAIR_MAX_FREQ_PCT = 0.1

Variant = tuple[str, int, str, str]


@dataclass
class Candidate:
    gene: str
    moi: str
    score: float
    variants: list[Variant]
    hgvs: list[str]
    panel: str = ""
    finding_type: str = "primary"
    notes: list[str] = field(default_factory=list)


def load_panel(path):
    with open(path) as fh:
        return {r["gene"]: r for r in csv.DictReader(fh, delimiter="\t")}


def _variant(r) -> Variant:
    chrom = r["CONTIG"] if r["CONTIG"].startswith("chr") else f"chr{r['CONTIG']}"
    return (chrom, int(r["START"]), r["REF"], r["ALT"])


def load_exomiser(path):
    """Contributing variants grouped by (gene, MOI), with the gene combined score.

    Also returns every retained (contributing or not) rare variant per gene, used to
    propose alternative compound-het pairings for the top recessive gene.
    """
    groups = defaultdict(lambda: {"score": 0.0, "vars": {}})
    all_rare = defaultdict(dict)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            info = (r["GENOTYPE"], r["HGVS"].split(":", 2)[-1], r["EXOMISER_ACMG_CLASSIFICATION"])
            if float(r["MAX_FREQ"] or 0) <= ALT_PAIR_MAX_FREQ_PCT:
                all_rare[r["GENE_SYMBOL"]][_variant(r)] = info
            if r["CONTRIBUTING_VARIANT"] != "1":
                continue
            g = groups[(r["GENE_SYMBOL"], r["MOI"])]
            g["score"] = max(g["score"], float(r["EXOMISER_GENE_COMBINED_SCORE"]))
            g["vars"][_variant(r)] = info
    return groups, all_rare


def build_candidates(groups, panel, exclude=frozenset()):
    cands = []
    for (gene, moi), g in groups.items():
        if gene.startswith(ARTEFACT_PRONE):
            continue
        vs = {v: info for v, info in g["vars"].items() if v not in exclude}
        if not vs:
            continue
        recessive = moi in ("AR", "XR")
        homs = [v for v, (gt, *_ ) in vs.items() if gt.replace("|", "/") in ("1/1", "1")]
        if recessive and not homs and len(vs) < 2:
            continue  # recessive gene needs two hits
        chosen = homs[:1] if homs else sorted(vs, key=lambda v: v[1])[:2 if recessive else 1]
        p = panel.get(gene)
        boost = PANEL_BOOST.get(p["panel"], 1.0) if p else 1.0
        c = Candidate(gene, moi, g["score"] * boost, chosen,
                      [f"{vs[v][1]} [{vs[v][2]}]" for v in chosen],
                      panel=p["panel"] if p else "")
        if len(chosen) == 2 and abs(chosen[0][1] - chosen[1][1]) < CIS_DISTANCE:
            c.score *= CIS_PENALTY
            c.notes.append(f"hets {abs(chosen[0][1] - chosen[1][1])} bp apart: likely one allele in cis")
        cands.append(c)
    # One hypothesis per gene: keep its best-scoring MOI
    best = {}
    for c in sorted(cands, key=lambda c: -c.score):
        best.setdefault(c.gene, c)
    return sorted(best.values(), key=lambda c: -c.score)


def alternative_pairs(top: Candidate, all_rare, n=1):
    """For a recessive top gene, pair its strongest allele with other rare hets in the gene.

    Hedge for the case where the true second allele is non-coding (e.g. a hypomorphic
    intronic/regulatory variant) rather than the one Exomiser chose.
    """
    if top.moi not in ("AR", "XR") or len(top.variants) != 2:
        return []
    anchor = top.variants[0]
    others = [(v, info) for v, info in sorted(all_rare.get(top.gene, {}).items())
              if v not in top.variants and info[0].replace("|", "/") == "0/1"]
    alts = []
    for v, info in others[:n]:
        c = Candidate(top.gene, top.moi, top.score * 0.05, [anchor, v],
                      [top.hgvs[0], f"{info[1]} [{info[2]}]"], panel=top.panel)
        c.notes.append("alternative second allele (hedge)")
        alts.append(c)
    return alts


def write_submission(cands, out_path):
    rows = []
    for i, c in enumerate(cands[:MAX_ROWS]):
        v1, v2 = c.variants[0], (c.variants[1] if len(c.variants) > 1 else ("", "", "", ""))
        note = f"{c.gene} {c.moi}; " + " + ".join(c.hgvs)
        if c.panel:
            note += f"; panel={c.panel}"
        note += "".join(f"; {n}" for n in c.notes)
        rows.append([PROBAND_ID, *v1, *v2, EPCR_LADDER[i], c.finding_type, note])
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
                    "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes"])
        w.writerows(rows)
    return rows


def main(exo, panel_path, out, secondary_genes=()):
    groups, all_rare = load_exomiser(exo)
    cands = build_candidates(groups, load_panel(panel_path))
    top, rest = cands[0], cands[1:]
    ordered = [top, *alternative_pairs(top, all_rare), *rest]
    for c in ordered:
        if c.gene in secondary_genes:
            c.finding_type = "secondary"
    rows = write_submission(ordered, out)
    for r in rows:
        print(*r, sep="\t")


if __name__ == "__main__":
    main(*sys.argv[1:4], secondary_genes=tuple(sys.argv[4].split(",")) if len(sys.argv) > 4 else ())
