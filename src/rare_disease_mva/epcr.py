"""EPCR scoring for the Track 1 submission, with an ordering check and a sensitivity analysis.

Primary candidates:
    EPCR = w_pheno * phenotype_score + w_acmg * ACMG tier
where phenotype_score is the Exomiser gene phenotype score (0-1) and the ACMG tier is
P 1.0, LP 0.75, VUS 0.4, LB 0.1, B 0. The phenotype term deliberately excludes variant
evidence: the Exomiser combined score already folds in the variant score, and the ACMG tier
scores the variant again, so using the combined score would count variant severity twice. For a compound het the weaker allele's tier
is used, since both alleles must be pathogenic. Two hets closer than CIS_DISTANCE_BP are
treated as one allele in cis (EPCR x CIS_PENALTY). Hedge and secondary rows get fixed
values and skip the formula. A manual adjustment (<= MAX_MANUAL_ADJUST) needs a reason.

Excluded before scoring: any allele classed B/LB, a recessive gene with a single het,
artefact-prone regions, and failed read review.

The mechanism panel enters only the control weight set (0.5/0.3/0.2) of the sensitivity
analysis; it is never part of the submitted EPCR. Two further sensitivity arms change the
phenotype term at 0.6/0.4: the combined score (the earlier definition), and the phenotype
score corrected for the no-match baseline, max(0, (score - b) / (1 - b)) with b = 0.5.

Candidate table (TSV), one row per candidate:
    candidate_id  gene  finding_type (primary|hedge|secondary)  inheritance (AR|AD|XR|XD|MT)
    hgvs_1  zygosity_1 (het|hom|hemi)  acmg_1  hgvs_2  acmg_2  distance_bp
    phenotype_score  [combined_score]  in_panel (0|1)  artefact (0|1)  read_fail (0|1)
    adjust  adjust_reason
Blank hgvs_2/acmg_2/distance_bp for single-variant candidates; acmg may be NA for a
hedge allele without a class. Tables from rank.py also carry chrom/pos/ref/alt_1/_2 and
notes; with --submission these are written as the <=10-row Track 1 CSV in EPCR order
(hedge rows are submitted as primary findings).

Usage:
    python -m rare_disease_mva.epcr <candidates.tsv> [--out DIR] [--submission CSV]
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from rare_disease_mva import config as C

RECESSIVE = ("AR", "XR")
MAX_ROWS = 10
SUBMISSION_COLS = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
                   "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes"]
FIXED = {"hedge": C.HEDGE_EPCR, "secondary": C.SECONDARY_EPCR}
PHENO_TERMS = ("phenotype", "combined", "baseline")
# (weights, phenotype term, role) for the sensitivity table
SENSITIVITY_ARMS = ([(ws, "phenotype", "control (panel bonus)" if ws == C.EPCR_CONTROL else "weights")
                     for ws in C.EPCR_SENSITIVITY]
                    + [(C.EPCR_WEIGHTS, "combined", "previous definition (combined score)"),
                       (C.EPCR_WEIGHTS, "baseline", "baseline-corrected phenotype score")])


def tavtigian_class(points: float) -> str:
    """Tavtigian et al. 2018/2020 point scale."""
    if points >= 10:
        return "P"
    if points >= 6:
        return "LP"
    if points >= 0:
        return "VUS"
    if points >= -6:
        return "LB"
    return "B"


def review_classes(path: Path) -> dict[str, dict]:
    """Sum the reviewed ACMG points per variant: hgvs -> {points, class, pending codes}."""
    out: dict[str, dict] = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            v = out.setdefault(r["hgvs"], {"points": 0.0, "pending": []})
            v["points"] += float(r["points"])
            if r["status"] != "FINAL":
                v["pending"].append(r["code"])
    for v in out.values():
        v["class"] = tavtigian_class(v["points"])
    return out


def load_candidates(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _flag(c: dict, key: str) -> bool:
    return str(c.get(key, "")).strip() in ("1", "true", "True", "yes")


def _classes(c: dict) -> list[str]:
    return [x for x in (c.get("acmg_1", ""), c.get("acmg_2", "")) if x and x != "NA"]


def exclusion_reason(c: dict) -> str | None:
    if any(cls in ("B", "LB") for cls in _classes(c)):
        return "ACMG B/LB"
    if _flag(c, "artefact"):
        return "artefact (prone region, clustered low-AB call or het call on non-PAR chrX)"
    if _flag(c, "read_fail"):
        return "failed read review"
    if (c["finding_type"] == "primary" and c["inheritance"] in RECESSIVE
            and not c.get("hgvs_2") and c.get("zygosity_1") not in ("hom", "hemi")):
        return "recessive gene with a single het"
    return None


def phenotype_term(c: dict, term: str = "phenotype") -> float:
    """Phenotype term of the formula: phenotype score (default), combined score, or baseline-corrected."""
    if term == "combined":
        return float(c["combined_score"])
    pheno = float(c["phenotype_score"])
    if term == "baseline":
        return max(0.0, (pheno - C.PHENO_BASELINE) / (1 - C.PHENO_BASELINE))
    return pheno


def score(c: dict, weights=C.EPCR_WEIGHTS, term: str = "phenotype") -> float | None:
    """EPCR for one candidate, or None if it is excluded."""
    if exclusion_reason(c):
        return None
    if c["finding_type"] in FIXED:
        return FIXED[c["finding_type"]]
    classes = _classes(c)
    if not classes:
        raise ValueError(f"{c['candidate_id']}: primary candidate without an ACMG class")
    w_pheno, w_acmg, w_panel = weights
    tier = min(C.ACMG_TIER[cls] for cls in classes)  # weaker allele for a compound het
    epcr = w_pheno * phenotype_term(c, term) + w_acmg * tier + w_panel * _flag(c, "in_panel")
    if c.get("distance_bp") and float(c["distance_bp"]) < C.CIS_DISTANCE_BP:
        epcr *= C.CIS_PENALTY
    adjust = float(c.get("adjust") or 0)
    if adjust:
        if abs(adjust) > C.MAX_MANUAL_ADJUST or not c.get("adjust_reason", "").strip():
            raise ValueError(f"{c['candidate_id']}: manual adjustment needs |adjust| <= "
                             f"{C.MAX_MANUAL_ADJUST} and a written reason")
        epcr += adjust
    return round(min(max(epcr, 0.0), 1.0), 3)


def score_all(cands: list[dict], weights=C.EPCR_WEIGHTS, term: str = "phenotype") -> list[dict]:
    """Scored rows, highest EPCR first; excluded candidates last with their reason."""
    rows = [dict(c, epcr=score(c, weights, term), excluded=exclusion_reason(c) or "") for c in cands]
    kept = sorted((r for r in rows if r["epcr"] is not None), key=lambda r: -r["epcr"])
    return kept + [r for r in rows if r["epcr"] is None]


def ordering_check(rows: list[dict]) -> None:
    """The top row's EPCR must be strictly higher than every other submitted row."""
    kept = [r for r in rows if r["epcr"] is not None]
    if len(kept) > 1 and kept[0]["epcr"] <= kept[1]["epcr"]:
        raise ValueError(f"ordering check failed: {kept[0]['candidate_id']} ({kept[0]['epcr']}) is not "
                         f"strictly above {kept[1]['candidate_id']} ({kept[1]['epcr']})")


def sensitivity(cands: list[dict], arms=SENSITIVITY_ARMS) -> list[dict]:
    """Top candidate under each (weights, phenotype term) arm; the panel-bonus set is the control."""
    base_top = score_all(cands)[0]["candidate_id"]
    out = []
    for ws, term, role in arms:
        if term == "combined" and not all(c.get("combined_score") for c in cands):
            continue  # older candidate tables carry no combined score
        kept = [r for r in score_all(cands, ws, term) if r["epcr"] is not None]
        top = kept[0]
        runner = kept[1] if len(kept) > 1 else None
        out.append({
            "weights": "/".join(f"{w:g}" for w in ws), "phenotype_term": term, "role": role,
            "top": top["candidate_id"], "top_epcr": top["epcr"],
            "runner_up": runner["candidate_id"] if runner else "",
            "runner_up_epcr": runner["epcr"] if runner else "",
            "gap": round(top["epcr"] - runner["epcr"], 3) if runner else "",
            "top_unchanged": top["candidate_id"] == base_top,
        })
    return out


def _write(path: Path, rows: list[dict], cols: list[str]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def submission_rows(rows: list[dict]) -> list[dict]:
    """Hedge and secondary rows plus the best primary rows, <= MAX_ROWS, in EPCR order.

    Fixed-value rows always keep their place; primary rows fill the remaining slots.
    """
    kept = [r for r in rows if r["epcr"] is not None]
    fixed = [r for r in kept if r["finding_type"] in FIXED]
    primary = [r for r in kept if r["finding_type"] not in FIXED][:MAX_ROWS - len(fixed)]
    kept = sorted(primary + fixed, key=lambda r: -r["epcr"])
    return [dict(r, proband_id=C.PROBAND_ID,
                 finding_type="secondary" if r["finding_type"] == "secondary" else "primary")
            for r in kept]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("candidates", type=Path)
    p.add_argument("--out", type=Path, default=C.WORK / "epcr")
    p.add_argument("--submission", type=Path, help="write the Track 1 CSV here")
    a = p.parse_args(argv)
    cands = load_candidates(a.candidates)
    rows = score_all(cands)
    ordering_check(rows)
    a.out.mkdir(parents=True, exist_ok=True)
    _write(a.out / "epcr_scored.tsv", rows,
           ["candidate_id", "gene", "finding_type", "epcr", "excluded", "phenotype_score",
            "combined_score", "acmg_1", "acmg_2", "distance_bp", "adjust", "adjust_reason"])
    if a.submission:
        with open(a.submission, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=SUBMISSION_COLS, extrasaction="ignore")
            w.writeheader()
            w.writerows(submission_rows(rows))
    sens = sensitivity(cands)
    _write(a.out / "epcr_sensitivity.tsv", sens, list(sens[0]))
    for r in rows:
        print(r["candidate_id"], r["epcr"], r["excluded"], sep="\t")
    print("\n" + "\n".join("\t".join(str(v) for v in s.values()) for s in sens))


if __name__ == "__main__":
    main()
