"""EPCR scoring for the Track 1 submission, with an ordering check and a sensitivity analysis.

Primary candidates:
    EPCR = w_pheno * phenotype_score + w_acmg * ACMG tier
where phenotype_score is the unweighted Exomiser gene combined score (0-1) and the ACMG
tier is P 1.0, LP 0.75, VUS 0.4, LB 0.1, B 0. For a compound het the weaker allele's tier
is used, since both alleles must be pathogenic. Two hets closer than CIS_DISTANCE_BP are
treated as one allele in cis (EPCR x CIS_PENALTY). Hedge and secondary rows get fixed
values and skip the formula. A manual adjustment (<= MAX_MANUAL_ADJUST) needs a reason.

Excluded before scoring: any allele classed B/LB, a recessive gene with a single het,
artefact-prone regions, and failed read review.

The mechanism panel enters only the control weight set (0.5/0.3/0.2) of the sensitivity
analysis; it is never part of the submitted EPCR.

Candidate table (TSV), one row per candidate:
    candidate_id  gene  finding_type (primary|hedge|secondary)  inheritance (AR|AD|XR|XD|MT)
    hgvs_1  zygosity_1 (het|hom|hemi)  acmg_1  hgvs_2  acmg_2  distance_bp
    phenotype_score  in_panel (0|1)  artefact (0|1)  read_fail (0|1)  adjust  adjust_reason
Blank hgvs_2/acmg_2/distance_bp for single-variant candidates; acmg may be NA for a
hedge allele without a class.

Usage:
    python -m rare_disease_mva.epcr <candidates.tsv> [--out DIR]
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from rare_disease_mva import config as C

RECESSIVE = ("AR", "XR")
FIXED = {"hedge": C.HEDGE_EPCR, "secondary": C.SECONDARY_EPCR}


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
        return "artefact-prone region"
    if _flag(c, "read_fail"):
        return "failed read review"
    if (c["finding_type"] == "primary" and c["inheritance"] in RECESSIVE
            and not c.get("hgvs_2") and c.get("zygosity_1") not in ("hom", "hemi")):
        return "recessive gene with a single het"
    return None


def score(c: dict, weights=C.EPCR_WEIGHTS) -> float | None:
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
    epcr = w_pheno * float(c["phenotype_score"]) + w_acmg * tier + w_panel * _flag(c, "in_panel")
    if c.get("distance_bp") and float(c["distance_bp"]) < C.CIS_DISTANCE_BP:
        epcr *= C.CIS_PENALTY
    adjust = float(c.get("adjust") or 0)
    if adjust:
        if abs(adjust) > C.MAX_MANUAL_ADJUST or not c.get("adjust_reason", "").strip():
            raise ValueError(f"{c['candidate_id']}: manual adjustment needs |adjust| <= "
                             f"{C.MAX_MANUAL_ADJUST} and a written reason")
        epcr += adjust
    return round(min(max(epcr, 0.0), 1.0), 3)


def score_all(cands: list[dict], weights=C.EPCR_WEIGHTS) -> list[dict]:
    """Scored rows, highest EPCR first; excluded candidates last with their reason."""
    rows = [dict(c, epcr=score(c, weights), excluded=exclusion_reason(c) or "") for c in cands]
    kept = sorted((r for r in rows if r["epcr"] is not None), key=lambda r: -r["epcr"])
    return kept + [r for r in rows if r["epcr"] is None]


def ordering_check(rows: list[dict]) -> None:
    """The top row's EPCR must be strictly higher than every other submitted row."""
    kept = [r for r in rows if r["epcr"] is not None]
    if len(kept) > 1 and kept[0]["epcr"] <= kept[1]["epcr"]:
        raise ValueError(f"ordering check failed: {kept[0]['candidate_id']} ({kept[0]['epcr']}) is not "
                         f"strictly above {kept[1]['candidate_id']} ({kept[1]['epcr']})")


def sensitivity(cands: list[dict], weight_sets=C.EPCR_SENSITIVITY) -> list[dict]:
    """Top candidate under each weight set; the panel-bonus set is marked as the control."""
    base_top = score_all(cands)[0]["candidate_id"]
    out = []
    for ws in weight_sets:
        kept = [r for r in score_all(cands, ws) if r["epcr"] is not None]
        top = kept[0]
        runner = kept[1] if len(kept) > 1 else None
        out.append({
            "weights": "/".join(f"{w:g}" for w in ws),
            "role": "control (panel bonus)" if tuple(ws) == C.EPCR_CONTROL else "weights",
            "top": top["candidate_id"], "top_epcr": top["epcr"],
            "runner_up": runner["candidate_id"] if runner else "",
            "gap": round(top["epcr"] - runner["epcr"], 3) if runner else "",
            "top_unchanged": top["candidate_id"] == base_top,
        })
    return out


def _write(path: Path, rows: list[dict], cols: list[str]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("candidates", type=Path)
    p.add_argument("--out", type=Path, default=C.WORK / "epcr")
    a = p.parse_args(argv)
    cands = load_candidates(a.candidates)
    rows = score_all(cands)
    ordering_check(rows)
    a.out.mkdir(parents=True, exist_ok=True)
    _write(a.out / "epcr_scored.tsv", rows,
           ["candidate_id", "gene", "finding_type", "epcr", "excluded", "phenotype_score",
            "acmg_1", "acmg_2", "distance_bp", "adjust", "adjust_reason"])
    sens = sensitivity(cands)
    _write(a.out / "epcr_sensitivity.tsv", sens, list(sens[0]))
    for r in rows:
        print(r["candidate_id"], r["epcr"], r["excluded"], sep="\t")
    print("\n" + "\n".join("\t".join(str(v) for v in s.values()) for s in sens))


if __name__ == "__main__":
    main()
