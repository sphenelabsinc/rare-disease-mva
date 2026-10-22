"""Genotype-level QC: a genome-wide funnel and a check of the submitted candidates.

Both commands only count and flag. Nothing is removed from the VCF that Exomiser
reads (Exomiser runs with PASS_ONLY), so a variant that misses a threshold here is
kept and flagged for read-level review. Thresholds come from config.py.

Usage:
    python -m rare_disease_mva.qc funnel [--vcf VCF]
    python -m rare_disease_mva.qc candidates <submission.csv> [--vcf VCF]

The VCF defaults to the normalised VCF from 20_prep_vcf.sh when present, otherwise
the raw VCF. Multi-allelic records in the raw VCF are handled per allele.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

from cyvcf2 import VCF

from rare_disease_mva import config as C

STEPS = [
    ("primary_contig", "chr1-22, X, Y, MT (drop decoys, unplaced and HLA contigs)"),
    ("filter_pass", "FILTER == PASS (GATK hard filters)"),
    ("carries_alt", "proband genotype contains a non-reference allele"),
    ("min_gq", f"FORMAT/GQ >= {C.MIN_GQ}"),
    ("min_dp", f"FORMAT/DP >= {C.MIN_DP}"),
    ("allele_balance", f"het AB in {C.HET_AB_RANGE[0]}-{C.HET_AB_RANGE[1]}; hom-alt AB >= {C.HOM_ALT_MIN_AB}"),
]
PRIMARY = set(C.PRIMARY_CONTIGS)


def default_vcf() -> Path:
    return C.NORM_VCF if C.NORM_VCF.exists() else C.RAW_VCF


def bare_contig(chrom: str) -> str:
    """'chr15' -> '15', 'chrM' -> 'MT', so raw and normalised VCFs compare the same."""
    c = chrom.removeprefix("chr")
    return "MT" if c == "M" else c


def fmt_int(v, name):
    """First sample's integer FORMAT value, or None when the field is absent or missing."""
    arr = v.format(name)
    if arr is None:
        return None
    x = int(arr[0][0])
    return x if x >= 0 else None


def alleles_and_ad(v):
    gt = [a for a in v.genotypes[0][:-1] if a >= 0]
    arr = v.format("AD")
    ad = None if arr is None else [int(x) for x in arr[0]]
    if ad is not None and (any(x < 0 for x in ad) or sum(ad) == 0):
        ad = None
    return gt, ad


def allele_balance_ok(gt: list[int], ad: list[int] | None) -> bool:
    """Ported from the Plan's 01_qc_filter. Handles 0/1, 1/1, hemizygous and 1/2 calls."""
    if not ad:
        return True  # no AD to judge by; GQ/DP already gate the call
    total = sum(ad)
    alts = {a for a in gt if a > 0}
    if any(a >= len(ad) for a in alts):
        return True
    if len(set(gt)) == 1:  # hom-alt, or hemizygous call on X/Y
        return ad[gt[0]] / total >= C.HOM_ALT_MIN_AB
    lo, hi = C.HET_AB_RANGE
    if 0 in gt:
        return all(lo <= ad[a] / total <= hi for a in alts)
    return all(ad[a] / total >= lo for a in alts)  # 1/2: every carried alt must look real


def first_failed_step(v) -> str | None:
    if bare_contig(v.CHROM) not in PRIMARY:
        return "primary_contig"
    if C.REQUIRE_PASS and v.FILTER is not None:  # cyvcf2: FILTER is None for PASS
        return "filter_pass"
    gt, ad = alleles_and_ad(v)
    if not any(a > 0 for a in gt):
        return "carries_alt"
    gq, dp = fmt_int(v, "GQ"), fmt_int(v, "DP")
    if gq is None or gq < C.MIN_GQ:
        return "min_gq"
    if dp is None or dp < C.MIN_DP:
        return "min_dp"
    if not allele_balance_ok(gt, ad):
        return "allele_balance"
    return None


def funnel(vcf_path: Path) -> None:
    failed, total = Counter(), 0
    for v in VCF(str(vcf_path)):
        total += 1
        step = first_failed_step(v)
        if step:
            failed[step] += 1
    C.QC_DIR.mkdir(parents=True, exist_ok=True)
    remaining = total
    with open(C.QC_FUNNEL, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["step", "rule", "flagged", "remaining"])
        w.writerow(["input", str(vcf_path.name), 0, total])
        for step, rule in STEPS:
            remaining -= failed[step]
            w.writerow([step, rule, failed[step], remaining])
    print(C.QC_FUNNEL.read_text(), end="")


def submitted_variants(csv_path: Path):
    """Yield (row, finding_type, gene, chrom, pos, ref, alt) for every variant in the CSV."""
    with open(csv_path) as fh:
        for i, r in enumerate(csv.DictReader(fh), start=1):
            gene = (r.get("notes") or "").split(" ", 1)[0]
            for k in ("1", "2"):
                if r.get(f"pos_{k}"):
                    yield (i, r["finding_type"], gene, r[f"chrom_{k}"], int(r[f"pos_{k}"]),
                           r[f"ref_{k}"], r[f"alt_{k}"])


def find_record(vcf: VCF, chrom: str, pos: int, ref: str, alt: str):
    """Exact match on POS/REF/ALT; returns (record, alt_index) or (None, None)."""
    name = chrom if chrom in vcf.seqnames else bare_contig(chrom)
    for v in vcf(f"{name}:{pos}-{pos}"):
        if v.POS == pos and v.REF == ref and alt in v.ALT:
            return v, v.ALT.index(alt) + 1
    return None, None


def candidates(csv_path: Path, vcf_path: Path) -> None:
    vcf = VCF(str(vcf_path))
    lo, hi = C.HET_AB_RANGE
    rows = []
    for row, ftype, gene, chrom, pos, ref, alt in submitted_variants(csv_path):
        v, idx = find_record(vcf, chrom, pos, ref, alt)
        out = {"row": row, "finding_type": ftype, "gene": gene, "variant": f"{chrom}:{pos} {ref}>{alt}"}
        if v is None:
            rows.append(out | {"match": "none", "qc": "check_manually"})
            continue
        gt, ad = alleles_and_ad(v)
        gq, dp = fmt_int(v, "GQ"), fmt_int(v, "DP")
        ab = round(ad[idx] / sum(ad), 3) if ad and idx < len(ad) else None
        flags = []
        if v.FILTER is not None:
            flags.append(f"filter={v.FILTER}")
        if gq is None or gq < C.MIN_GQ:
            flags.append("min_gq")
        if dp is None or dp < C.MIN_DP:
            flags.append("min_dp")
        if ab is not None:
            hom = gt.count(idx) == len(gt)
            if (hom and ab < C.HOM_ALT_MIN_AB) or (not hom and not lo <= ab <= hi):
                flags.append("allele_balance")
        rows.append(out | {
            "match": "exact", "filter": v.FILTER or "PASS",
            "gt": "/".join(map(str, gt)), "gq": gq, "dp": dp,
            "ad": ",".join(map(str, ad)) if ad else "", "ab": ab,
            "flags": ",".join(flags), "qc": "flag" if flags else "pass",
        })
    cols = ["row", "finding_type", "gene", "variant", "match", "filter", "gt", "gq", "dp", "ad", "ab", "flags", "qc"]
    C.QC_DIR.mkdir(parents=True, exist_ok=True)
    with open(C.CANDIDATE_QC, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)
    print(C.CANDIDATE_QC.read_text(), end="")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("funnel", help="genome-wide QC funnel counts")
    f.add_argument("--vcf", type=Path, default=None)
    c = sub.add_parser("candidates", help="GQ/DP/AB check of the submitted variants")
    c.add_argument("csv", type=Path)
    c.add_argument("--vcf", type=Path, default=None)
    a = p.parse_args(argv)
    vcf = a.vcf or default_vcf()
    print(f"# VCF: {vcf}")
    if a.cmd == "funnel":
        funnel(vcf)
    else:
        candidates(a.csv, vcf)


if __name__ == "__main__":
    main()
