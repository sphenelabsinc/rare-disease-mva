"""Build the Track 1 candidate table from Exomiser output.

Each candidate is a gene-level hypothesis:
  * AR / XR genes -> the two best contributing heterozygous variants (compound het) or
    one homozygous / hemizygous variant.
  * AD / XD genes -> the single best variant.
"Best" means highest Exomiser variant score, then lowest MAX_FREQ.

Before candidates are built, variants at or above the rarity threshold for the MOI are
dropped (gnomAD popmax < 0.1% recessive, < 0.01% dominant; secondary genes exempt). After
selection, each chosen variant is checked in the VCF for the clustered low-allele-balance
pattern of reads mismapped from a paralog (>= 3 calls in a 10 bp window and alt allele
balance < 0.2); a match flags the candidate as an artefact. In an XY sample, a heterozygous
call on non-PAR chrX is also flagged as an artefact (one X copy cannot be heterozygous).
A variant marked `fail` in the read-level review (config/read_review.tsv, with a written
reason) sets read_fail; this is the only way a candidate is removed by hand.

Candidates are ranked by the unweighted Exomiser gene combined score, one MOI per gene; the
top N by that score are passed to epcr.py. The table carries the Exomiser gene phenotype
score as phenotype_score (the EPCR phenotype term) and the combined score for reference.
Membership of the mechanism panel (config/gene_panels_merged.tsv) is recorded for the
agreement check and the EPCR control weights only; it never changes the order.

ACMG classes come from the manual review (config/acmg_review.tsv) where a variant was
reviewed, otherwise from Exomiser's automatic classification; the notes label each class
"manual" or "automatic". The notes also carry any genotype QC flag of a chosen variant (from
qc.py; flagged, never removed) and the reason for each secondary finding (config.py).
EPCR is assigned by epcr.py.

Outputs (in --out):
    raw_ranking.tsv   every gene-level hypothesis, highest score first
    candidates.tsv    the top primary candidates, a hedge row for the top recessive gene and
                      any secondary genes, in the epcr.py table format plus coordinates/notes
    exclusions.tsv    variants removed by the rarity rule (genes within the top-N score range)
                      and candidate variants flagged by the cluster or chrX ploidy rule or
                      failed in the read review

Usage:
    python -m rare_disease_mva.rank <exomiser.variants.tsv> <gene_panels_merged.tsv>
        [--out DIR] [--top N] [--secondary GENE,...] [--vcf VCF]
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from cyvcf2 import VCF

from rare_disease_mva import config as C
from rare_disease_mva.epcr import review_classes
from rare_disease_mva.qc import find_record, genotype_flags

# Genes whose calls are dominated by reference/paralog artefacts in short-read WGS.
# Olfactory receptors are matched as OR<family><subfamily> so ORC1/4/6 and ORAI1 are kept.
ARTEFACT_PRONE = re.compile(r"^(HLA-|KIR|PRSS1|MUC|OR\d+[A-Z])")
# Hedge second alleles for the top recessive gene must be at least this rare (percent, Exomiser units).
ALT_PAIR_MAX_FREQ_PCT = 0.1
EXOMISER_ACMG = {"PATHOGENIC": "P", "LIKELY_PATHOGENIC": "LP", "UNCERTAIN_SIGNIFICANCE": "VUS",
                 "LIKELY_BENIGN": "LB", "BENIGN": "B"}
RECESSIVE = ("AR", "XR")

Variant = tuple[str, int, str, str]


@dataclass
class Call:
    gt: str
    hgvs: str
    acmg: str            # Exomiser automatic class (P/LP/VUS/LB/B/NA)
    score: float         # Exomiser variant score
    freq: float          # MAX_FREQ, percent

    @property
    def zygosity(self) -> str:
        gt = self.gt.replace("|", "/")
        return "hom" if gt == "1/1" else "hemi" if gt == "1" else "het"

    def key(self):
        return (-self.score, self.freq)


@dataclass
class Candidate:
    gene: str
    moi: str
    score: float         # Exomiser gene combined score
    pheno: float         # Exomiser gene phenotype score
    variants: list[Variant]
    calls: list[Call]
    in_panel: bool = False
    panel_sources: str = ""
    finding_type: str = "primary"
    notes: list[str] = field(default_factory=list)
    cluster: bool = False
    ploidy: bool = False
    read_fail: bool = False

    @property
    def artefact(self) -> bool:
        return self.cluster or self.ploidy or bool(ARTEFACT_PRONE.match(self.gene))

    @property
    def distance(self) -> int | None:
        return abs(self.variants[0][1] - self.variants[1][1]) if len(self.variants) == 2 else None


def load_panel(path):
    with open(path) as fh:
        return {r["gene"]: r for r in csv.DictReader(fh, delimiter="\t")}


def _variant(r) -> Variant:
    chrom = r["CONTIG"] if r["CONTIG"].startswith("chr") else f"chr{r['CONTIG']}"
    return (chrom, int(r["START"]), r["REF"], r["ALT"])


def _hgvs(raw: str) -> str:
    """'GENE:TX:c.1A>G:p.(Met1?)' -> 'c.1A>G p.(Met1?)' (the form used in acmg_review.tsv)."""
    c, _, p = raw.split(":", 2)[-1].partition(":")
    return c if p in ("", "p.(=)") else f"{c} {p}"


def _call(r) -> Call:
    return Call(r["GENOTYPE"], _hgvs(r["HGVS"]),
                EXOMISER_ACMG.get(r["EXOMISER_ACMG_CLASSIFICATION"], "NA"),
                float(r["EXOMISER_VARIANT_SCORE"] or 0), float(r["MAX_FREQ"] or 0))


def load_exomiser(path):
    """Contributing variants grouped by (gene, MOI), with the gene combined and phenotype scores.

    Also returns every retained (contributing or not) rare variant per gene, used to
    propose a hedge compound-het pairing for the top recessive gene.
    """
    groups = defaultdict(lambda: {"score": 0.0, "pheno": 0.0, "vars": {}})
    all_rare = defaultdict(dict)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            call = _call(r)
            if call.freq <= ALT_PAIR_MAX_FREQ_PCT:
                all_rare[r["GENE_SYMBOL"]][_variant(r)] = call
            if r["CONTRIBUTING_VARIANT"] != "1":
                continue
            g = groups[(r["GENE_SYMBOL"], r["MOI"])]
            g["score"] = max(g["score"], float(r["EXOMISER_GENE_COMBINED_SCORE"]))
            g["pheno"] = max(g["pheno"], float(r["EXOMISER_GENE_PHENO_SCORE"]))
            g["vars"][_variant(r)] = call
    return groups, all_rare


def max_freq_pct(moi: str) -> float:
    return C.MAX_FREQ_PCT_RECESSIVE if moi in RECESSIVE else C.MAX_FREQ_PCT_DOMINANT


def rarity_filter(groups, exempt=frozenset()):
    """Drop variants at or above the MOI's rarity threshold; returns (groups, exclusions)."""
    kept, excluded = {}, []
    for (gene, moi), g in groups.items():
        limit = max_freq_pct(moi)
        vs = g["vars"] if gene in exempt else {v: c for v, c in g["vars"].items() if c.freq < limit}
        excluded += [{"gene": gene, "moi": moi, "variant": "-".join(map(str, v)), "hgvs": c.hgvs,
                      "rule": "rarity", "score": g["score"],
                      "detail": f"MAX_FREQ {c.freq:.4g}% >= {limit:g}% ({moi})"}
                     for v, c in g["vars"].items() if v not in vs]
        kept[(gene, moi)] = dict(g, vars=vs)
    return kept, excluded


def cluster_hit(pos: int, positions: list[int], size=C.CLUSTER_SIZE, window=C.CLUSTER_WINDOW_BP):
    """Largest number of calls in any `window`-bp interval that contains `pos`, if >= size."""
    best = max(sum(s <= p < s + window for p in positions) for s in range(pos - window + 1, pos + 1))
    return best if best >= size else 0


class ClusterCheck:
    """Clustered low-allele-balance calls around candidate variants, read from the VCF."""

    def __init__(self, path: Path):
        self.vcf = VCF(str(path))
        self.chr_prefix = any(s.startswith("chr") for s in self.vcf.seqnames)

    def _records(self, chrom: str, start: int, end: int):
        chrom = chrom if self.chr_prefix else chrom.removeprefix("chr")
        return [r for r in self.vcf(f"{chrom}:{start}-{end}")
                if start <= r.POS <= end and any(a > 0 for a in r.genotypes[0][:2])]

    def check(self, v: Variant) -> dict | None:
        """Details if `v` sits in a cluster and its own allele balance is below the cutoff."""
        chrom, pos, ref, alt = v
        w = C.CLUSTER_WINDOW_BP
        recs = self._records(chrom, pos - w + 1, pos + w - 1)
        n = cluster_hit(pos, [r.POS for r in recs])
        own = next((r for r in recs if r.POS == pos and r.REF == ref and alt in r.ALT), None)
        if not n or own is None:
            return None
        ad = own.format("AD")[0].tolist()
        ab = ad[own.ALT.index(alt) + 1] / sum(ad) if sum(ad) else 0.0
        if ab >= C.CLUSTER_MAX_AB:
            return None
        calls = ", ".join(f"{r.POS} {r.REF}>{','.join(r.ALT)} AD {','.join(map(str, r.format('AD')[0]))}"
                          for r in recs)
        return {"variant": "-".join(map(str, v)), "rule": "cluster",
                "detail": f"{n} calls within {w} bp, AB {ab:.3f}; calls {pos - w + 1}-{pos + w - 1}: {calls}"}


def ploidy_hit(v: Variant, call: Call, sex=C.SAMPLE_SEX) -> dict | None:
    """Details if `v` is a heterozygous call on non-PAR chrX in an XY sample."""
    chrom, pos = v[0].removeprefix("chr"), v[1]
    if sex != "XY" or chrom != "X" or call.zygosity != "het":
        return None
    if any(start <= pos <= end for start, end in C.CHRX_PAR):
        return None
    return {"variant": "-".join(map(str, v)), "rule": "ploidy",
            "detail": f"GT {call.gt} on non-PAR chrX in an {sex} sample"}


def load_read_review(path: Path) -> dict[Variant, dict]:
    """Read-level review decisions keyed by variant; a `fail` needs a written reason."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pos, change = r["variant"].split(" ")
            chrom, pos = pos.split(":")
            ref, alt = change.split(">")
            if r["decision"] not in ("pass", "fail"):
                raise ValueError(f"{r['variant']}: decision must be pass or fail")
            if r["decision"] == "fail" and not r["reason"].strip():
                raise ValueError(f"{r['variant']}: a failed read review needs a written reason")
            out[(chrom, int(pos), ref, alt)] = r
    return out


def choose_variants(vs: dict[Variant, Call], moi: str) -> list[Variant]:
    """Best single variant (dominant), or a hom/hemi or the two best hets (recessive)."""
    ranked = sorted(vs, key=lambda v: (vs[v].key(), v))
    if not ranked:
        return []
    if moi not in RECESSIVE:
        return ranked[:1]
    if vs[ranked[0]].zygosity != "het":
        return ranked[:1]
    hets = [v for v in ranked if vs[v].zygosity == "het"]
    if len(hets) >= 2:
        return hets[:2]
    homs = [v for v in ranked if vs[v].zygosity != "het"]
    return homs[:1]  # empty: recessive gene needs two hits


def build_candidates(groups, panel):
    """One hypothesis per gene (best-scoring MOI), highest Exomiser combined score first."""
    cands = []
    for (gene, moi), g in groups.items():
        chosen = choose_variants(g["vars"], moi)
        if not chosen:
            continue
        p = panel.get(gene)
        c = Candidate(gene, moi, g["score"], g["pheno"], chosen, [g["vars"][v] for v in chosen],
                      in_panel=p is not None, panel_sources=p["sources"] if p else "")
        if c.distance is not None and c.distance < C.CIS_DISTANCE_BP:
            c.notes.append(f"hets {c.distance} bp apart: likely one allele in cis")
        cands.append(c)
    best = {}
    for c in sorted(cands, key=lambda c: (-c.score, c.gene, c.moi)):
        best.setdefault(c.gene, c)
    return sorted(best.values(), key=lambda c: (-c.score, c.gene))


def hedge_pair(top: Candidate, all_rare) -> Candidate | None:
    """Pair the top recessive gene's strongest allele with its rarest other het.

    Hedge for the case where the true second allele is non-coding (e.g. a hypomorphic
    intronic/regulatory variant) rather than the one Exomiser chose. Ties on MAX_FREQ go
    to the higher Exomiser variant score.
    """
    if top.moi not in RECESSIVE or len(top.variants) != 2:
        return None
    others = [(v, call) for v, call in all_rare.get(top.gene, {}).items()
              if v not in top.variants and call.zygosity == "het"]
    if not others:
        return None
    v, call = min(others, key=lambda vc: (vc[1].freq, -vc[1].score, vc[0]))
    c = Candidate(top.gene, top.moi, top.score, top.pheno, [top.variants[0], v], [top.calls[0], call],
                  in_panel=top.in_panel, panel_sources=top.panel_sources, finding_type="hedge")
    c.notes.append("alternative second allele (hedge)")
    return c


def _acmg(call: Call, reviewed: dict) -> tuple[str, str]:
    r = reviewed.get(call.hgvs)
    return (r["class"], "manual") if r else (call.acmg, "automatic")


CANDIDATE_COLS = ["candidate_id", "gene", "finding_type", "inheritance", "hgvs_1", "zygosity_1",
                  "acmg_1", "hgvs_2", "acmg_2", "distance_bp", "phenotype_score", "combined_score",
                  "in_panel",
                  "artefact", "read_fail", "adjust", "adjust_reason",
                  "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2", "ref_2", "alt_2",
                  "panel_sources", "notes"]


def candidate_row(c: Candidate, reviewed: dict) -> dict:
    acmg = [_acmg(call, reviewed) for call in c.calls]
    v2 = c.variants[1] if len(c.variants) > 1 else ("", "", "", "")
    hgvs = " + ".join(f"{call.hgvs} [not classified]" if cls == "NA" else f"{call.hgvs} [{cls}, {src}]"
                      for call, (cls, src) in zip(c.calls, acmg))
    return {
        "candidate_id": f"{c.gene} hedge" if c.finding_type == "hedge" else c.gene,
        "gene": c.gene, "finding_type": c.finding_type, "inheritance": c.moi,
        "hgvs_1": c.calls[0].hgvs, "zygosity_1": c.calls[0].zygosity, "acmg_1": acmg[0][0],
        "hgvs_2": c.calls[1].hgvs if len(c.calls) > 1 else "",
        "acmg_2": acmg[1][0] if len(acmg) > 1 else "",
        "distance_bp": c.distance if c.distance is not None else "",
        "phenotype_score": f"{c.pheno:.4f}", "combined_score": f"{c.score:.4f}",
        "in_panel": int(c.in_panel),
        "artefact": int(c.artefact), "read_fail": int(c.read_fail), "adjust": "", "adjust_reason": "",
        "chrom_1": c.variants[0][0], "pos_1": c.variants[0][1], "ref_1": c.variants[0][2],
        "alt_1": c.variants[0][3], "chrom_2": v2[0], "pos_2": v2[1], "ref_2": v2[2], "alt_2": v2[3],
        "panel_sources": c.panel_sources,
        "notes": "; ".join([f"{c.gene} {c.moi}; {hgvs}", *c.notes]),
    }


def _write(path: Path, rows: list[dict], cols: list[str]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("exomiser", type=Path)
    p.add_argument("panel", type=Path)
    p.add_argument("--out", type=Path, default=C.WORK / "rank")
    p.add_argument("--top", type=int, default=20, help="primary candidates passed to epcr.py")
    p.add_argument("--secondary", default="", help="comma-separated genes reported as secondary")
    p.add_argument("--review", type=Path, default=C.ROOT / "config" / "acmg_review.tsv")
    p.add_argument("--vcf", type=Path, default=C.NORM_VCF if C.NORM_VCF.exists() else C.RAW_VCF,
                   help="VCF for the cluster check (normalised VCF if present, else the raw VCF)")
    p.add_argument("--read-review", type=Path, default=C.ROOT / "config" / "read_review.tsv")
    a = p.parse_args(argv)

    secondary = {g for g in a.secondary.split(",") if g}
    groups, all_rare = load_exomiser(a.exomiser)
    panel = load_panel(a.panel)
    unfiltered = build_candidates(groups, panel)
    groups, rarity = rarity_filter(groups, exempt=secondary)
    cands = build_candidates(groups, panel)
    reviewed = review_classes(a.review) if a.review.exists() else {}

    a.out.mkdir(parents=True, exist_ok=True)
    raw = [dict(candidate_row(c, reviewed), rank=i) for i, c in enumerate(cands, 1)]
    _write(a.out / "raw_ranking.tsv", raw, ["rank", *CANDIDATE_COLS])

    chosen = cands[:a.top]
    chosen += [c for c in cands[a.top:] if c.gene in secondary]
    for c in chosen:
        if c.gene in secondary:
            c.finding_type = "secondary"
            reason = C.SECONDARY_REASONS.get(c.gene)
            c.notes.append(f"secondary finding; {reason}" if reason else "secondary finding")
    hedge = hedge_pair(cands[0], all_rare)

    # Rarity exclusions are listed for genes that would otherwise reach the top-N score range.
    floor = unfiltered[min(a.top, len(unfiltered)) - 1].score
    exclusions = [e for e in rarity if e["score"] >= floor]
    checker = ClusterCheck(a.vcf)
    read_review = load_read_review(a.read_review) if a.read_review.exists() else {}
    for c in [*chosen, *([hedge] if hedge else [])]:
        for v, call in zip(c.variants, c.calls):
            hit = checker.check(v)
            if hit:
                c.cluster = True
                c.notes.append(f"clustered low-AB call ({hit['detail'].split(';')[0]}): likely mapping artefact")
                exclusions.append(dict(hit, gene=c.gene, moi=c.moi, hgvs=call.hgvs, score=c.score))
            hit = ploidy_hit(v, call)
            if hit:
                c.ploidy = True
                c.notes.append(f"het call on non-PAR chrX in an {C.SAMPLE_SEX} sample: likely mapping artefact")
                exclusions.append(dict(hit, gene=c.gene, moi=c.moi, hgvs=call.hgvs, score=c.score))
            rec, idx = find_record(checker.vcf, *v)
            flags, ab = genotype_flags(rec, idx) if rec else ([], None)
            if flags:
                shown = [f"AB {ab:.2f}" if f == "allele_balance" else f for f in flags]
                c.notes.append(f"genotype QC flag ({', '.join(shown)}), low confidence")
            review = read_review.get(v)
            if review and review["decision"] == "fail":
                c.read_fail = True
                c.notes.append(f"failed read review: {review['reason']}")
                exclusions.append({"variant": "-".join(map(str, v)), "rule": "read_review",
                                   "detail": review["reason"], "gene": c.gene, "moi": c.moi,
                                   "hgvs": call.hgvs, "score": c.score})
    _write(a.out / "exclusions.tsv", exclusions, ["gene", "moi", "variant", "hgvs", "rule", "score", "detail"])
    rows = [candidate_row(c, reviewed) for c in [chosen[0], *([hedge] if hedge else []), *chosen[1:]]]
    _write(a.out / "candidates.tsv", rows, CANDIDATE_COLS)
    for r in rows:
        print(r["candidate_id"], r["finding_type"], r["phenotype_score"], r["acmg_1"], r["acmg_2"],
              sep="\t")


if __name__ == "__main__":
    main()
