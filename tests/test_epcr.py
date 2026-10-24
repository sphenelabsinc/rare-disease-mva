import pytest

from rare_disease_mva import config as C
from rare_disease_mva.epcr import (
    exclusion_reason,
    ordering_check,
    phenotype_term,
    review_classes,
    score,
    score_all,
    sensitivity,
    tavtigian_class,
)


def cand(cid, **kw):
    base = {"candidate_id": cid, "gene": cid, "finding_type": "primary", "inheritance": "AD",
            "hgvs_1": "c.1A>G", "zygosity_1": "het", "acmg_1": "VUS", "hgvs_2": "", "acmg_2": "",
            "distance_bp": "", "phenotype_score": "0.5", "in_panel": "0", "artefact": "0",
            "read_fail": "0", "adjust": "", "adjust_reason": ""}
    base.update({k: str(v) for k, v in kw.items()})
    return base


@pytest.mark.parametrize("points,cls", [(10, "P"), (9, "LP"), (6, "LP"), (5, "VUS"), (0, "VUS"),
                                        (-1, "LB"), (-6, "LB"), (-7, "B")])
def test_tavtigian_class(points, cls):
    assert tavtigian_class(points) == cls


def test_formula_uses_weaker_allele_for_compound_het():
    c = cand("pair", inheritance="AR", acmg_1="P", hgvs_2="c.2T>C", acmg_2="VUS",
             distance_bp=10911, phenotype_score=0.982)
    assert score(c) == round(0.6 * 0.982 + 0.4 * C.ACMG_TIER["VUS"], 3) == 0.749


def test_cis_pair_is_penalised():
    far = cand("far", inheritance="AR", acmg_1="LP", hgvs_2="x", acmg_2="LP", distance_bp=500)
    near = cand("near", inheritance="AR", acmg_1="LP", hgvs_2="x", acmg_2="LP", distance_bp=40)
    assert score(near) == round(score(far) * C.CIS_PENALTY, 3)


@pytest.mark.parametrize("kw,reason", [
    ({"acmg_1": "LB"}, "ACMG B/LB"),
    ({"acmg_1": "P", "hgvs_2": "x", "acmg_2": "B"}, "ACMG B/LB"),
    ({"artefact": 1}, "artefact (prone region, clustered low-AB call or het call on non-PAR chrX)"),
    ({"read_fail": 1}, "failed read review"),
    ({"inheritance": "AR"}, "recessive gene with a single het"),
])
def test_hard_filters(kw, reason):
    c = cand("x", **kw)
    assert exclusion_reason(c) == reason
    assert score(c) is None


def test_recessive_homozygote_is_kept():
    assert exclusion_reason(cand("hom", inheritance="AR", zygosity_1="hom")) is None


def test_hedge_and_secondary_are_fixed():
    assert score(cand("h", finding_type="hedge", acmg_2="NA", hgvs_2="x")) == C.HEDGE_EPCR
    assert score(cand("s", finding_type="secondary", phenotype_score=0.99)) == C.SECONDARY_EPCR


def test_manual_adjustment_needs_reason_and_limit():
    with pytest.raises(ValueError):
        score(cand("a", adjust=0.05))
    with pytest.raises(ValueError):
        score(cand("a", adjust=0.2, adjust_reason="too big"))
    assert score(cand("a", adjust=0.05, adjust_reason="published phenotype match")) == round(
        score(cand("a")) + 0.05, 3)


def test_panel_only_counts_in_control_weights():
    c = cand("p", in_panel=1)
    assert score(c) == score(cand("p", in_panel=0))
    assert score(c, C.EPCR_CONTROL) > score(cand("p", in_panel=0), C.EPCR_CONTROL)


def test_ordering_check_rejects_a_tie_at_the_top():
    rows = score_all([cand("a", phenotype_score=0.7), cand("b", phenotype_score=0.7)])
    with pytest.raises(ValueError):
        ordering_check(rows)
    ordering_check(score_all([cand("a", phenotype_score=0.8), cand("b", phenotype_score=0.7)]))


def test_excluded_rows_go_last():
    rows = score_all([cand("lb", acmg_1="LB", phenotype_score=0.99), cand("ok")])
    assert [r["candidate_id"] for r in rows] == ["ok", "lb"]


def test_sensitivity_reports_control_and_stability():
    cands = [cand("top", acmg_1="P", phenotype_score=0.9),
             cand("panel", acmg_1="VUS", phenotype_score=0.6, in_panel=1)]
    sens = sensitivity(cands)
    assert all(s["top"] == "top" and s["top_unchanged"] for s in sens)
    assert [s["role"] for s in sens].count("control (panel bonus)") == 1


def test_phenotype_term_uses_phenotype_not_combined_score():
    # A truncating VUS in a gene with no phenotype match: the combined score is inflated by the
    # variant score, which the ACMG tier already counts.
    c = cand("x", phenotype_score=0.5, combined_score=0.75)
    assert score(c) == round(0.6 * 0.5 + 0.4 * 0.4, 3)
    assert score(c, term="combined") == round(0.6 * 0.75 + 0.4 * 0.4, 3)


@pytest.mark.parametrize("pheno,expected", [(0.5, 0.0), (0.3, 0.0), (0.75, 0.5), (1.0, 1.0)])
def test_baseline_corrected_phenotype_term(pheno, expected):
    assert phenotype_term(cand("x", phenotype_score=pheno), "baseline") == pytest.approx(expected)


def test_sensitivity_includes_previous_definition_and_baseline_arms():
    cands = [cand("top", acmg_1="P", phenotype_score=0.8, combined_score=0.98),
             cand("vus", phenotype_score=0.5, combined_score=0.75)]
    sens = {s["role"]: s for s in sensitivity(cands)}
    assert sens["previous definition (combined score)"]["top_epcr"] == round(0.6 * 0.98 + 0.4, 3)
    assert sens["baseline-corrected phenotype score"]["runner_up_epcr"] == 0.16
    assert all(s["top_unchanged"] for s in sens.values())


def test_sensitivity_skips_combined_arm_without_combined_scores():
    roles = [s["role"] for s in sensitivity([cand("a", acmg_1="P"), cand("b")])]
    assert "previous definition (combined score)" not in roles


def test_review_classes(tmp_path):
    p = tmp_path / "review.tsv"
    p.write_text("variant\thgvs\tcode\tstrength\tpoints\tin_exomiser\tdecision\tstatus\treason\n"
                 "v\tc.1\tPVS1\tvery_strong\t8\tyes\tkeep\tPENDING\tx\n"
                 "v\tc.1\tPP4\tsupporting\t1\tyes\tkeep\tFINAL\tx\n"
                 "v\tc.1\tPM2_Supporting\tsupporting\t1\tyes\tkeep\tFINAL\tx\n"
                 "w\tc.2\tPM2_Supporting\tsupporting\t1\tyes\tkeep\tFINAL\tx\n")
    r = review_classes(p)
    assert r["c.1"]["class"] == "P" and r["c.1"]["pending"] == ["PVS1"]
    assert r["c.2"]["class"] == "VUS" and r["c.2"]["pending"] == []
