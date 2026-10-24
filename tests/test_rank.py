import pytest

from rare_disease_mva.epcr import submission_rows
from rare_disease_mva.rank import (
    ARTEFACT_PRONE,
    Call,
    Candidate,
    _hgvs,
    build_candidates,
    candidate_row,
    choose_variants,
    cluster_hit,
    hedge_pair,
    load_read_review,
    ploidy_hit,
    rarity_filter,
)


def call(gt="0/1", score=0.5, freq=0.0, hgvs="c.1A>G", acmg="VUS"):
    return Call(gt, hgvs, acmg, score, freq)


@pytest.mark.parametrize("gene,artefact", [("OR8U1", True), ("OR2T4", True), ("HLA-DRB1", True),
                                           ("ORC1", False), ("ORC6", False), ("ORAI1", False),
                                           ("BUB1B", False)])
def test_artefact_regex_keeps_orc_and_orai(gene, artefact):
    assert bool(ARTEFACT_PRONE.match(gene)) is artefact


def test_hgvs_matches_review_table_form():
    assert _hgvs("BUB1B:ENST00000287598.11:c.2210T>G:p.(Leu737*)") == "c.2210T>G p.(Leu737*)"
    assert _hgvs("BUB1B:ENST00000287598.11:c.2679-1026A>G:p.(=)") == "c.2679-1026A>G"


def test_recessive_takes_two_best_hets_not_first_by_position():
    vs = {("chr1", 100, "A", "G"): call(score=0.2),
          ("chr1", 200, "C", "T"): call(score=0.9),
          ("chr1", 300, "G", "A"): call(score=0.9, freq=0.01),
          ("chr1", 400, "T", "C"): call(score=0.9, freq=0.5)}
    assert [v[1] for v in choose_variants(vs, "AR")] == [200, 300]


def test_recessive_single_het_has_no_candidate_and_hom_is_kept():
    assert choose_variants({("chr1", 1, "A", "G"): call()}, "AR") == []
    assert choose_variants({("chr1", 1, "A", "G"): call(gt="1/1")}, "AR") == [("chr1", 1, "A", "G")]


def test_dominant_takes_best_variant():
    vs = {("chr1", 1, "A", "G"): call(score=0.3), ("chr1", 2, "A", "G"): call(score=0.8)}
    assert choose_variants(vs, "AD") == [("chr1", 2, "A", "G")]


def test_panel_membership_does_not_change_order():
    groups = {("PANELGENE", "AD"): {"score": 0.7, "pheno": 0.5, "vars": {("chr1", 1, "A", "G"): call()}},
              ("OTHER", "AD"): {"score": 0.8, "pheno": 0.5, "vars": {("chr2", 1, "A", "G"): call()}}}
    cands = build_candidates(groups, {"PANELGENE": {"gene": "PANELGENE", "sources": "x"}})
    assert [c.gene for c in cands] == ["OTHER", "PANELGENE"]
    assert cands[1].in_panel and not cands[0].in_panel


def test_one_hypothesis_per_gene_keeps_best_moi():
    groups = {("G", "AD"): {"score": 0.6, "pheno": 0.5, "vars": {("chr1", 1, "A", "G"): call()}},
              ("G", "AR"): {"score": 0.9, "pheno": 0.5, "vars": {("chr1", 1, "A", "G"): call(gt="1/1")}}}
    assert [(c.gene, c.moi) for c in build_candidates(groups, {})] == [("G", "AR")]


def test_hedge_takes_rarest_other_het_then_best_score():
    a, b = ("chr15", 10, "T", "G"), ("chr15", 20, "T", "G")
    top = Candidate("G", "AR", 0.9, 0.8, [a, b], [call(score=1.0), call(score=0.9)])
    rare = {a: call(), b: call(),
            ("chr15", 5, "A", "G"): call(freq=0.05, score=0.9),
            ("chr15", 30, "A", "G"): call(freq=0.0, score=0.1),
            ("chr15", 40, "A", "G"): call(freq=0.0, score=0.3),
            ("chr15", 50, "A", "G"): call(gt="1/1")}
    h = hedge_pair(top, {"G": rare})
    assert h.variants == [a, ("chr15", 40, "A", "G")] and h.finding_type == "hedge"


def test_submission_keeps_hedge_and_secondary_within_ten_rows():
    rows = [{"candidate_id": f"p{i}", "finding_type": "primary", "epcr": 0.9 - i / 100}
            for i in range(12)]
    rows += [{"candidate_id": "h", "finding_type": "hedge", "epcr": 0.08},
             {"candidate_id": "s", "finding_type": "secondary", "epcr": 0.1}]
    out = submission_rows(rows)
    assert len(out) == 10
    assert [r["candidate_id"] for r in out[-2:]] == ["s", "h"]
    assert out[-1]["finding_type"] == "primary" and out[-2]["finding_type"] == "secondary"


def test_rarity_filter_uses_moi_threshold_and_exempts_secondary():
    groups = {("DOM", "AD"): {"score": 0.8, "pheno": 0.5,
                              "vars": {("chr1", 1, "C", "T"): call(freq=0.77),
                                       ("chr1", 2, "C", "T"): call(freq=0.005)}},
              ("REC", "AR"): {"score": 0.9, "pheno": 0.8,
                              "vars": {("chr2", 1, "T", "G"): call(freq=0.00998),
                                       ("chr2", 2, "T", "G"): call(freq=0.05)}},
              ("SEC", "AD"): {"score": 0.7, "pheno": 0.5, "vars": {("chr3", 1, "C", "G"): call(freq=0.5)}}}
    kept, excluded = rarity_filter(groups, exempt={"SEC"})
    assert list(kept[("DOM", "AD")]["vars"]) == [("chr1", 2, "C", "T")]
    assert len(kept[("REC", "AR")]["vars"]) == 2      # both below the recessive 0.1% threshold
    assert len(kept[("SEC", "AD")]["vars"]) == 1
    assert [(e["gene"], e["rule"]) for e in excluded] == [("DOM", "rarity")]


@pytest.mark.parametrize("pos,positions,expected", [
    (749, [735, 742, 743, 749], 3),          # three calls within 8 bp
    (749, [735, 749], 0),
    (100, [91, 100, 109], 0),                # spread over 19 bp: no 10 bp window holds three
    (100, [100, 101, 109, 120], 3),
])
def test_cluster_hit_counts_calls_in_a_ten_bp_window(pos, positions, expected):
    assert cluster_hit(pos, positions) == expected


@pytest.mark.parametrize("v,gt,sex,flagged", [
    (("chrX", 119469743, "A", "AG"), "0/1", "XY", True),    # non-PAR het in XY
    (("chrX", 119469743, "A", "AG"), "1/1", "XY", False),   # hemizygous call
    (("chrX", 119469743, "A", "AG"), "0/1", "XX", False),
    (("chrX", 1_000_000, "A", "G"), "0/1", "XY", False),    # PAR1
    (("chr1", 1000, "A", "G"), "0/1", "XY", False),
])
def test_ploidy_hit_flags_non_par_chrx_het_in_xy(v, gt, sex, flagged):
    assert bool(ploidy_hit(v, call(gt=gt), sex=sex)) is flagged


def test_read_review_keys_by_variant_and_needs_a_reason_to_fail(tmp_path):
    p = tmp_path / "read_review.tsv"
    p.write_text("variant\tgene\tdecision\treason\nchr5:140807399 GC>G\tPCDHA4\tfail\tno reads\n")
    assert load_read_review(p)[("chr5", 140807399, "GC", "G")]["decision"] == "fail"
    p.write_text("variant\tgene\tdecision\treason\nchr5:140807399 GC>G\tPCDHA4\tfail\t\n")
    with pytest.raises(ValueError):
        load_read_review(p)


def test_candidate_table_carries_phenotype_and_combined_scores():
    c = Candidate("G", "AD", 0.75, 0.5, [("chr1", 1, "A", "G")], [call()])
    row = candidate_row(c, {})
    assert (row["phenotype_score"], row["combined_score"]) == ("0.5000", "0.7500")
    c.cluster = True
    assert candidate_row(c, {})["artefact"] == 1
