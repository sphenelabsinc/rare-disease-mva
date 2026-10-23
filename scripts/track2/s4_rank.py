"""Stage 4: scenario-aware ranking with safety gates and stress tests.

Reads config/track2/candidates.json. For each candidate:

  score(M) = c1 + c2[M] + c3 + c4 + c5 + c6          (out of 100)
  expected = sum over M of prior(M) * score(M)

Gates are applied before ranking: "exclude" removes a candidate and "review"
keeps it but marks it as needing specialist review. Stress tests:

  * per-scenario ranks (M1..M4) and robustness (worst rank across scenarios);
  * leave-one-source-out: drop every c1/c3 item from one source and record
    the largest rank change;
  * speculative-edge ablation: halve c2 for candidates whose mechanistic fit
    rests on a speculative edge;
  * prior uncertainty: 5,000 Dirichlet draws around the scenario priors
    (concentration 10), reporting P(top 3) and a 90% interval on expected score.

Writes reports/track2/tables/s4_ranking.tsv and s4_stress.tsv.

Usage: python -I scripts/track2/s4_rank.py
"""

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CFG = ROOT / "config/track2/candidates.json"
OUT = ROOT / "reports/track2/tables"
SCEN = ["M1", "M2", "M3", "M4"]
CAPS = dict(c1=25, c2=20, c3=15, c4=15, c5=15, c6=10)


def scores(c, drop_src=None, halve_spec=False):
    c1 = min(CAPS["c1"], sum(e["pts"] for e in c["c1"] if e["src"] != drop_src))
    c3 = min(CAPS["c3"], sum(e["pts"] for e in c["c3"] if e["src"] != drop_src))
    fixed = c1 + c3 + c["c4"]["pts"] + c["c5"]["pts"] + c["c6"]["pts"]
    out = {}
    for m in SCEN:
        c2 = c["c2"][m] / 2 if (halve_spec and c["c2"]["speculative"]) else c["c2"][m]
        out[m] = fixed + min(CAPS["c2"], c2)
    return out


def expected(s, priors):
    w = np.array([priors[m] for m in SCEN])
    w = w / w.sum()
    return float(sum(w[i] * s[m] for i, m in enumerate(SCEN)))


def rank(values):
    order = sorted(values, key=lambda k: -values[k])
    return {k: i + 1 for i, k in enumerate(order)}


def main():
    cfg = json.loads(CFG.read_text())
    priors = cfg["scenario_priors"]
    cands = [c for c in cfg["candidates"] if c["gate"] != "exclude"]
    excluded = [c["drug"] for c in cfg["candidates"] if c["gate"] == "exclude"]
    for c in cfg["candidates"]:
        for k in ("c4", "c5", "c6"):
            assert 0 <= c[k]["pts"] <= CAPS[k], (c["drug"], k)
        assert all(0 <= c["c2"][m] <= CAPS["c2"] for m in SCEN), c["drug"]

    base = {c["drug"]: scores(c) for c in cands}
    exp = {d: expected(s, priors) for d, s in base.items()}
    r_exp = rank(exp)
    r_scen = {m: rank({d: s[m] for d, s in base.items()}) for m in SCEN}

    # Leave-one-source-out, over each candidate's own sources
    loso = {d: (0, "") for d in base}
    for c in cands:
        d = c["drug"]
        for src in sorted({e["src"] for e in c["c1"] + c["c3"]}):
            r = rank({x["drug"]: expected(scores(x, drop_src=src), priors) for x in cands})
            delta = abs(r[d] - r_exp[d])
            if delta > loso[d][0]:
                loso[d] = (delta, src)

    # Speculative-edge ablation
    r_spec = rank({c["drug"]: expected(scores(c, halve_spec=True), priors) for c in cands})

    # Dirichlet prior uncertainty
    rng = np.random.default_rng(20261018)
    alpha = 10 * np.array([priors[m] for m in SCEN]) / sum(priors.values())
    draws = rng.dirichlet(alpha, size=5000)
    mat = np.array([[base[d][m] for m in SCEN] for d in base])  # cands x scen
    ev = draws @ mat.T                                            # draws x cands
    order = np.argsort(-ev, axis=1)
    names = list(base)
    p_top3 = {names[j]: float(np.mean((order[:, :3] == j).any(axis=1))) for j in range(len(names))}
    lo, hi = np.percentile(ev, [5, 95], axis=0)

    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for c in sorted(cands, key=lambda c: r_exp[c["drug"]]):
        d = c["drug"]
        j = names.index(d)
        worst = max(r_scen[m][d] for m in SCEN)
        rows.append(dict(
            rank=r_exp[d], drug=d, branch=c["branch"], gate=c["gate"],
            expected_score=round(exp[d], 1), ci90=f"{lo[j]:.1f}-{hi[j]:.1f}",
            **{f"score_{m}": base[d][m] for m in SCEN},
            **{f"rank_{m}": r_scen[m][d] for m in SCEN},
            worst_rank=worst, p_top3=round(p_top3[d], 2),
            hypothesis=c["hypothesis"], endpoint=c["endpoint"], counter=c["counter"]))
    with open(OUT / "s4_ranking.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    with open(OUT / "s4_stress.tsv", "w") as fh:
        fh.write("drug\tbase_rank\tmax_rank_shift_loso\tmost_influential_source\trank_speculative_halved\tsingle_source_flag\n")
        for r in rows:
            d = r["drug"]
            shift, src = loso[d]
            fh.write(f"{d}\t{r['rank']}\t{shift}\t{src}\t{r_spec[d]}\t{'FLAG' if shift >= 2 else ''}\n")

    print(f"Excluded by gate: {excluded}")
    hdr = f"{'rank':>4} {'drug':42s} {'gate':7s} {'E[score]':>8} {'90%':>11} " + " ".join(f"{m:>4}" for m in SCEN) + "  worst P(top3) LOSO spec"
    print(hdr)
    for r in rows:
        d = r["drug"]
        print(f"{r['rank']:>4} {d[:42]:42s} {r['gate']:7s} {r['expected_score']:>8} {r['ci90']:>11} "
              + " ".join(f"{r[f'rank_{m}']:>4}" for m in SCEN)
              + f"  {r['worst_rank']:>5} {r['p_top3']:>7} {loso[d][0]:>4} {r_spec[d]:>4}")


if __name__ == "__main__":
    main()
