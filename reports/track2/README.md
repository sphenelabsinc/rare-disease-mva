# Track 2 working files

The plan is in `docs/04-track2-plan.md`. The submitted report
(`sphenelabs_track2_report.md`/`.pdf`) and the Stage 1 dossier
(`01_variant_dossier.md`) contain patient-specific variant detail, so they are
kept local and gitignored; see "Data policy" in the top-level README.

| Stage | Write-up | Script | Inputs / outputs |
|---|---|---|---|
| 0 Literature | — | `scripts/track2/s0_pubmed.py` | `tables/s0_pubmed_hits.tsv` (592 PMIDs, fixed queries) |
| 1 Variant dossier | `01_variant_dossier.md` | `s1_structure.py`, `s1_conservation.py` | `tables/s1_*.tsv` |
| 2 Causal network | `02_causal_network.md` | `s2_network.py` | `config/track2/causal_edges.tsv` → `tables/s2_*.tsv`, `s2_network.mmd` |
| 3 Drug longlist | `03_longlist.md` | `s3_drugs.py` | `config/track2/intervention_targets.tsv` → `tables/s3_longlist.tsv` |
| 4 Ranking | `04_ranking.md` | `s4_rank.py` | `config/track2/candidates.json` → `tables/s4_*.tsv` |
| 5 Validation | `05_validation_plan.md` | — | — |
| 6 Limitations | `06_limitations.md` | — | self-critique; method limits and criticisms of each result |

Run the stages in order with `.venv/bin/python -I scripts/track2/<script>`.
Downloads from public APIs are cached in `data/work/track2/cache/`, which is
gitignored. No patient-level data is written under `reports/track2/`.
