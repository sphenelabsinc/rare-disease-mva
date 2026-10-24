# Track 2 pitch video script (3 minutes)

Team: sphenelabs · Rare Disease Real Kid: MVA Hackathon 2026, Track 2 (drug
repurposing). Public video: written to be safe for YouTube.

**What this script deliberately leaves out (data policy).** No genomic
coordinates, no exact variant names (HGVS / residue number / substituted amino
acid), no read counts, allele frequencies or other sample metrics, no
secondary or carrier findings, and no clinical detail beyond the coded features
in the organisers' public case summary. The patient's missense allele is
described only as "a novel missense in the C-terminal pseudokinase domain."
Named variants (I909T, R727C, L1012P, L844F) are previously published MVA
variants from other patients, not this child's.

Pace: about 150 words a minute, so roughly 450 spoken words. Timings are
targets; the bracketed lines are on-screen visuals.

---

## 0:00–0:20 · Hook

> [Title card: "One gene, one fold, one assay: a mechanism-first repurposing
> hypothesis for MVA1" · team name · repo URL]

A standard "kill aneuploid cells" repurposing search ranks HSP90 inhibitors
near the top. For this child, that's likely the one class of drug that would
make the underlying defect worse. Here's why, and what we'd test instead.

## 0:20–0:50 · The case and the gene

> [Coded features as icons: rhabdomyosarcoma, growth restriction, short
> stature, nephrocalcinosis, family history of pregnancy loss → arrow → BUB1B]

The child's features fit mosaic variegated aneuploidy syndrome type 1. Our
Track 1 analysis found two variants in BUB1B, the gene for BUBR1, a core
component of the spindle assembly checkpoint. One is a known pathogenic
truncating variant, expected to make no protein. The other is a novel missense
in the C-terminal pseudokinase domain, with no functional data. So all of the
BUBR1 this child has comes from that one missense allele. The whole drug
question turns on what that change does.

## 0:50–1:25 · Mechanism, kept as competing hypotheses

> [BUBR1 domain map: pseudokinase domain, published MVA missense sites marked;
> AlphaFold structure with the region highlighted, no residue label]
> [Table: M1–M5 with priors 0.50 / 0.05 / 0.25 / 0.10 / 0.10]

We didn't assume a mechanism; we weighed five. The residue is buried in the
pseudokinase fold, about twelve ångströms from L1012, a published MVA variant,
and conserved in 117 of 120 vertebrate orthologs. In other patients' cells, every
comparable buried kinase-domain MVA variant tested is turned over about twice
as fast and depends on the chaperone HSP90 (Suijkerbuijk 2010). So our leading
model, at fifty percent, is an unstable fold: low BUBR1, weak checkpoint,
mosaic aneuploidy. We kept a folded-but-defective model and a near-null model
open, and ranked drugs separately under each.

## 1:25–2:00 · Method

> [Workflow figure: PubMed queries → structure + conservation → 35-edge
> causal network → ChEMBL / Open Targets longlist → gated, scenario-weighted
> ranking → stress tests]

The pipeline is reproducible from public data: a 35-edge causal network where
every edge cites a PMID and is labelled demonstrated, inferred or speculative;
ChEMBL and Open Targets for approved drugs, with required directions; cancer
and kidney safety as gates before scoring; then a 100-point rubric scored
under each mechanism and stress-tested.

Our key methods finding: no approved drug in those databases hits any target
that would *raise* BUBR1. Substrate drugs like niacin are invisible to
target-based search. We had to add them by mechanism, and we say so.

## 2:00–2:35 · Candidates and the avoid list

> [Ranking table: niacin 57.8, rapalog 46.7, dasatinib ± quercetin 42.1;
> arimoclomol as M1-dependent alternate; red "avoid" panel]

First: niacin, the approved NAD-plus precursor. NAD-plus feeds SIRT2, which
stabilises BUBR1; in mice that axis raises BubR1 and extends lifespan (North
2014), and more BUBR1 protects against tumours (Baker 2013). It's top-ranked
under three of four mechanisms. Second: a low-dose rapalog, for
mTOR-linked tissue decline. Third: dasatinib with or without quercetin, a
senolytic. If the unstable-fold model is confirmed, arimoclomol, an approved
heat-shock co-inducer, moves up.

Avoid: HSP90 inhibitors, which strip mutant BUBR1; checkpoint kinase
inhibitors; p53 suppression. And flag vincristine and radiotherapy, standard
in rhabdomyosarcoma care, for dose review.

## 2:35–3:00 · Honest limits and how to confirm

> [Staged plan: Tier 0 → Tier 1 gate → Tier 2 drug tests, with go/no-go
> criteria]

Each lead rests on one key mouse study. So the next step is one decisive
experiment: in the child's cells,
measure BUBR1 level and half-life, with and without HSP90 and proteasome
inhibition, and test whether the mutant protein works once its level is
restored. If it does, niacin goes first. If not, we pivot to the downstream
candidates. Everything is open on GitHub. These are hypotheses to test, not
treatments. Thank you.

---

## On-screen text and description box (copy into YouTube)

> A mechanism-first, scenario-aware drug-repurposing hypothesis for BUB1B-related
> mosaic variegated aneuploidy syndrome type 1 (MVA1), submitted to the Rare
> Disease Real Kid: MVA Hackathon 2026, Track 2. These are hypotheses for
> laboratory follow-up, not treatment recommendations. No patient-level data
> is shown.
>
> Code, configs and public-data tables:
> https://github.com/sphenelabsinc/rare-disease-mva (`scripts/track2/`,
> `config/track2/`, `reports/track2/`)
>
> AI use: analysis assisted by Claude Code (Anthropic, Claude Opus 5.5), Claude
> Pro plan, model-training setting off, no feedback submitted.

## Methods detail (presenter notes and backup slides)

**Stage 0, literature.** `s0_pubmed.py` runs a fixed query set against PubMed
(592 records), so the corpus is reproducible. An LLM helped retrieve and read
papers; it never assigned uncited evidence and never set a score.

**Stage 1, variant mechanism.** `s1_structure.py` (AlphaFold DB v6, UniProt
O60566): relative solvent accessibility, pLDDT, residue contacts and distances
to published MVA missense sites. Every published MVA1 missense position in the
pseudokinase domain (R814, L844, I909, Q921, L1012) is also buried.
`s1_conservation.py` (Ensembl Compara): 120 one-to-one vertebrate orthologs
plus the BUB1 paralog. The result is a partial loss-of-function, dosage
insufficiency classification (not gain-of-function, not dominant-negative),
with five mechanism models:

| Model | Description | Prior |
|---|---|---|
| M1 | Unstable fold, low abundance (HSP90-dependent, proteasomal turnover) | 0.50 |
| M2 | Normal abundance, checkpoint-specific defect | 0.05 |
| M3 | Folded, but defective KARD phosphorylation / PP2A-B56 recruitment | 0.25 |
| M4 | Near-null | 0.10 |
| M5 | Not causal, or variants in cis | 0.10 (carried as a caveat) |

**Stage 2, causal network.** `s2_network.py` reads
`config/track2/causal_edges.tsv`: 35 curated edges (27 demonstrated, 6
inferred, 2 speculative), each with a PMID and an evidence tier (E1 = rescue in
this variant, down to E5 = computational only). STRING (40 partners, score ≥
0.7) and Reactome (10 pathways) add context. Modules are kept separate:
primary defect (low BUBR1 → weak checkpoint, premature chromatid separation,
missegregation), secondary responses (p53, cGAS-STING, senescence/SASP,
mTORC1) and clinical outcomes.

**Stage 3, approved-drug longlist.** `s3_drugs.py`: ChEMBL curated mechanisms,
max_phase 4, cross-checked against Open Targets (76 of 78 rows confirmed).
Result: 78 drug–target rows, 66 approved drugs; 49 direction-matched
(proteasome, senolytic, JAK/IFN), 19 on the avoid side (tubulin agents,
including vincristine). Zero hits on SIRT2, NAMPT, HSP90/HSF1 or BUB1B itself.
A ChEMBL search for any compound with measured BUBR1 activity found no
approved drug at relevant potency.

**Stage 4, ranking.** `s4_rank.py` reads `config/track2/candidates.json`, in
which every point has a written reason. Rubric: direct BUB1B/MVA evidence 25,
mechanistic fit to scenario 20, downstream phenotype evidence 15, paediatric
and regulatory suitability 15, oncology and kidney compatibility 15,
feasibility 10. Gates: exclude (bortezomib: neuropathy with vincristine) or
specialist review.

| Rank | Candidate | Branch | E[score] (90% band) | Rank M1/M2/M3/M4 | P(top 3) |
|---|---|---|---|---|---|
| 1 | Niacin | A: raise BUBR1 via NAD+/SIRT2 | 57.8 (55.0–60.1) | 1/1/1/2 | 1.00 |
| 2 | Sirolimus / everolimus | B/C: mTORC1, SASP | 46.7 (46.3–47.2) | 2/2/2/1 | 1.00 |
| 3 | Dasatinib ± quercetin | B: senolytic | 42.1 (41.5–42.8) | 4/3/3/3 | 0.92 |
| 4 | Hydroxychloroquine | A/B: autophagy | 40.0 | 5/4/4/4 | 0.00 |
| 5 | Arimoclomol | A: chaperone capacity | 39.3 (36.4–41.8) | 3/5/5/6 | 0.08 |

Stress tests: leave-one-source-out flags all three leads (niacin rests on
North 2014, rapalog on Sieben 2020, dasatinib on Baker 2011). Removing
speculative edges changes no ranks. The best absolute score is 62/100: these
are hypotheses, not a shortlist to act on.

**Counterarguments we raise ourselves.**
- Niacin: raises uric acid, a concern with nephrocalcinosis; NAD+ may support
  tumour metabolism; SIRT2 does not extend wild-type mouse lifespan (Wu 2023);
  Baker 2013 used *wild-type* BubR1, and raising a defective mutant could
  worsen missegregation if M2 holds.
- Rapalog: the RMS trial ARST1431 showed feasibility but no event-free survival
  benefit; induced autophagy may lower BUBR1 (Goutas 2023).
- Dasatinib ± quercetin: healthspan, not lifespan; senescence is a tumour
  barrier in a cancer-predisposed child; the D+Q data in this mouse model come
  from a patent.

## Future steps and ways to confirm

Proposal only; none of this has been done.

1. **Tier 0, confirm the genotype (weeks).** Parental testing to show the two
   variants are on opposite copies; karyotype with ≥100 metaphases, premature
   chromatid separation scoring and interphase FISH; RNA analysis to confirm
   the truncating allele undergoes NMD and the missense allele splices normally.
2. **Tier 1, resolve the mechanism (months).** Quantitative BUBR1 western blot;
   cycloheximide chase with geldanamycin and MG132 arms (as in Suijkerbuijk
   2010); kinetochore BUBR1 and PP2A-B56 immunofluorescence and KARD
   phosphorylation; nocodazole challenge and live imaging of mitotic timing;
   isogenic knock-in cells with an L1012P control; rescue with wild-type vs
   mutant BUBR1 at matched expression.
   **Gate:** if the mutant works once its level is restored (M1), niacin and
   arimoclomol go first; if not (M3/M4), pivot to the downstream candidates.
3. **Tier 2, drug tests in patient cells.** Clinically achievable doses, vehicle
   and positive controls, and every candidate combined with vincristine to
   check for added mitotic-arrest toxicity. Go criteria: BUBR1 up ≥30% *and*
   chromatid separation or missegregation down (Branch A), or the target
   phenotype down (Branch B); no rise in aneuploid-cell survival; no synergy
   with vincristine toxicity.
4. **Tier 3, translation.** Only after Tier 2, and only with paediatric
   oncology, nephrology, clinical genetics and ethics review. Blood biomarkers
   that could be followed: lymphocyte BUBR1 level, premature chromatid
   separation rate and blood NAD+.

**Scalability.** The pipeline is gene-agnostic: structure and conservation
from a UniProt ID, fixed literature queries, a curated edge schema with tiers,
target-direction lists, ChEMBL/Open Targets retrieval and a JSON rubric with
scenario priors and gates. It applies directly to CEP57-, TRIP13- and
MAD1L1-related aneuploidy syndromes, and to any missense VUS whose mechanism
is uncertain. The same null-plus-unstable-missense pattern covers a large
share of reported biallelic MVA1, so a confirmed niacin/chaperone result would
matter beyond this child.

## Reference materials

PMIDs as recorded in the report; verify against PubMed before publishing the
description.

| Reference | PMID | Used for |
|---|---|---|
| Hanks et al. 2004, *Nat Genet* | 15475955 | BUB1B mutations cause MVA; embryonal RMS in original families |
| Matsuura et al. 2006, *Am J Med Genet A* | 16411201 | MVA with premature chromatid separation |
| Bohers et al. 2008, *Hum Genet* | 18932004 | Checkpoint dosage sensitivity |
| Suijkerbuijk et al. 2010, *Cancer Res* | 20516114 | MVA missense classes; fast turnover, HSP90 dependence |
| Suijkerbuijk et al. 2012, *Dev Cell* | 22698286 | Pseudokinase fold and conformational stability |
| Gama Braga et al. 2020, *Cell Rep* | 33207204 | Pseudokinase domain supports KARD and PP2A-B56 |
| Baker et al. 2004, *Nat Genet* | 15208629 | BubR1-hypomorphic mice are progeroid |
| Baker et al. 2011, *Nature* | 22048312 | Clearing senescent cells in BubR1 H/H mice |
| Baker et al. 2013, *Nat Cell Biol* | 23242215 | BubR1 overexpression: longer healthspan, fewer tumours |
| North et al. 2014, *EMBO J* | 24825348 | NAD+/SIRT2 stabilises BubR1 |
| Wu et al. 2023, *Aging Cell* | 38009412 | SIRT2 and wild-type lifespan (counterevidence) |
| Sieben et al. 2020, *J Clin Invest* | 31738183 | MVA allele-combination mice; mTORC1 |
| Goutas et al. 2023, *Redox Biol* | 37094517 | Autophagy and BUBR1 loss (rapalog counterargument) |
| Tang et al. 2011, *Cell* | 21315436 | Aneuploidy-selective compounds (17-AAG, chloroquine) |
| Malumbres et al. 2024, *Nat Rev Genet* | 39169218 | Aneuploidy review |
| Nakano et al. 2024, *Clin Cancer Res* | 39264246 | Treatment hypersensitivity in instability disorders |
| Pirinen et al. 2020, *Cell Metab* | — | Niacin raises human NAD+ |
| Gupta et al. 2024, *Lancet Oncol* (ARST1431) | — | Temsirolimus with RMS chemotherapy |
| Miplyffa (arimoclomol) US label | — | Approval, age ≥2, kidney monitoring |

## Recording checklist

- [ ] No slide shows coordinates, the HGVS name, the residue number or the
      substituted amino acid (crop AlphaFold screenshots; no residue labels).
- [ ] No screenshot of the submitted report, `01_variant_dossier.md`, Track 1
      files, IGV or any VCF/BAM.
- [ ] Network and ranking figures come from the committed tables only.
- [ ] Closing card repeats "hypotheses for follow-up, not treatment
      recommendations."
- [ ] Upload as unlisted first; check the transcript and captions before making
      it public.
