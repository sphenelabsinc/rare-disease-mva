# Track 2 report: mechanism-first, scenario-aware drug repurposing for BUB1B-related MVA1

> Draft, 2026-10-18. Rename to `<hf-username>_track2_report.md` before
> submitting. These are hypotheses for follow-up, **not** treatment
> recommendations. Working files: `reports/track2/01–05_*.md`, `tables/`,
> `scripts/track2/` and `config/track2/`.

## 1. Executive summary

PROBAND01 is a 46,XY child whose features are:

- rhabdomyosarcoma
- intrauterine growth restriction (IUGR) and small for gestational age
- short stature and failure to thrive
- nephrocalcinosis
- a family history of recurrent pregnancy loss

The child carries compound-heterozygous **BUB1B** variants: p.(Leu737\*), which is
ClinVar Pathogenic, and **p.(Asn1002Lys)**, a novel VUS in the C-terminal
pseudokinase domain. Together they are consistent with MVA1.

**Mechanism.** Allele 1 is expected to be null (NMD), so **all of the child's
BUBR1 is the N1002K protein**. We judge the most likely mechanism (prior 0.50)
to be a **destabilised, HSP90-dependent pseudokinase fold that is turned over
quickly**, leaving low BUBR1, a weak spindle checkpoint and mosaic
aneuploidy. The evidence:

- Asn1002 is buried (relative accessibility 0.16, pLDDT 91).
- It is conserved in 117 of 120 vertebrate orthologs and in the BUB1 paralog.
  Lys is never observed at this position.
- It lies 12 Å from L1012, a known MVA variant.
- In patient cells, every analogous buried kinase-domain MVA missense variant
  tested (I909T, R727C, L1012P) has about twice the normal turnover and
  depends on HSP90 (Suijkerbuijk 2010).

The alternatives are a folded protein with a defective pseudokinase/KARD
function (0.25) and a near-null protein (0.10). They are kept open, and the
drug rankings are computed separately under each.

**Candidates.** We ran a reproducible pipeline (variant → causal network →
approved-drug longlist → gated, scenario-weighted ranking with stress tests).
It prioritises:

1. **Niacin (nicotinic acid).** This is the approved NAD+ precursor. It targets
   the NAD+ → SIRT2 → BUBR1 stabilisation axis, which is the only
   demonstrated in vivo route for raising BUBR1. It ranks first in three of
   the four mechanism scenarios.
2. **Low-dose rapalog (sirolimus/everolimus).** It targets mTORC1-linked
   tissue decline, and has a track record of being combined with RMS
   chemotherapy.
3. **Dasatinib ± quercetin.** A senolytic, aimed at progeroid tissue decline.

Alternates: **arimoclomol**, an approved heat-shock co-inducer, which we
promote if the destabilised-fold model is confirmed; and hydroxychloroquine.

**Avoid list:** HSP90 inhibitors, SAC/Aurora/PLK1 inhibitors and p53
inhibition. Microtubule agents and radiotherapy need oncology dose review.

**Honest limits.** Each top candidate rests on a single key mouse study (the
leave-one-source-out test flags all three). None has drug evidence in
BUB1B-mutant human cells. Phase is inferred. A staged validation plan with
go/no-go criteria is in §8.

## 2. Clinical and genomic evidence

| | Allele 1 | Allele 2 |
|---|---|---|
| Variant (GRCh38) | chr15:40209701 T>G | chr15:40220612 T>G |
| HGVS (ENST00000287598) | c.2210T>G p.(Leu737\*) | c.3006T>G p.(Asn1002Lys) |
| Class | Nonsense; ClinVar P/LP for MVA1 | Missense, last exon; ACMG VUS |
| Population | gnomAD 115/1.46M, 0 homozygotes | gnomAD 1/1.46M |
| In silico | — | AlphaMissense 0.92, REVEL 0.47 |
| Reads | 25/46 alt, both strands | 13/28 alt, both strands |

- **Prioritisation:** Exomiser ranks BUB1B first genome-wide. CEP57 and TRIP13
  are clean.
- **Copy number:** no SV or CNV in the MVA genes.
- **Aneuploidy:** bulk WGS shows no whole-chromosome aneuploidy above about
  10% of cells. This doesn't exclude MVA.
- **Phase:** the variants are 10.9 kb apart, so phase can't be observed. Trans
  is **inferred**.
- **Phenotype fit:** the features fit MVA1. Embryonal RMS was in the original
  BUB1B-MVA families (Hanks 2004). Recurrent pregnancy loss is plausible from
  carrier-level PCS (Schmid 2014; Wei 2026).

## 3. Variant-specific mechanism (Stage 1; details in `01_variant_dossier.md`)

**Structure** (AlphaFold v6, O60566; `s1_structure.py`):

- N1002 is buried in the C-lobe of the pseudokinase domain (relative SASA
  0.16, pLDDT 91.1).
- It contacts residues 973–999.
- It is 12.3 Å from L1012 (MVA1 L1012P) and 10.7 Å from T1008, a PLK1
  phosphosite.
- Every known MVA1 missense position in this domain (R814, L844, I909, Q921,
  L1012) is also buried.

**Conservation** (`s1_conservation.py`):

- 97% identity across 120 one-to-one vertebrate orthologs. The only
  substitutes are D, S and E, never K.
- The residue is part of a `W-N-x-F-F-x-x-x-L-N` core shared with BUB1.

**Literature** (Suijkerbuijk 2010, full text):

- MVA missense variants fall into two classes. The *low-abundance* class lies
  in or near the kinase domain: I909T, R727C, L1012P and L844F.
- Turnover is about 2× normal.
- Geldanamycin (an HSP90 inhibitor) plus cycloheximide removed "virtually all
  mutant BUBR1", and MG132 prevented this.
- Adding BUBR1 back restored the checkpoint.
- BUBR1's pseudokinase fold depends on its ATP-site residues for
  conformational stability (Suijkerbuijk 2012). It also supports KARD
  phosphorylation and PP2A-B56 recruitment (Gama Braga 2020).
- Checkpoint function is steeply dosage-dependent: PCS appears with any loss,
  and aneuploidy below about 50% (Bohers 2008).

**Classification.** This is partial loss of function through dosage
insufficiency. It is not gain of function and not dominant-negative.

| Model | Prior |
|---|---|
| M1: unstable fold, low abundance | 0.50 |
| M2: normal abundance, checkpoint-specific defect | 0.05 |
| M3: folded but defective KARD/alignment support | 0.25 |
| M4: near-null | 0.10 |
| M5: not causal / in cis | 0.10 |

**Caution.** The mouse "BubR1 L1002P" in Sieben 2020 is the mouse equivalent
of human L1012P. It is not the patient's N1002.

## 4. Causal pathway analysis (Stage 2; `02_causal_network.md`)

We built a 35-edge curated network: 27 edges demonstrated, 6 inferred and 2
speculative, each with a PMID and an evidence tier. STRING (40 partners) and
Reactome (10 pathways) add context.

```
PRIMARY    p.N1002K -> unstable fold (HSP90 client; proteasome) ┐
           p.L737*  -> NMD                                       ├-> low BUBR1 -> weak SAC, PCS,
           NAD+ -> SIRT2 -| K668-Ac (CBP) -> BUBR1 stability    ┘   missegregation (+ KARD/PP2A-B56
                                                                     alignment defect if M3)
SECONDARY  missegregation -> micronuclei -> cGAS-STING (direction unclear)
           aneuploidy -> p53/p21 (net protective) | p16 senescence -> SASP | apoptosis
           mTORC1 hyperactivity (correlative)
CLINICAL   growth failure/IUGR · embryonal tumours (RMS) · pregnancy loss · treatment hypersensitivity
```

**Key implications:**

1. The only levers on the primary defect that *raise* BUBR1 are NAD+/SIRT2
   and the chaperone/proteasome balance.
2. Every lever on the mitotic machinery itself (MPS1, Aurora B, PLK1, APC/C)
   goes the wrong way.
3. Higher BUBR1 is tumour-protective in mice, even with oncogenic Ras (Baker
   2013 *NCB*). Branch A therefore passes the cancer-safety gate *in
   principle*.
4. Downstream levers can trade tissue health against tumour suppression: for
   example, senescence is a tumour barrier and p53 is protective.

## 5. Computational drug-discovery methods (Stage 3)

**Longlist** (`s3_drugs.py`):

- **Targets:** 26 curated intervention targets
  (`config/track2/intervention_targets.tsv`), each with a desired direction
  (increase / decrease / avoid), plus BUB1B's 40 STRING partners.
- **Sources:** ChEMBL curated mechanisms, restricted to max_phase 4, and
  cross-checked against Open Targets known drugs.
- **Direction check:** each drug–target pair is labelled direction-matched,
  mismatched, AVOID (the drug acts against a protective node) or context-only.
- **Result:** 78 approved drug–target rows (66 drugs). 49 are direction-matched, mostly proteasome, senolytic and JAK/IFN drugs. 19 are AVOID (tubulin agents, including vincristine). **None hit a 'raise BUBR1' target** (SIRT2, NAMPT, HSP90/HSF1, BUB1B): target-based databases can't see substrate-type or network-level drugs. Details in `03_longlist.md`.
- **Additions outside these databases:** we added approved drugs whose
  mechanism isn't a single protein target (niacin, arimoclomol,
  phenylbutyrate, hydroxychloroquine). We also searched ChEMBL for any
  approved drug with measured BUBR1 binding (the pharmacological-chaperone
  idea). There are none at relevant potency: dasatinib leaves >75% activity
  at 1 µM, and ceritinib has an IC50 of 13 µM.

**Ranking** (`s4_rank.py`):

- 100-point rubric, scored per mechanism scenario.
- Weighted by the Stage 1 priors.
- Gates are applied before ranking.
- Stress tests: leave-one-source-out, speculative-edge ablation, and 5,000
  Dirichlet draws over the priors.

**Provenance and use of AI.** Literature retrieval is a fixed set of PubMed
queries (`s0_pubmed.py`, 592 records). Every claim in the evidence files cites
a PMID or a named label or trial. An LLM assistant helped with retrieval,
reading and drafting. Every score is written down in
`config/track2/candidates.json` with its reason, so it can be audited.

## 6. Ranked candidates and counter-evidence (Stage 4; `04_ranking.md`)

| Rank | Candidate (approval) | Hypothesis | Key evidence (tier) | Main counter-evidence / safety |
|---|---|---|---|---|
| 1 | **Niacin** (approved Rx) | NAD+ repletion → SIRT2 deacetylates K668 → more stable BUBR1-N1002K → better SAC | NMN raises BubR1 in vivo; SIRT2 extends BubR1 H/H lifespan (E3, North 2014). Niacin raises human NAD+ up to 8× (Pirinen 2020). More BUBR1 is tumour-protective (E3, Baker 2013) | No direct test on BUBR1. Needs N1002K to work once stabilised. Uric acid and kidney stones (nephrocalcinosis). NAD+ may fuel tumour metabolism. SIRT2 doesn't extend WT lifespan (Wu 2023) |
| 2 | **Sirolimus/everolimus** (approved, paediatric) | Damp mTORC1-driven tissue decline and SASP | mTORC1 hyperactivity tracks sarcopenia in BubR1 MVA mice (E3, correlative). Rapalog combined with RMS chemotherapy was feasible (ARST1431) | No EFS benefit in RMS. Immunosuppression. **Autophagy induction may degrade BUBR1** (37094517) |
| 3 | **Dasatinib ± quercetin** (dasatinib approved, paediatric) | Remove p16+ senescent cells that drive progeroid decline | INK-ATTAC clearance in BubR1 H/H (E3, genetic). D+Q in BubR1 H/H fat (patent only) | Healthspan, not lifespan. Senescence is a tumour barrier. Myelosuppression |
| 4 | Hydroxychloroquine (approved) | Block autophagic BUBR1 loss; aneuploidy-selective stress | Chloroquine is aneuploidy-selective (E4, Tang 2011) | MVA mutants are cleared by the proteasome, not autophagy. Killing aneuploid cells in a mosaic child is untested |
| 5 | Arimoclomol (approved 2024, age ≥2) | Raise chaperone capacity for the HSP90-dependent N1002K fold | HSP90 dependence of analogous MVA mutants (E2) | No BUBR1 data. OCT2-mediated creatinine rise confounds kidney monitoring. HSF1 is pro-tumour in some models |
| — | Bortezomib | (probe) proteasome block stabilises the mutant | MG132 rescue (E2) | **Excluded:** neuropathy with vincristine; chronic toxicity |

**Mechanism dependence:**

- Niacin and arimoclomol are worth most if M1 holds. Arimoclomol ranks 3rd
  under M1 and 6th under M4.
- Rapalogs and senolytics don't depend on the mechanism, but they only treat
  consequences.

**Avoid / flag:**

- **HSP90 inhibitors:** they remove mutant BUBR1, which overrides 17-AAG's
  "aneuploidy-selective" appeal.
- **MPS1, Aurora B and PLK1 inhibitors.**
- **p53 inhibition.**
- **Vincristine and radiotherapy:** standard RMS care, but flag for dose
  review. A PCS-MVA RMS case was treated with reduced-intensity chemotherapy
  (Nishitani-Isa 2019; full text not reviewed). Genomic-instability disorders
  show treatment hypersensitivity (Nakano 2024).

## 7. Potential impact

- **For this child:**
  - A testable Branch A hypothesis that a cheap, approved drug could raise
    the level of the only BUBR1 protein the child has.
  - An explicit avoid list relevant to the child's oncology care.
  - A cell-based assay panel (BUBR1 level, half-life, PCS) that can act as a
    pharmacodynamic biomarker.
- **For other MVA patients:** the same pattern (one null allele plus an
  unstable kinase-domain missense: I909T, R727C, L1012P, L844F) covers a
  large share of reported biallelic MVA1. So the NAD+/chaperone hypothesis
  would apply to that whole subgroup.

## 8. Validation and safety (Stage 5; `05_validation_plan.md`)

- **Tier 0:** confirm the genotype.
  - Parental Sanger testing for phase.
  - Karyotype, PCS and FISH.
  - RNA, to confirm NMD of allele 1 and normal splicing of allele 2.
- **Tier 1:** resolve the mechanism.
  - BUBR1 western blot.
  - Cycloheximide chase with geldanamycin and MG132.
  - KARD/B56 immunofluorescence.
  - Nocodazole and live imaging.
  - Isogenic N1002K/– cells with L1012P/– controls.
  - A matched-expression rescue experiment.
- **Gate:** if N1002K works once its abundance is restored, go with Branch A.
- **Tier 2:** drug testing in patient cells, with vincristine-combination
  safety checks and aneuploid-fraction monitoring.
- **Tier 3:** multidisciplinary review before any clinical consideration.

## 9. Limitations, scalability and conclusions

**Limitations.** A full self-critique, including the criticisms we expect a
reviewer to raise against each result, is in `06_limitations.md`. The ones
that most constrain the conclusions:

- *Phase is inferred, not measured* (the variants are 10.9 kb apart). If they
  are in *cis*, the premise fails. M5 is carried as a caveat rather than
  weighted into the expected scores, so those scores are optimistic by ~10%.
- *N1002K has no functional data.* M1 rests on analogy to other buried MVA
  missense variants (I909T, R727C, L1012P) — and Suijkerbuijk 2010 shows that
  other variants fall into different classes. The mechanism priors are
  elicited judgements; the Dirichlet sampling tests sensitivity around them,
  not whether they are centred correctly.
- *The rank-1 candidate was inserted by hand.* No approved drug hits a Branch
  A target, so niacin, arimoclomol, phenylbutyrate and HCQ were added by
  mechanism reasoning. That is our methodological finding, but it leaves the
  shortlist without a recall denominator, and the key step is not reproducible
  from the scripts.
- *The rubric is one unblinded rater, scored once.* No inter-rater agreement.
  The reported intervals reflect prior uncertainty only, not scoring error,
  which is likely larger and would reorder ranks 3–7.
- *The drug-level evidence is mouse-only and single-source.* All three leads
  are flagged by leave-one-source-out, none has E1/E2 evidence, and the best
  score is 62/100. The ranking is an order of hypotheses to test, not a
  shortlist to act on.
- *The cancer-safety argument is weaker than it reads.* Baker 2013 concerns
  overexpressed **wild-type** BubR1, whereas this would stabilise a mutant of
  unknown checkpoint competence. If N1002K is abundant but defective (M2,
  prior 0.05 — our least defensible number), raising its level could worsen
  missegregation.
- *Sirolimus ranks 2nd on a trial that missed its endpoint* (ARST1431:
  feasible, no EFS benefit), and rapalog-induced autophagy may lower BUBR1
  (37094517) — it could help Branch B while harming Branch A.
- *Dasatinib's D+Q data in this model comes from a patent*, not peer review,
  and senolysis removes a tumour barrier in a cancer-predisposed child.
- *No pharmacokinetic, dosing or tissue-exposure modelling*, and
  nephrocalcinosis — which constrains both Branch A leads — is absent from the
  network.
- *Bulk WGS can't detect low-level mosaicism.*
- *Two load-bearing sources were read only through secondary channels*
  (Pirinen 2020; Nishitani-Isa 2019, for which no abstract is available).
- *Publication bias runs toward our conclusion*: negative intervention results
  in this mouse model are largely unpublished.

The most robust single result is a negative one: HSP90 inhibitors should be
avoided, because they strip mutant BUBR1 (demonstrated experimentally in MVA
patient variants). A naive "aneuploidy-selective drug" pipeline would instead
have returned 17-AAG as a hit.

**Scalability.** The pipeline is gene-agnostic. Its parts are:

- structure and conservation scripts that take a UniProt ID;
- fixed literature queries;
- a curated edge schema with tier and status;
- target-direction lists;
- ChEMBL/Open Targets retrieval;
- a JSON rubric with scenario priors and gates.

It applies directly to CEP57-, TRIP13- and MAD1L1-related aneuploidy
syndromes, and more generally to any missense VUS whose mechanism is
uncertain. The scenario-aware ranking makes explicit how drug priorities
depend on what the variant actually does.

**Conclusion.** The clearest evidence chain runs from an unstable,
HSP90-dependent N1002K pseudokinase fold to low BUBR1 to a weak checkpoint.
That points to raising BUBR1 abundance as the disease-modifying strategy, with
niacin as the approved candidate to test first. One assay (a BUBR1 half-life
and rescue experiment in patient cells) decides whether this branch survives.

## Use of AI assistance

Analysis was run with Claude Code (Anthropic, Claude Opus 5.5).
**[Fill in: plan/tier and data-handling setting, e.g. "Claude <tier>,
training opt-out enabled, no feedback provided".]**

## Key references

- Hanks 2004 *Nat Genet* 15475955
- Matsuura 2006 *AJMG A* 16411201
- Bohers 2008 *Hum Genet* 18932004
- Suijkerbuijk 2010 *Cancer Res* 20516114
- Suijkerbuijk 2012 *Dev Cell* 22698286 and 23079597
- Xu 2013 *Biol Open* 23789096
- Gama Braga 2020 *Cell Rep* 33207204
- Matsumura 2007 *JBC* 17376779
- Baker 2004 *Nat Genet* 15208629
- Baker 2011 *Nature* 22048312
- Baker 2013 *Cell Rep* 23602569
- Baker 2013 *Nat Cell Biol* 23242215
- Wijshake 2012 *PLoS Genet* 23300461
- North 2014 *EMBO J* 24825348
- Wu 2023 *Aging Cell* 38009412
- Sieben 2020 *JCI* 31738183
- Goutas 2023 *Redox Biol* 37094517
- Tang 2011 *Cell* 21315436
- Cohen-Sharir 2021 *Nature* 33505028
- Mackenzie 2017 *Nature* 28738408
- Yost 2017 *Nat Genet* 28553959
- Malumbres 2024 *Nat Rev Genet* 39169218
- Nakano 2024 *Clin Cancer Res* 39264246
- Nishitani-Isa 2019 *Pediatr Int* 31184400
- Ochiai 2019 *Pediatr Blood Cancer* 31081598
- Schmid 2014 *Chromosome Res* 24981203
- Wei 2026 *Front Endocrinol* 42434306
- Pirinen 2020 *Cell Metab* (niacin and NAD+)
- ARST1431 (Gupta 2024 *Lancet Oncol*)
- Miplyffa (arimoclomol) US label
