# Track 2 plan: MVA1 drug repurposing

The requirements come from `03-drug-repurposing.md` and the case details from
`reports/track1_methods.md`. Revision 2 (2026-10-18) brings in reviewer feedback:

- The second allele's mechanism is now something to test, not something assumed.
- The causal hierarchy separates the primary defect from what follows from it.
- Every candidate gets a cancer-safety check.
- Evidence is graded into tiers.
- Drug rankings depend on which mechanism scenario is true.

## Case summary

The proband is a 46,XY child with:

- rhabdomyosarcoma
- intrauterine growth restriction, small for gestational age
- short stature and failure to thrive
- nephrocalcinosis
- a family history of recurrent pregnancy loss

Proband-only WGS found compound-heterozygous _BUB1B_ variants (MANE
ENST00000287598). These fit mosaic variegated aneuploidy syndrome 1 (MVA1,
OMIM 257300).

- **c.2210T>G p.(Leu737\*)**: nonsense, upstream of the C-terminal domain.
  ClinVar Pathogenic for MVA1; ACMG Pathogenic. Expected to be null via NMD.
  NMD hasn't been tested.
- **c.3006T>G p.(Asn1002Lys)**: novel missense in the C-terminal
  pseudokinase domain, in the last exon (so no NMD is expected). gnomAD 1 in
  about 1.46M alleles, AlphaMissense 0.92, REVEL 0.47; ACMG VUS. It lies 10
  residues from the known MVA missense L1012P. Its functional effect is
  **unknown**: "hypomorphic" is a hypothesis.

Bulk WGS shows no whole-chromosome aneuploidy above about 10% of cells, which
doesn't exclude MVA. Phase is inferred (the variants are 10.9 kb apart and no
parents were sequenced).

## Guiding principles

1. **Resolve the mechanism before choosing treatments.** Drug priorities depend
   on what p.Asn1002Lys does, so rank candidates separately under each
   plausible mechanism.
2. **Keep the primary defect separate from what follows it.** Restoring
   accurate chromosome segregation is a different goal from reducing
   senescence or inflammation.
3. **Cancer safety is a gate, not a weight.** The child already has
   rhabdomyosarcoma. Anything that helps genomically unstable cells survive, or
   interacts with RMS treatment (vincristine and other microtubule agents in a
   checkpoint-deficient background), is excluded or sent for specialist review
   whatever its score. Kidney safety (nephrocalcinosis) is a gate in the same
   way.
4. **Approved means approved.** Keep approved pharmaceuticals separate from
   experimental compounds and supplements (e.g. NAD precursors). Record the
   jurisdiction and any paediatric labelling.
5. **Falsifiability.** Every hypothesis states what result would refute it.
6. **Provenance.** Every evidence edge cites a primary source. An LLM can
   suggest papers and pull out claims, but it never assigns uncited evidence
   and never decides the ranking.

## Causal hierarchy

```
BUB1B biallelic dysfunction        (working hypothesis; phase and allele 2 function unconfirmed)
  -> reduced functional BUBR1 and/or defective kinetochore signalling
       primary: SAC weakening, KARD / PP2A-B56 recruitment, alignment errors,
                abnormal mitotic timing, chromosome missegregation, PCS
  -> secondary, context-dependent responses
       p53, cGAS-STING, proteotoxic stress, mTORC1, senescence (p16Ink4a)
  -> clinical manifestations
       growth failure, tissue dysfunction, cancer predisposition (RMS, Wilms, leukaemia)
```

Links are not strictly linear, and the secondary responses can either protect
or harm depending on context. Label every edge as **demonstrated** (in human or
BUB1B models), **inferred** (from a related model) or **speculative**.

## Competing mechanisms for p.Asn1002Lys

| ID | Model | Predicted observation | Drug implication |
|---|---|---|---|
| M1 | Reduced protein stability | Low BUBR1 abundance | Stabilise or raise BUBR1 (e.g. SIRT2/NAD axis) |
| M2 | Protein present, checkpoint function impaired | Near-normal abundance, abnormal mitosis | Raising levels may not be enough; target the cellular phenotype |
| M3 | Impaired kinetochore localisation or KARD signalling | Less BUBR1 or PP2A-B56 at kinetochores, alignment defects | Rescue the interaction or the downstream phosphatase balance |
| M4 | Near-null | Very little functional BUBR1 | Only downstream mitigation is realistic |
| M5 | Genotype doesn't explain the case | Functional tests negative, or variants in cis | Revisit phase and alternative diagnoses |

Structural question to test: does Asn1002Lys destabilise the pseudokinase
fold, which supports KARD phosphorylation and PP2A-B56 recruitment? Test it
against AlphaFold or experimental structures, conservation, and the location
of L1012P and other MVA missense variants.

## Evidence tiers

| Tier | Description | Example |
|---|---|---|
| E1 | Rescue in a model carrying this patient's variant | Drug reduces missegregation in p.Asn1002Lys cells |
| E2 | Rescue in another BUB1B-deficient human cell model | Improved checkpoint phenotype in BUB1B-mutant cells |
| E3 | Rescue in a closely related animal model | BubR1-hypomorphic or MVA-allele mice |
| E4 | Mechanistic evidence in another setting | Drug modulates a relevant node in an unrelated system |
| E5 | Computational only | Network proximity, signature reversal |

There is no E1 evidence yet, and the report must not imply there is. Signature
reversal alone (E5) never puts a drug in the top tier.

## Workflow

### Stage 1: Variant evidence dossier (decision gate)

- Re-check transcript numbering, read-level quality, splice effects
  (SpliceAI on both variants), NMD assumptions and the phase inference.
- Analyse p.Asn1002Lys structurally: AlphaFold or PDB models of the BUBR1
  C-terminal domain, conservation, its distance to L1012P and other known MVA
  missense variants, and nearby interaction surfaces.
- Literature: MVA1 genotype–phenotype; allele-specific mouse models;
  pseudokinase-domain function.
- **Output:** a dossier plus the M1–M5 matrix, with a prior weight for each
  model.
- **Gate:** these weights set the scenario weights in Stage 4. If M5 can't be
  ruled out, say so up front in the report.

### Stage 2: Evidence-annotated causal network

- Nodes come from Reactome (mitotic spindle checkpoint, kinetochore
  signalling), STRING neighbours of BUB1B, and curated primary papers.
- Each edge records its source, its tier and whether it is demonstrated,
  inferred or speculative.
- Primary-defect, secondary-response, cancer and organ-specific modules are
  kept separate.
- **Output:** a machine-readable graph (TSV or JSON with a fixed schema) and a
  mechanism figure.

### Stage 3: Reproducible approved-drug longlist

- Query DrugBank/ChEMBL, Open Targets and DGIdb against the network nodes,
  then confirm approval status from regulatory labels (FDA, EMA).
- Perturbation signatures (CMap/LINCS) count as E5 support only.
- Retrieval runs from versioned dataset snapshots, so it is deterministic.
- **Output:** a longlist listing targets, mechanism, approval and provenance,
  with approved, experimental and supplement entries in separate tables.

### Stage 4: Scenario-aware ranking and stress tests

Score each candidate under each mechanism scenario M1–M4, using this rubric
(out of 100):

| Criterion | Weight |
|---|---|
| Direct BUB1B/MVA functional evidence | 25 |
| Mechanistic fit to the variant scenario | 20 |
| Relevant downstream phenotype evidence | 15 |
| Paediatric safety and regulatory suitability | 15 |
| Oncology and renal compatibility | 15 |
| Experimental and clinical feasibility | 10 |

- Safety gates (principle 3) are applied before scoring.
- The **robust candidates** are those that stay highly ranked across the
  scenarios, weighted by the Stage 1 priors.
- Sensitivity checks: leave one source out, and drop the speculative edges.
  A drug whose rank depends on a single paper or a single speculative edge is
  flagged.
- **Output:** the top 3–5 candidates. Each has at least one primary source,
  at least one counterargument or safety limitation, a falsifiable hypothesis
  and a measurable endpoint.

### Stage 5: Validation proposal (proposal only)

- Baseline assays in patient lymphocytes or fibroblasts: BUBR1 abundance
  (western blot), kinetochore BUBR1 and PP2A-B56 (immunofluorescence),
  mitotic timing (live imaging), missegregation, premature chromatid
  separation and karyotype or interphase FISH.
- Isogenic CRISPR knock-in of p.Asn1002Lys in a near-diploid line or iPSCs,
  with an L1012P control.
- Drug-response endpoints tied to each candidate's hypothesis, with go/no-go
  criteria.
- Parental Sanger testing to confirm the variants are on opposite copies.
- **Output:** a staged experimental plan that distinguishes M1–M5.

### Stage 6: Package the submission

- Report as `<hf-username>_track2_report.md`/`.pdf`.
- GitHub repo with scripts, configs and dataset versions. No patient-level
  data, per the data policy.
- Methods description, including the AI disclosure: provider, plan or tier,
  and data-handling setting.
- 3-minute pitch video.

## Report outline (about 17 pages plus references)

| Section | Pages |
|---|---|
| Executive summary and top hypotheses | 1 |
| Clinical and genomic evidence | 2 |
| Variant-specific mechanisms (M1–M5) | 3 |
| Causal pathway analysis | 2 |
| Computational discovery methods | 2 |
| Ranked candidates and counterevidence | 4 |
| Validation and safety | 2 |
| Limitations, scalability, conclusions | 1 |

Figures: (1) a BUBR1 domain map with the variants and known MVA missense
variants; (2) the evidence-labelled causal network; (3) the discovery and
ranking workflow; (4) a candidate matrix of hypothesis, tier, scenario
robustness and safety constraints.

## Initial candidate classes (leads, not recommendations)

- **SIRT2/NAD axis.** SIRT2 stabilises BubR1 and improves phenotypes in mice
  (E3). NAD precursors are supplements, not approved drugs. Check whether any
  approved drug acts on this axis. *Strongest mechanism-first lead, valid
  under M1.*
- **mTOR inhibitors (sirolimus, everolimus).** There is an mTORC1 rationale
  in some BubR1-mutant mice. No rescue of chromosome segregation has been
  shown. Concerns are immunosuppression, growth, interactions with RMS
  treatment, and kidney effects. *Exploratory.*
- **Senolytics.** In BubR1-hypomorphic mice, clearing senescent cells delayed
  some ageing phenotypes but did not extend lifespan. Navitoclax is not
  approved; dasatinib is approved for other indications. The cancer-safety
  analysis is critical here. *Exploratory.*
- **cGAS-STING and other inflammatory modulators.** The therapeutic
  direction isn't clear and depends on context (acute vs adapted chromosomal
  instability). *Exploratory.*
- **Avoid list.** Agents that weaken the checkpoint or cause aneuploidy
  (MPS1 or Aurora inhibitors). Flag microtubule poisons used in RMS treatment
  as an interaction concern.

## Starting papers (verify each citation before use)

- Hanks et al. 2004, *Nat Genet*: BUB1B mutations in MVA, with embryonal RMS.
- Matsuura et al. 2006, *Am J Med Genet A*: MVA with BUB1B variants and
  premature chromatid separation in Japanese families.
- Baker et al. 2004, *Nat Genet*: BubR1-hypomorphic mice develop progeroid
  features.
- Baker et al. 2011, *Nature*: clearing senescent cells in BubR1 H/H mice
  delays ageing phenotypes.
- Baker et al. 2013, *Nat Cell Biol*: BubR1 overexpression extends healthy
  lifespan.
- North et al. 2014, *EMBO J* (PMID 24825348; corrected from *Nat Commun*):
  SIRT2 stabilises BubR1 and extends lifespan.
- Sieben et al. 2020, *J Clin Invest*: BubR1 allele combinations, including
  L1012P, drive phenotypic heterogeneity in MVA mice.
- Gama Braga et al. 2020, *Cell Rep* (PMID 33207204, verified): the domain supports
  KARD phosphorylation and PP2A-B56 recruitment.
- Malumbres et al. 2024, *Nat Rev Genet* (PMID 39169218, verified): review of aneuploidy in
  development, ageing and cancer.

## Research questions, in priority order

1. What does p.Asn1002Lys do to BUBR1, and how confident can we be?
2. Which molecular defects could an existing approved drug plausibly reverse?
3. Which drugs have credible experimental evidence (E2/E3 or better) on those
   defects?
4. Could any benefit be offset by tumour risk, kidney toxicity or interference
   with RMS treatment?
5. Which one or two hypotheses could a lab realistically test first?

## SciSpace prompt

```
I'm building a mechanism-first drug-repurposing hypothesis for mosaic variegated aneuploidy syndrome type 1 (MVA1, OMIM 257300).

Case: a 46,XY child with rhabdomyosarcoma, IUGR/SGA, short stature, failure to thrive and nephrocalcinosis, with familial recurrent pregnancy loss. Compound-heterozygous BUB1B (ENST00000287598):
- p.(Leu737*): nonsense, ClinVar Pathogenic.
- p.(Asn1002Lys): novel missense in the C-terminal pseudokinase domain, 10 residues from the known MVA variant L1012P. AlphaMissense 0.92, REVEL 0.47; function unknown.
Phase is inferred, not observed.

Find and summarise primary literature to answer these questions, in priority order:

1. Variant mechanism. Compare these competing models for a C-terminal BUBR1 missense: (M1) reduced protein stability; (M2) normal abundance but impaired checkpoint function; (M3) impaired kinetochore localisation or KARD phosphorylation / PP2A-B56 recruitment; (M4) near-null. Look for functional studies of pseudokinase-domain missense variants (especially L1012P, and the Sieben 2020 JCI allele-combination mouse models), structural work on the BUBR1 pseudokinase fold, and MVA genotype–phenotype series.
2. Causal hierarchy. Separate what is demonstrated, inferred or speculative for: primary defects (SAC weakening, alignment errors, missegregation, PCS); secondary responses (p53, cGAS-STING, proteotoxic stress, mTORC1, senescence); and clinical outcomes (growth failure, RMS/Wilms/leukaemia predisposition).
3. Interventions with experimental evidence in BUB1B/BubR1-deficient systems: SIRT2/NAD-mediated BubR1 stabilisation, BubR1 overexpression, senescent-cell clearance, mTOR inhibition. For each, report whether it improved chromosome segregation or only downstream phenotypes, its effect on tumour incidence, and any lifespan effect.
4. Cancer and kidney safety. Evidence on whether these interventions, or suppression of stress responses, allow aneuploid or tumour cells to persist. Also look at sensitivity to microtubule agents (vincristine) and radiation in checkpoint-deficient cells or MVA patients.
5. Approved drugs only. For any candidate, confirm whether it is an approved pharmaceutical (vs experimental compound or supplement), and its paediatric and renal safety.

For each paper, give: citation, year, study type (patient, cell, mouse or review), evidence tier (E2 = BUB1B-deficient human cells; E3 = related animal model; E4 = mechanistic, other setting; E5 = computational), key finding, which mechanism model or intervention it bears on, and any evidence against it. Flag contested findings and gaps. Mark the 10–15 must-read papers.
```
