# Stage 6: Limitations and criticisms

A reviewer's-eye read of our own work. Written to be adversarial: where the
method is weak, and where the results do not support the weight that could be
put on them. Nothing here is rhetorical hedging — each item names what would
have to change for the conclusion to move.

---

## 1. Limitations of the method, by stage

### 1.1 Variant interpretation (Stage 1)

- **Phase is inferred, not measured.** The two variants are 10.9 kb apart,
  which is beyond short-read phasing. We assumed *trans* because the phenotype
  fits MVA1 and because recurrent pregnancy loss in the family is consistent
  with both parents carrying one allele. If they are in *cis*, the entire
  Track 2 premise is wrong (scenario M5, prior 0.10). We carried M5 as a
  caveat rather than as a weighted column in the ranking, which understates it:
  under M5 **every** candidate scores zero, so the honest expected score of
  the top candidate is ~0.9 × 57.8 ≈ 52, not 57.8. We did not propagate this.
- **Structural inference is from a predicted model.** AlphaFold pLDDT 91.1 at
  N1002 is high, but there is no experimental structure of the BUBR1
  pseudokinase domain with this region resolved at side-chain accuracy. RSA
  0.16 ("buried") and the measured 10.7 Å to the T1008 PLK1 site are
  properties of a model, not of a protein. A 2–3 Å error in side-chain
  placement would not change "buried", but it would change any claim about
  specific contacts — so we made none.
- **Conservation is suggestive, not mechanistic.** 97% conservation across 120
  orthologs with D/S/E as the only substitutions tells us the position is
  under selection. It does not tell us *which* function fails. Asn→Lys is a
  charge-introducing substitution at a buried position, which is why we
  weighted instability (M1) over a pure checkpoint defect (M2) — but that is
  a biophysical prior, not data on this variant.
- **The ortholog set is biased.** Compara one-to-one orthologs skew to
  well-annotated genomes. A set that over-represents mammals inflates apparent
  conservation. We did not weight by phylogenetic distance.
- **Mechanism priors (M1 0.50, M3 0.25, …) are elicited, not estimated.** They
  come from our reading of Suijkerbuijk 2010 plus the position's burial. A
  different analyst could justify M1 0.30 / M3 0.45 from the same papers. The
  Dirichlet sampling tests sensitivity *around* our priors; it does not test
  whether the priors are centred correctly. That is the single largest
  unquantified uncertainty in the pipeline.
- **Guilt by analogy.** The instability evidence is from I909T, R727C and
  L1012P — other MVA missense variants, not N1002K. We are asserting that
  N1002K behaves like its structural neighbours. Suijkerbuijk themselves show
  this is not safe: R550Q, Q921H and Y155C fall into *different* classes. One
  Tier 1 western blot could overturn M1.
- **A near-miss worth recording:** the mouse allele written "L1002P" in
  Sieben 2020 is human L1012P, not this patient's residue 1002. We caught it;
  a pipeline that string-matched positions would not have.

### 1.2 Causal network (Stage 2)

- **The network is hand-curated, so it inherits our reading.** 35 edges, 27
  demonstrated / 6 inferred / 2 speculative. Curation gives provenance (every
  edge has a PMID and tier) but guarantees selection bias: we added the edges
  we went looking for. There is no systematic recall estimate — we cannot say
  what fraction of the real BUB1B interventional network we captured.
- **Direction and context are flattened.** An edge "NAD+ → SIRT2 → BUBR1" is
  annotated as a single direction of effect, but these relationships are
  tissue-, dose- and age-dependent. Most supporting work is in aged mouse
  tissue or cancer lines, not in a growing child's skeletal muscle or kidney.
- **STRING and Reactome contribute association, not causation.** The 40 STRING
  partners (≥700) are co-mention and co-expression driven; they informed
  target *enumeration* only, and none produced an approved drug anyway.
- **Speculative-edge ablation is weaker than it looks.** Halving c2 for
  candidates resting on speculative edges changed no ranks. That is partly
  because the top three rest on demonstrated edges — but also because c2 is
  capped at 20 of 100, so halving it can move a score by at most 10. The
  ablation is under-powered by construction.

### 1.3 Drug longlist (Stage 3)

- **The target-based search has a structural blind spot, and we patched it by
  hand.** No approved drug has a curated mechanism on any Branch A
  "raise BUBR1" target. We therefore added niacin, arimoclomol,
  phenylbutyrate and hydroxychloroquine by mechanism reasoning. This is our
  headline methodological finding, but it is also the pipeline's weakest
  seam: **the eventual rank-1 candidate was inserted by hand, not discovered
  by the search.** Reproducibility of the automated part does not extend to
  that decision.
- **Hand-insertion has no denominator.** We cannot state a recall. There may
  be other approved substrate-type or network-level drugs we simply did not
  think of; nothing in the method would reveal them. A systematic alternative
  (transcriptomic signature reversal, e.g. LINCS/CMap, or a chemical-genetics
  screen) would give a denominator. We did not run one.
- **max_phase 4 is a coarse filter.** It takes no account of paediatric
  licensing, formulation, route, or whether approval is for an indication that
  makes chronic use plausible. We applied those judgements later, by hand,
  in the rubric.
- **Annotation routing errors are likely under-detected.** We found two
  (rapalogs → FKBP12, metformin → complex I) and fixed one by adding FKBP1A;
  everolimus and temsirolimus still had to be added by name. If two surfaced
  by accident, others probably remain.
- **ChEMBL activity screening for a BUBR1 pharmacological chaperone was
  shallow.** 15 compounds with measured activity, none approved at relevant
  potency. Absence of recorded activity is not absence of activity — BUBR1 is
  not a routine screening target, so the database is close to empty by
  default. The TG101209 observation (BUB1B protein 2.2× in CAKI-1, IC50
  67 µM) is a single cell-line result at a concentration with no clinical
  relevance.
- **No pharmacokinetic filter.** Nothing in the pipeline asks whether a drug
  reaches the tissues that matter at an achievable dose. For an intracellular
  mitotic protein in muscle, bone and kidney, that is a serious omission.

### 1.4 Ranking (Stage 4)

- **The rubric is a documented opinion with arithmetic on top.** The six
  weights (25/20/15/15/15/10) were chosen by us, before scoring, which
  protects against tuning to a desired winner — but they are not derived from
  anything. Shifting weight from "direct BUB1B/MVA evidence" (25) toward
  "paediatric suitability" (15) would favour niacin further; shifting it the
  other way would favour nothing, because no candidate has direct evidence.
- **Point allocations are not blinded or duplicated.** One analyst assigned
  every number. There is no second rater and therefore no inter-rater
  agreement statistic. In a rubric-based method, that is the standard
  reliability check and we do not have it.
- **The precision is false.** "57.8 (55.0–60.1)" invites a reading as a
  measurement. The interval reflects only prior uncertainty over M1–M4 with
  concentration 10 — a parameter we also chose. It does **not** include
  scoring error, which is almost certainly larger. A ±5-point rubric
  disagreement would reorder ranks 3–7.
- **Leave-one-source-out measures fragility, not correctness.** It tells us
  niacin drops 2 ranks without North 2014. It cannot tell us whether North
  2014 is right.
- **Scores are not comparable across branches.** Branch A candidates are
  scored on "would raise BUBR1"; Branch B/C candidates on "would mitigate
  downstream phenotype". These are different claims with different
  falsifiability, summed on one scale. A high Branch B score is a weaker
  statement than the same Branch A score.
- **Gates are binary and applied before ranking**, so a "review" candidate and
  a "pass" candidate sit in the same list. Sirolimus (immunosuppression in a
  child with a cancer predisposition) and niacin (flushing, uric acid) are not
  comparable risks, and the single 15-point safety column does not express
  that gap.

### 1.5 Validation plan (Stage 5)

- **Tier 0/1 feasibility is assumed, not confirmed.** The plan needs patient
  fibroblasts or lymphoblasts and parental samples. We do not know these exist
  or are consented for this use. If they do not, every go/no-go gate is
  unreachable and the ranking stays untested.
- **Cell-based readouts may not reflect the affected tissues.** A fibroblast
  western blot for BUBR1 level answers M1, but IUGR, short stature,
  nephrocalcinosis and RMS risk are not fibroblast phenotypes.
- **No power calculation.** "Reduces PCS" and "raises BUBR1" are stated as
  directions without effect sizes, replicate counts or the variability of
  these assays. The gates are not yet decidable.

### 1.6 Cross-cutting

- **Literature access.** Two load-bearing sources were read only through
  secondary channels: Pirinen 2020 (niacin raising human NAD+ up to 8×) and
  Nishitani-Isa 2019 (reduced-intensity chemotherapy in MVA-associated RMS).
  Both are cited in support of our conclusions.
- **A patent is cited as data.** The dasatinib + quercetin result in BubR1 H/H
  mice comes from US11517572, which is not peer-reviewed and whose assignee
  has an interest in the result.
- **One analyst, one pass, LLM-assisted.** No independent replication of any
  stage. Scripts are deterministic given the caches, which makes the
  *computation* reproducible; it does not make the *judgements* reproducible.
- **Publication bias runs toward our conclusion.** "Intervention X rescues
  BubR1 mouse" is publishable; "X does nothing" largely is not. Every Branch A
  and B lead rests on positive-result mouse literature.

---

## 2. Where the results are open to criticism

### 2.1 The top three each rest on one mouse paper

This is the central criticism and it is in the report. Niacin → North 2014
(*EMBO J*); sirolimus → Sieben 2020 (*JCI*); dasatinib → Baker 2011
(*Nature*). None has E1 or E2 evidence — nothing in BUB1B-mutant human cells,
nothing in this variant. Best score is 62/100 under M1, and the rubric's
ceiling is unreachable by design because the direct-evidence column is empty
for every candidate. **The correct reading of the ranking is "ordered
hypotheses for Tier 1 testing", not "a treatment shortlist".**

### 2.2 Niacin's evidence chain has a species and model gap at every link

The chain is: NAD+ declines → SIRT2 deacetylates BubR1 K668 → BubR1 is
stabilised → phenotype improves. North 2014 shows this with **NMN in aged
mice**. Our candidate is **niacin in a child with a missense BUB1B allele**.
Three substitutions at once: different compound, different species, different
mechanism of low BubR1 (age-related decline vs. a destabilising mutation).
Pirinen 2020 bridges only the pharmacology of raising NAD+ in humans, not the
BubR1 link. And Wu 2023 (SIRT2 overexpression does not extend WT lifespan)
cuts against the pathway's potency. A reviewer would be right to say niacin
ranks first partly because it is cheap, approved and low-risk — i.e. on
rubric columns 4–6 — rather than on mechanistic evidence.

### 2.3 The safety reasoning for the cancer gate is thinner than it reads

We used Baker 2013 (*NCB*: BubR1 overexpression is tumour-protective) to argue
that raising BUBR1 is safe in a child with RMS. That paper is about
**overexpressing wild-type** BubR1 in mice. Our candidate stabilises a
**mutant** protein whose checkpoint competence is unknown. If N1002K is
abundant but checkpoint-defective (M2, or a partial M3), then raising its
level could increase missegregation, not reduce it — the opposite of the
intended effect. We assigned M2 a prior of only 0.05, which makes this risk
look negligible; that prior is the least defensible number in the matrix,
since nothing we measured excludes normal abundance with a functional defect.

### 2.4 Sirolimus is ranked on a trial that failed its endpoint

ARST1431 shows temsirolimus + chemotherapy is *feasible* in rhabdomyosarcoma,
with **no EFS benefit**. We scored that as safety/feasibility evidence, which
is legitimate, but the rank-2 position can read as efficacy support. Separately,
there is a direct counterargument we logged and did not resolve: rapalogs
induce autophagy, and autophagy may *lower* BUBR1 (PMID 37094517). Sirolimus
could work on Branch B while actively harming Branch A. Adding chronic
immunosuppression for a child with a cancer-predisposition syndrome is also a
larger risk than a 15-point column conveys.

### 2.5 Dasatinib's supporting data is partly non-peer-reviewed

Baker 2011 supports senolysis in BubR1 H/H mice; the specific D+Q result in
that model is from a patent. Senescence is also a tumour-suppressive barrier,
so clearing senescent cells in a cancer-predisposed child is a real hazard,
and our own falsification criterion ("falsified if it enriches aneuploid
survivors") describes a mechanism by which the intervention could do harm. A
reviewer could reasonably argue dasatinib should have been gated to `exclude`
alongside bortezomib rather than `review`.

### 2.6 Arimoclomol's rank is the most informative number, and it is unstable

It sits 5th overall but 3rd under M1, and 6th under M4 — it moves on exactly
the question the Tier 1 experiment answers. That is the pipeline working as
designed. But it also means the overall ordering is partly an artefact of
priors we elicited rather than measured, and the 8% P(top 3) is a statement
about our Dirichlet concentration parameter as much as about biology. It also
raises creatinine 10–20% via OCT2, which matters with nephrocalcinosis;
"review" understates that.

### 2.7 The HSP90 conclusion is the most robust result here

Worth stating, since little else is: HSP90 inhibitors should be avoided
because they strip mutant BUBR1 (Suijkerbuijk 2010, geldanamycin +
cycloheximide, MG132-reversible). This is a negative, mechanism-level,
experimentally demonstrated result in MVA patient variants, and it overrides
17-AAG's "aneuploidy-selective" appeal (Tang 2011). Negative findings like
this are the output we should have most confidence in — and note that a
naive "aneuploidy-selective drug" pipeline would have returned 17-AAG as a
*hit*.

### 2.8 Scope criticisms a reviewer will raise

- **No dosing, no exposure modelling, no interaction analysis** beyond the
  vincristine flag.
- **No consideration of what is actually actionable** for this family:
  surveillance, nutrition and growth management are not drug repurposing but
  may dominate outcome.
- **Nephrocalcinosis is under-modelled.** It constrains niacin (uric acid) and
  arimoclomol (creatinine) — our two Branch A leads — and renal handling is
  not in the network at all.
- **Phenotype weighting is implicit.** We never stated which phenotype we are
  trying to modify. Cancer risk, growth and the progeroid features have
  different time courses and would give different rankings.

---

## 3. What would most change the conclusions

In order of expected information gained per unit effort:

1. **Parental testing** → confirms or destroys the *trans* assumption (M5).
2. **BUBR1 western blot + half-life in patient cells, ± HSP90 and proteasome
   inhibition** → collapses M1/M2/M4 from a prior into an observation. This is
   one experiment and it is upstream of nearly every criticism above.
3. **A second, blinded rater re-scoring the rubric** → the missing reliability
   check; cheap, no wet lab.
4. **A systematic signature-based search (LINCS/CMap)** → gives the hand-
   insertion step a denominator.
5. **Checkpoint-competence assay on N1002K** → decides whether raising the
   mutant's level is therapeutic or harmful (§2.3).

Nothing here is a treatment recommendation. All candidates are hypotheses
requiring the Tier 0–1 work in `05_validation_plan.md` before any clinical
discussion, and any such discussion belongs with the treating team.
