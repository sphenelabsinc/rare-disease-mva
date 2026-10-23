# Stage 5: Validation proposal (proposal only)

This is a staged plan for a qualified lab. Each tier has go/no-go criteria,
and earlier tiers change the weights used in later ones. None of this work has
been done.

## Tier 0: Confirm the genotype and resolve M5 (weeks)

| Test | Material | Go criterion |
|---|---|---|
| Parental Sanger for c.2210T>G and c.3006T>G | Parental blood | Variants inherited one from each parent (in trans). If they are in cis, stop and revisit the diagnosis |
| Karyotype (≥100 metaphases), PCS scoring, interphase FISH (e.g. chr 7/8/18/X) | Proband lymphocytes; buccal cells | Mosaic variegated aneuploidy above lab control levels; PCS present |
| Parental PCS scoring | Parental lymphocytes | Raised PCS in the carrier parent would support carrier-level effects and explain the pregnancy loss |
| RT-PCR / RNA-seq of BUB1B | Proband lymphocytes or LCL | Allele-1 transcript reduced (NMD); allele-2 transcript normal, which excludes a splice effect of c.3006T>G |

## Tier 1: Mechanism of p.Asn1002Lys (M1 vs M3 vs M4), months

| Test | Distinguishes | Expected under M1 / M3 / M4 |
|---|---|---|
| BUBR1 quantitative western vs parental and control LCLs | Abundance | low / near-normal / very low |
| Cycloheximide chase (0–8 h); add geldanamycin and MG132 arms, as in Suijkerbuijk 2010 | Turnover; HSP90 and proteasome dependence | about 2× turnover, GA-sensitive, MG132-rescued / normal / — |
| Kinetochore BUBR1 and B56 (PPP2R5) immunofluorescence; pKARD (S676) | KARD / PP2A-B56 function | proportional to abundance / reduced per molecule / absent |
| Nocodazole challenge (mitotic index), live imaging of NEBD→anaphase | SAC strength | shortened / shortened / severely shortened |
| Isogenic knock-in: RPE1 or iPSC BUB1B^N1002K/–, with WT/– and L1012P/– controls | Removes background effects | phenotype tracks L1012P under M1 |
| Rescue: lentiviral WT vs N1002K BUBR1 at matched expression | Whether the mutant works once abundance is restored | N1002K rescues at high expression under M1, does not rescue under M3 |

**Gate:** if N1002K at matched expression rescues the SAC and alignment
(M1), Branch A goes first. If it doesn't (M3 or M4), deprioritise Branch A
and move to Branch B/C candidates.

## Tier 2: Drug testing in patient cells, months

The shortlist is from `04_ranking.md`. Doses stay within clinically
achievable concentrations, and each arm includes a vehicle and a positive
control.

| Candidate | Primary endpoint | Secondary endpoint | Safety readout |
|---|---|---|---|
| Niacin (also NAM/NMN as tool compounds) | Steady-state BUBR1 and half-life | PCS rate, mitotic timing, acetyl-K668 | Proliferation of aneuploid vs euploid cells (FISH) |
| Arimoclomol | BUBR1 half-life (HSF1-dependent?) | HSP70 induction | — |
| Sirolimus (low dose) | pS6, SASP panel, SA-β-gal | BUBR1 level, which tests the autophagy counterargument | Aneuploid fraction |
| Dasatinib ± quercetin | Senescent-cell clearance in patient fibroblasts | SASP | Enrichment of aneuploid survivors |
| Hydroxychloroquine | LC3-II flux, BUBR1 under oxidative stress | Aneuploid fraction over passages | Viability of euploid cells |
| Bortezomib / MG132 | **Probe only**: BUBR1 rescue, which confirms M1 | — | — |

**Combination check:** test every candidate together with vincristine at
clinically relevant exposure, to check for synergy in mitotic-arrest toxicity
in checkpoint-weak cells.

**Go criteria for a candidate:**

1. BUBR1 goes up by ≥30% *and* PCS or missegregation goes down (Branch A),
   or the target phenotype goes down without the aneuploid fraction rising
   (Branch B).
2. No increase in aneuploid-cell survival.
3. No synergy with vincristine toxicity.

## Tier 3: Translation (requires clinicians and ethics)

- Only after Tier 2: a specialist review by paediatric oncology, nephrology
  and clinical genetics, with N-of-1 design considerations.
- Branch A biomarkers that could be followed in blood: lymphocyte BUBR1
  level, PCS rate and blood NAD+ (for niacin).
- These are hypotheses for investigation. Nothing here is a treatment
  recommendation.
