# Stage 1: Variant evidence dossier

Plan: `docs/04-track2-plan.md`, Stage 1. Generated 2026-10-18.

Inputs:

- `reports/track1_methods.md` for the variant calls and read-level QC.
- `scripts/track2/s1_structure.py` produces `tables/s1_structure*.tsv`.
- `scripts/track2/s1_conservation.py` produces `tables/s1_conservation*.tsv`.
- `scripts/track2/s0_pubmed.py` produces `tables/s0_pubmed_hits.tsv`.

## 1. Genotype recap (from Track 1)

| | Variant | Consequence | Read support |
|---|---|---|---|
| Allele 1 | chr15:40209701 T>G, c.2210T>G p.(Leu737\*) | Nonsense. ClinVar P/LP for MVA1 (VCV 533901) | 25/46 alt reads, both strands |
| Allele 2 | chr15:40220612 T>G, c.3006T>G p.(Asn1002Lys) | Missense, last exon (no NMD). Novel; gnomAD 1/1.46M | 13/28 alt reads, both strands |

- **Phase:** the variants are 10.9 kb apart, so phase is not observable from
  short reads, read-backed phasing or the 1000 Genomes panel. Trans is
  **inferred**.
- **Splice effect of allele 2:** c.3006T>G is in the last exon. There was no
  genome-wide SpliceAI score at Track 1. Recorded as **not assessed**; it would
  be resolved by RNA from patient cells.
- **Allele 1 truncation:** p.Leu737\* removes the whole C-terminal domain
  (766–1050). It is expected to trigger NMD, since the stop codon is far
  upstream of the last exon junction. That fits Suijkerbuijk et al. 2010
  (PMID 20516114), where transcripts of MVA truncating alleles were absent.

## 2. Structural context of Asn1002 (AlphaFold v6, O60566)

| Residue | Role | pLDDT | Relative SASA | Burial |
|---|---|---|---|---|
| **N1002** | **patient** | 91.1 | 0.16 | buried |
| L1012 | MVA1 L1012P | 95.3 | 0.00 | buried |
| T1008 | PLK1 phosphosite | 92.1 | 0.06 | buried |
| R814 | MVA1 R814H | 90.2 | 0.16 | buried |
| L844 | MVA1 L844F | 85.9 | 0.01 | buried |
| I909 | MVA1 I909T | 87.0 | 0.03 | buried |
| Q921 | MVA1 Q921H | 83.2 | 0.01 | buried |

- N1002 sits in a confidently predicted part of the pseudokinase C-lobe
  (pLDDT 91). It is mostly buried, making side-chain contacts with residues
  973, 977, 978, 998 and 999.
- N1002 Cα is 12.3 Å from L1012 and 10.7 Å from the PLK1 site T1008.
- Every known MVA1 missense variant in the C-terminal domain is at a buried
  position. These variants lie throughout the domain, not at an interaction
  surface, which fits a **fold-stability** class of variants.
- Asn→Lys puts a long, positively charged side chain into a buried position
  where the native Asn forms side-chain contacts. This is very likely to
  disrupt local packing or hydrogen bonding. *Inference from the model, not
  measured.*

## 3. Conservation

From Ensembl Compara one-to-one orthologs (about 120 vertebrates), aligned
with MAFFT:

| Position | Identity to human | Residues seen | Lys ever seen? | BUB1 paralog |
|---|---|---|---|---|
| **N1002** | **0.97** (117/120) | N, D, S, E (1 each) | **no** | **N** |
| L1012 (MVA1) | 0.97 | L (V, R, T, E once each) | no | L |
| T1008 (PLK1) | 0.94 | T, S, R, P | no | L |
| R814 (MVA1) | 0.97 | R | no | R |
| L737 (nonsense site) | 0.92 | — | — | V |

- The local motif `W-N-K-F-F-V-R-I-L-N` (residues 993–1002) is conserved from
  fish to mammals. The human BUB1 paralog has `W-N-E-F-F-H-V-M-L-N`. Asn1002
  is the final Asn of the shared `W-N-x-F-F-x-x-x-L-N` core, so it is
  conserved across the **BUB1/BUBR1 family**.
  That is evidence of structural constraint that predates the split between
  the two paralogs.
- The only substitutions ever observed are small polar or acidic residues
  (D, S, E). A basic Lys has never been observed.

## 4. Literature on how C-terminal-domain missense variants act

| Finding | Source | Tier | Bearing |
|---|---|---|---|
| In biallelic MVA, a missense allele pairs with a truncating allele. MVA mutations fall into two classes: (i) specific checkpoint or attachment defects, (ii) **low BUBR1 abundance**. The low-abundance class is "absent transcripts from truncating mutants combined with **high protein turnover of missense mutants**", and "the amino acid change consistently occurs **in or near the BUBR1 kinase domain**". Adding BUBR1 back restored checkpoint activity. | Suijkerbuijk 2010 *Cancer Res* (20516114) | E2 | Supports **M1** for kinase-domain missense; rescue by more BUBR1 is proven in patient cells |
| Full text (PMC2887387): kinase-domain substitutions **I909T, R727C and L1012P** lower BUBR1 levels, and L844F was poorly expressed. With cycloheximide, "the protein turnover of the substitution mutant proteins I909T and R727C was increased ~2-fold". **HSP90 inhibition** (geldanamycin) plus cycloheximide "removed virtually all mutant BUBR1 protein", and the proteasome inhibitor MG132 prevented this. By contrast, R550Q (central region) and Q921H were normal in level and function, and Y155C (TPR) is a separation-of-function allele. | Suijkerbuijk 2010 (full text) | E2 | **Strong support for M1**: the closest-matching variants are HSP90-dependent clients that are cleared by the proteasome. **Safety:** HSP90 inhibitors are likely harmful |
| Residues that bind ATP in conventional kinases are essential for BUBR1 **conformational stability**. Catalysis is dispensable. | Suijkerbuijk 2012 *Dev Cell* (22698286) | E2 | The domain acts as a stability scaffold, so a buried substitution plausibly destabilises it |
| The pseudokinase domain is required for KARD phosphorylation, then PP2A-B56 recruitment, checkpoint silencing and chromosome alignment. | Gama Braga 2020 *Cell Rep* (33207204) | E2 | Supports **M3**: folded protein with a defective domain could still disrupt alignment |
| The PLK1 sites in the kinase domain are needed for chromosome alignment. | Matsumura 2007 *JBC* (17376779) | E2 | T1008 is 10.7 Å away; M3 is possible |
| Function is steeply dosage-dependent: PCS appears with any loss; aneuploidy appears below about 50% BUBR1. | Bohers 2008 *Hum Genet* (18932004); Matsuura 2006 (16411201) | E2 | Even partial restoration could matter |
| Mouse BubR1 L1002P (= human **L1012P**): H/L1002P mice have a total-BUBR1 reduction similar to H/H, but some progeroid features are **milder**. X753/L1002P is lethal. | Sieben 2020 *JCI* (31738183) | E3 | Allele-specific effects beyond abundance (argues for keeping M2/M3 open). Caution: "L1002P" in this paper is mouse numbering, not the patient's N1002 |
| Truncating MVA alleles reduce lifespan even in heterozygous mice. | Wijshake 2012 *PLoS Genet* (23300461) | E3 | Background on allele 1 |

## 5. Competing-mechanism matrix and prior weights

The weights are subjective priors, justified by the evidence above. Stage 4
uses them as scenario weights and tests how sensitive the ranking is to them.

| ID | Model | For | Against | Prior |
|---|---|---|---|---|
| M1 | Unstable (HSP90-dependent) protein, low abundance | Buried, family-conserved residue; Lys never seen; 3 of 3 tested buried kinase-domain MVA missense variants (I909T, R727C, L1012P) have high turnover, and L844F is poorly expressed (Suijkerbuijk 2010); the pseudokinase fold depends on conformational stability | No measured abundance for this allele; Q921H (also buried in the model) was normal | **0.50** |
| M2 | Normal abundance, impaired checkpoint | Suijkerbuijk 2010 class (i) exists | Class (i) examples are N-terminal or central (Y155C, R550Q), not C-terminal | 0.05 |
| M3 | Folded but defective KARD / PP2A-B56 support or alignment | Pseudokinase-domain role (Gama Braga 2020); near T1008 PLK1 site; L1012P allele-specific effects (Sieben 2020) | Residue is buried, not surface-exposed | 0.25 |
| M4 | Near-null | AlphaMissense 0.92 | A null/null genotype (X753/L1002P-like) is lethal in mice; the child is alive | 0.10 |
| M5 | Not causal or in cis | Phase not observed; VUS | Exomiser rank 1; textbook MVA1 phenotype; ROH-negative | 0.10 |

**Gate decision:** proceed with all of M1–M4, weighted by these priors. M5 is
not excluded, and the report will say so up front. The priority experiment
that would shift the weights most is a **BUBR1 western blot plus
cycloheximide chase** in patient lymphocytes or fibroblasts (M1 vs M2/M3),
together with **parental Sanger testing** (M5).

## 6. ACMG notes (for the report, not reclassified here)

- PM2_Supporting (gnomAD) and PP4 (phenotype) were already applied in Track 1.
- PM1 (critical domain without benign variation) is now **arguably
  applicable**: the buried C-terminal domain holds most MVA1 missense
  variants. This needs a gnomAD missense-constraint check for the region
  before it is claimed.
- PM3 (in trans with a pathogenic variant) **cannot be applied** until phase
  is confirmed.
- The variant stays a VUS for now. PS3 needs functional data.
