# Track 1 methods: PROBAND01, Model 2 (exomiser-acmg-epcr)

Team: sphenelabs · Model 2 (`exomiser-acmg-epcr`), an upgrade of Model 1 (`exomiser-rank`,
[track1_methods.md](track1_methods.md)).
This document gives the methods of the Model 2 report and methods description form.
The predictions and the full report were submitted through the challenge.

## Summary

The top hypothesis is compound heterozygosity in **BUB1B** (MANE Select ENST00000287598.11),
consistent with mosaic variegated aneuploidy syndrome 1 (MVA1, OMIM 257300):

| | HGVS | Manual ACMG/AMP (Tavtigian points) |
|---|---|---|
| Allele 1 | c.2210T>G p.(Leu737\*) | Pathogenic, 10: PVS1 (8), PM2_Supporting (1), PP4 (1) |
| Allele 2 | c.3006T>G p.(Asn1002Lys) | VUS, 3: PM2_Supporting (1), PM3_Supporting (1), PP4 (1) |

Exomiser ranks BUB1B first of all genes (combined score 0.982). The phenotype-only LIRICAL ranking
places it 4th of 4,983 genes, within the threshold of 20 set before the run. EPCR for the pair is
0.648, and the pair is the top row in all six sensitivity arms.

## Changes from Model 1

| Component | Model 1 | Model 2 |
|---|---|---|
| EPCR | Fixed ladder by rank (0.85, 0.06, 0.03, …) | Formula from the Exomiser gene phenotype score and the ACMG class, with an ordering check and a sensitivity analysis |
| Mechanism gene panel | Multiplied the Exomiser score (×1.0–1.5) | Not used in ranking or EPCR; used for the agreement check and as a sensitivity control |
| Genotype QC | FILTER=PASS | PASS plus GQ / DP / allele-balance thresholds, reported as a genome-wide funnel and checked for every submitted variant; flags only |
| Phenotype evidence | Exomiser only | Exomiser plus a phenotype-only LIRICAL ranking (HPO 2026-09-01) |
| ACMG classification | Exomiser automatic | Manual review of the two BUB1B variants; other candidates keep the automatic class |
| Candidate rules | First two variants by position; `OR*` prefix excluded | Two best variants by Exomiser variant score, then frequency; rarest other het as the hedge allele; `^OR\d+[A-Z]` olfactory-receptor match; rarity thresholds by inheritance mode; clustered-call artefact rule |

## Data

- Proband-only WGS: Sentieon-called VCF (GRCh38 no_alt + hs38d1, GATK hard filters; 5,012,204
  records) and paired-end FASTQ (4 lanes, 2 × 149 bp).
- Phenotype: 8 HPO terms: HP:0002859 rhabdomyosarcoma, HP:0000121 nephrocalcinosis, HP:0004322 short
  stature, HP:0001508 failure to thrive, HP:0003202 skeletal muscle atrophy, HP:0001622 premature
  birth, HP:0001518 small for gestational age, HP:0200067 recurrent spontaneous abortion (in the
  parents).
- No parental samples, so segregation and phase cannot be observed directly.

## Methods

All steps are scripted: `scripts/pipeline/` (run in numeric order), `scripts/run_lirical.sh` and the
modules in `src/rare_disease_mva/` (`qc.py`, `panel.py`, `rank.py`, `epcr.py`, `convergence.py`).
Thresholds and weights are kept in one file, `src/rare_disease_mva/config.py`. Everything ran
locally, and no variant coordinates were sent to external services.

### 1. Normalisation and genotype QC

bcftools 1.21: chr-style contigs, primary assembly, split multiallelics, left-align, remove `*`
alleles. Sex check XY; no autosomal ROH > 3 Mb. Genotype QC thresholds: FILTER=PASS, GQ ≥ 20,
DP ≥ 10, heterozygous allele balance 0.2–0.8, homozygous-alternate allele balance ≥ 0.8. The
thresholds count and flag variants; they do not remove variants from the Exomiser input, because hard
allele-balance filtering could discard a true mosaic variant.

| Step | Flagged | Remaining |
|---|---:|---:|
| Input records | — | 5,012,204 |
| Primary contigs | 116,485 | 4,895,719 |
| FILTER = PASS | 233,846 | 4,661,873 |
| GQ ≥ 20 | 28,480 | 4,633,393 |
| DP ≥ 10 | 33,479 | 4,599,914 |
| Allele balance | 44,012 | 4,555,902 |

### 2. Annotation and genotype-plus-phenotype ranking

ClinVar 2026-10-04, Ensembl VEP 113 and Exomiser 14.0.0 (data release 2406): hiPHIVE + OMIM
prioritisers with all 8 HPO terms; gnomAD v4 and UK10K frequencies; REVEL, MVP, AlphaMissense and
SpliceAI; all modes of inheritance. Exomiser analyses PASS records and removes variants above 2%
population frequency before ranking.

### 3. Phenotype-only ranking (LIRICAL)

LIRICAL v2.4.1 in phenotype-only mode with the HPO ontology and annotations of release 2026-09-01.
Diseases were mapped to genes with the same release, keeping each gene's best disease rank. We set
the reading before the run: BUB1B within the top 20 genes counts as agreement with Exomiser. A
sensitivity run omits HP:0200067, which describes the parents. Both rankings draw on HPO annotations,
so agreement shows robustness across tools rather than independent evidence.

### 4. Mechanism gene panel

116 genes merged from four sources, each recorded per gene: the in-repository differential list; the
PanelApp green genes for R457 (panel 1637 v1.3) and Childhood solid tumours (panel 243 v5.12); every
gene of a disease named "mosaic variegated aneuploidy" in the HPO 2026-09-01 annotations; and SLF2
and SMC5 from the literature. The panel does not affect ranking or EPCR.

### 5. Candidate construction

- One hypothesis per gene from the Exomiser contributing variants, ranked by the unweighted Exomiser
  combined score. The top 20 genes are scored for EPCR, together with the hedge row.
- Recessive genes need two heterozygous hits or one homozygous/hemizygous variant. They take the two
  variants with the highest Exomiser variant score, with ties broken by lower population frequency.
  Dominant genes take the single best variant.
- Heterozygous pairs < 100 bp apart are treated as one allele in cis (EPCR × 0.3).
- Artefact-prone gene families (HLA, KIR, PRSS1, MUC, olfactory receptors matched as `^OR\d+[A-Z]`)
  are excluded.
- Hedge row: for the top recessive gene, the stronger allele is paired with the rarest other
  heterozygous variant (MAX_FREQ ≤ 0.1%) in the same gene.
- **Rarity:** before candidates are built, variants at or above a frequency threshold for the mode
  of inheritance are removed (gnomAD popmax < 0.1% for AR/XR, < 0.01% for AD/XD). The
  secondary-finding gene is exempt. The thresholds were written in [docs/02-plan.md](../docs/02-plan.md)
  before the first submission.
- **Clustered calls:** a candidate is flagged as a likely mapping artefact if a 10 bp window
  containing it holds three or more calls (GATK VariantFiltration convention, cluster size 3,
  window 10) and its own allele balance is below 0.2 (the gnomAD high-quality heterozygous genotype
  standard). Several low-fraction calls carried by the same reads within a few base pairs is the
  pattern of reads from a paralogous or unrepresented sequence placed on the wrong locus. This is
  the only genotype-level exclusion. An isolated low-allele-balance call is kept, because it may be
  a true mosaic variant.

### 6. Read-level and structural review

FASTQ realigned with minimap2 2.28. Read-level review of candidates (MAPQ ≥ 20, BQ ≥ 20), Manta SV
calling, and mosdepth depth in 100 kb bins (MAPQ ≥ 20) for whole-chromosome copy number. Read-backed
phasing was attempted, and the 1000 Genomes high-coverage phased panel was checked.

### 7. Manual ACMG/AMP review

The two BUB1B variants were reviewed code by code, starting from Exomiser's automatic codes, using
ACMG/AMP 2015 with ClinGen SVI refinements and the Tavtigian point scale (P ≥ 10, LP 6–9,
VUS 0–5). Decisions and reasons are in `config/acmg_review.tsv`.

- **PVS1:** kept at full strength after checking NMD on the MANE Select transcript (stop in exon 17
  of 23, 747 nt before the last junction).
- **PP5_Strong:** removed, following the ClinGen SVI recommendation against using database
  assertions as evidence.
- **PM2_Supporting:** kept for both alleles. The maximum credible allele frequency for a recessive
  disorder with prevalence < 1 in 10⁶ is about 5 × 10⁻⁴.
- **PP3:** REVEL was chosen in advance; 0.47 lies between the BP4 (≤ 0.290) and PP3 (≥ 0.644)
  thresholds, so neither applies.
- **PM3_Supporting:** added for p.(Asn1002Lys), which is observed with a pathogenic allele, phase
  unknown. This scores 0.5 on the ClinGen in-trans scale, which corresponds to PM3_Supporting
  (1 Tavtigian point).
- **PP4:** kept at supporting strength.

Other candidates keep Exomiser's automatic class and are labelled "automatic".

### 8. EPCR

```
EPCR = 0.6 × phenotype_score + 0.4 × ACMG_tier
```

- **`phenotype_score`** is the Exomiser gene phenotype score (hiPHIVE + OMIM, 0–1, no panel
  weighting). It is used instead of the combined score because the combined score already includes
  the Exomiser variant score, and the ACMG tier scores the same variant again. With the combined
  score, variant severity is counted twice: a truncating variant in a gene with no phenotype match
  (phenotype score ≈ 0.50, the median over all genes in this run) reaches EPCR ≈ 0.61. This was
  identified while reviewing the ranking before submission. The previous definition is kept as a
  sensitivity arm.
- **PP4 overlap:** ACMG criterion PP4 overlaps the phenotype term conceptually, but removing it
  changes neither BUB1B tier used in the formula.
- **`ACMG_tier`** is P 1.0, LP 0.75, VUS 0.4, LB 0.1, B 0. For a compound heterozygote the weaker
  allele is used.
- **Exclusions:** candidates with a B/LB allele, a single heterozygous variant in a recessive gene,
  an artefact-prone locus, a clustered low-allele-balance call or a failed read review.
- **Fixed values:** hedge rows 0.08 and secondary findings 0.10. A manual adjustment (≤ 0.1) needs
  a written reason; none was used.
- **Ordering check:** F-max sweeps EPCR as a confidence threshold, so only the order matters. The top
  row must be strictly highest, and `epcr.py` enforces this.
- **Sensitivity:** weights 0.5/0.5 and 0.75/0.25; a control that adds a mechanism-panel bonus
  (0.5 / 0.3 / 0.2); the combined-score phenotype term; and a baseline-corrected phenotype term,
  max(0, (score − 0.5) / 0.5).

| Weights (phenotype / ACMG / panel) | Phenotype term | BUB1B pair | Runner-up | Gap |
|---|---|---:|---:|---:|
| 0.6 / 0.4 / 0 (submitted) | phenotype score | 0.648 | 0.485 | 0.163 |
| 0.5 / 0.5 / 0 | phenotype score | 0.607 | 0.471 | 0.136 |
| 0.75 / 0.25 / 0 | phenotype score | 0.710 | 0.507 | 0.203 |
| 0.5 / 0.3 / 0.2 (panel control) | phenotype score | 0.727 | 0.391 | 0.336 |
| 0.6 / 0.4 / 0 | combined score (previous) | 0.749 | 0.653 | 0.096 |
| 0.6 / 0.4 / 0 | baseline-corrected | 0.536 | 0.211 | 0.325 |

## Human inputs and manual review

The CSV is produced by `rank.py` and `epcr.py` from the Exomiser results. Human inputs are part of
the run configuration and are in the repository:

- the manual ACMG review of the two BUB1B variants (`config/acmg_review.tsv`)
- the artefact-gene list
- the rarity and clustered-call thresholds (fixed values in `config.py`, applied by code)
- LZTR1 as a secondary-finding gene

No row was added, removed, reordered or re-scored by hand outside these inputs.

## Compound heterozygotes, hedge and secondary findings

- **Compound heterozygotes:** recessive genes output a compound-heterozygous pair or a single
  homozygous/hemizygous variant, and the EPCR uses the weaker allele's ACMG class.
- **Hedge row:** pairs c.2210T>G with the non-coding c.2679-1026A>G, as a fallback for a non-coding
  second allele. Its second allele has no ACMG class, so it keeps the fixed value 0.08 and sorts below
  the formula-scored rows. Its position reflects that rule, not lower biological plausibility.
- **Secondary finding:** LZTR1 c.2244C>G p.(Tyr748\*), heterozygous, is reported as secondary.
  LZTR1 loss of function is associated with LZTR1-related schwannomatosis and, when biallelic, with
  AR Noonan syndrome; it does not explain the phenotype.

## Discussion and limitations

- **Other MVA genes.** MVA is genetically heterogeneous, and an affected child is expected to carry
  biallelic variants in one MVA gene. CEP57, TRIP13 and MAD1L1 carry no rare coding or splice
  variant, and no exon-level deletion or SV was found in BUB1B, CEP57 or TRIP13.
- **Artefact or mosaic.** MTCH2 c.229A>T lies in a cluster of four calls carried by the same 11
  reads and is excluded. CACNB1 c.1797A>G is an isolated call (allele balance 0.19) and is kept as a
  low-confidence row.
- **Lower rows are background.** A typical genome carries about 100 genuine loss-of-function
  variants, so a rare truncating VUS in a gene without a phenotype match is expected by chance.
  Rows 2–8 are near-tied. FLG p.(Arg2447\*) is a common allele (gnomAD popmax 0.77%), a known
  ichthyosis vulgaris / atopic dermatitis allele, and is excluded by the dominant rarity threshold.
- **Singleton data.** Trans phase and inheritance are inferred, not observed.
- **Novel second allele.** p.(Asn1002Lys) has no functional data and remains a VUS.
- **Scoring.** The EPCR weights are a stated heuristic, not fitted, and manual ACMG review covers
  only the two BUB1B variants.
- **Shared phenotype knowledge.** Exomiser and LIRICAL both use HPO annotations, so their agreement
  is not independent evidence. Phen2Gene was trialled and replaced because its 2021 knowledge base
  lacks newer MVA genes.
- **Mosaic aneuploidy.** Bulk short-read WGS cannot detect low-level mosaic aneuploidy; karyotyping
  or interphase FISH would be needed.

## Public data and run time

- **Public data:** GRCh38 no_alt + hs38d1; GENCODE v46 / MANE Select; ClinVar 2026-10-04; gnomAD v4
  and UK10K; REVEL, MVP, AlphaMissense and SpliceAI via Exomiser 2406; OMIM and Orphanet disease
  annotations, MGI mouse and ZFIN zebrafish phenotypes and STRING interactions used by hiPHIVE; HPO
  release 2026-09-01; PanelApp panels 1637 v1.3 and 243 v5.12; 1000 Genomes high-coverage phased
  panel. No proprietary data.
- **Run time:** on one 20-core / 31 GB workstation, VCF QC, ClinVar, Exomiser and ranking take about
  1 h, realignment 5–8 h, and Manta and mosdepth 1–2 h. Local compute only, no cloud cost.

## Use of AI assistance

Anthropic Claude Code (VS Code extension), model Claude Opus 5.5; Claude Pro individual plan,
Anthropic Consumer Terms. "Help improve Claude" (use of chats for model training) was set to OFF,
and no feedback (ratings) was provided on outputs. It was used interactively to write and run the
analysis code, review intermediate results and draft the report. Raw FASTQ/BAM/VCF files were
processed locally and never uploaded as files. All variant interpretation was checked by the team
against the primary sources cited.

## Key references

- Smedley D et al. Next-generation diagnostics and disease-gene discovery with the Exomiser. *Nat Protoc* 2015;10:2004–15. doi:10.1038/nprot.2015.124
- Bone WP et al. Computational evaluation of exome sequence data using human and model organism phenotypes improves diagnostic efficiency. *Genet Med* 2016;18:608–17. doi:10.1038/gim.2015.137
- Robinson PN et al. Interpretable clinical genomics with a likelihood ratio paradigm. *Am J Hum Genet* 2020;107:403–417. doi:10.1016/j.ajhg.2020.06.021
- Richards S et al. Standards and guidelines for the interpretation of sequence variants. *Genet Med* 2015;17:405–24. doi:10.1038/gim.2015.30
- Tavtigian SV et al. Fitting a naturally scaled point system to the ACMG/AMP variant classification guidelines. *Hum Mutat* 2020;41:1734–1737. doi:10.1002/humu.24088
- Hanks S et al. Constitutional aneuploidy and cancer predisposition caused by biallelic mutations in BUB1B. *Nat Genet* 2004;36:1159–61. doi:10.1038/ng1449
- Karczewski KJ et al. The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature* 2020;581:434–443. doi:10.1038/s41586-020-2308-7
- Li H. Toward better understanding of artifacts in variant calling from high-coverage samples. *Bioinformatics* 2014;30:2843–51. doi:10.1093/bioinformatics/btu356
- Acuna-Hidalgo R et al. Post-zygotic point mutations are an underrecognized source of de novo genomic variation. *Am J Hum Genet* 2015;97:67–74. doi:10.1016/j.ajhg.2015.05.008
- MacArthur DG et al. A systematic survey of loss-of-function variants in human protein-coding genes. *Science* 2012;335:823–8. doi:10.1126/science.1215040
- Broad Institute. GATK VariantFiltration documentation: https://gatk.broadinstitute.org/hc/en-us/articles/360036460072-VariantFiltration

The full reference list (39 sources) is in the submitted report.
