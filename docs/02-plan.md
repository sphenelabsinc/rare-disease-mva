# Track 1 plan: variant prediction

Goal: a ranked list of candidate causal variant(s) for proband `WGS_EX2312012`,
plus a GitHub link and a methods writeup. Deadline **24 Oct 2026**, with
**6 submissions** allowed. See `01-variant-analysis.md` for the scoring rules.

## 1. What we have

| Input | Notes |
|---|---|
| `WGS_EX2312012_HGWCNDSX7.vcf.gz` (+ `.tbi`) | One sample (proband only, no parents). GRCh38 `no_alt_analysis_set` + decoys, contigs named `1..22,X,Y` with no `chr` prefix. Called with Sentieon 202308 (Haplotyper → GVCFtyper), then GATK hard filters. 5,012,204 records, of which 4,740,790 PASS (4.08M SNVs, 0.93M indels/other). Unannotated: no VEP/CSQ fields. Has `PGT`/`PID` physical phasing. |
| 8 × FASTQ (L001–L004, R1/R2), ~85 GB | Raw reads, roughly 30–40× WGS. Needed for anything the small-variant VCF misses: SVs, CNVs, mosaicism, read-backed phasing. |
| `Challenge_Clinical_Phenotype_1.docx` | HPO terms, listed below. |

Phenotype (HPO):

- HP:0002859 Rhabdomyosarcoma
- HP:0000121 Nephrocalcinosis
- HP:0004322 Short stature
- HP:0001508 Failure to thrive
- HP:0003202 Skeletal muscle atrophy
- HP:0001622 Premature birth
- HP:0001518 Small for gestational age
- HP:0200067 Recurrent spontaneous abortion (parental history)

The organisers' notes point toward a **chromosomal instability disorder** and
say to treat parental miscarriages as a real signal.

Machine: 20 cores, 31 GB RAM, ~1 TB free, Docker available. Not installed:
bcftools, samtools, aligners or Java, so we'll use Docker images (or the
`pysam` wrappers).

## 2. Working hypothesis and differential

**Leading hypothesis: mosaic variegated aneuploidy (MVA)**, which is also the
hackathon's name. The phenotype fits closely: prenatal growth restriction,
short stature, childhood cancer (rhabdomyosarcoma and Wilms tumour are both
reported in MVA), and reproductive loss consistent with aneuploidy. MVA is
autosomal recessive:

- **BUB1B** (MVA1): often compound heterozygous, with one truncating allele
  plus a hypomorphic allele. The literature reports hypomorphic alleles that
  are missense, splice-region, or **non-coding / low-expression haplotypes**,
  which a coding-only filter would miss. Verify the specific reported alleles
  against ClinVar and the primary papers; don't rely on memory.
- **CEP57** (MVA2)
- **TRIP13** (MVA3): Wilms tumour predisposition

The scoring rubric mentions "one of the two variants in a heterozygous pair",
which also suggests a compound-heterozygous answer is plausible.

**Don't anchor too early.** Score the whole genome, then check these
alternatives explicitly:

- Chromosome breakage / DNA repair: Bloom (BLM), Fanconi anaemia genes
  (FANCA…, which also fit growth restriction and renal anomalies), NBN, LIG4,
  RECQL4
- Paediatric sarcoma predisposition: TP53 (Li-Fraumeni), CMMRD (biallelic
  MLH1/MSH2/MSH6/PMS2), DICER1, NF1, HRAS (Costello), PTCH1/SUFU
- Primordial dwarfism with cancer: TRIM37 (Mulibrey nanism)
- **Nephrocalcinosis as a second or incidental finding**: it may have a
  separate cause from the main diagnosis (or be iatrogenic from prematurity).
  Screen CYP24A1, SLC34A1 and other hypercalciuria / tubulopathy genes.
  Secondary findings don't hurt the automated score.

## 3. Pipeline

Each step gets a script or module in `src/` and a notebook in `notebooks/`
for inspection.

### Phase A: small-variant prioritisation (VCF only; do first, cheap)

1. **QC**: check sex (X heterozygosity, Y calls), Ti/Tv, het/hom ratio, and
   runs of homozygosity (`bcftools roh`). Consanguinity would make recessive
   homozygous candidates more likely.
2. **Normalise**: split multiallelics and left-align (`bcftools norm`).
   Keep PASS records, but keep near-threshold FILTER calls in the candidate
   genes for manual review.
3. **Annotate** with Ensembl VEP (Docker `ensemblorg/ensembl-vep`, GRCh38
   cache): consequence, gnomAD v4 allele frequencies, ClinVar, plus the
   SpliceAI, REVEL, AlphaMissense and CADD plugins.
4. **Exomiser** (Docker) using the 8 HPO terms gives an independent,
   phenotype-driven ranking. It's the standard baseline for this kind of
   task, so we should at least match it.
5. **Our own scoring**:
   - rarity (gnomAD popmax AF; < 0.1% for recessive, < 0.01% for dominant)
   - predicted impact (LoF, splice, deleterious missense)
   - ClinVar evidence
   - phenotype match (HPO similarity, with extra weight for the
     chromosomal-instability gene set above)
   - inheritance fit: homozygous or **two hits in the same gene** for
     recessive genes
6. **Phase compound-het pairs** using `PGT`/`PID` where available and
   read-backed phasing (WhatsHap) once a BAM exists. Two variants on the same
   haplotype (cis) don't count as biallelic.

### Phase B: from the FASTQs (run in the background, started early)

7. **Align** to the same reference (GRCh38 no_alt + hs38d1), using bwa-mem2
   with a prebuilt index (building one needs more RAM than we have) or bwa,
   then sort, mark duplicates and index. Expect several hours on 20 cores,
   and plan for about 100 GB of intermediates.
8. **SV/CNV calling**: Manta for SVs; read-depth CNVs with CNVnator or
   cn.mops-style binning. Look for exon-level deletions in BUB1B, CEP57 and
   TRIP13, which the small-variant VCF can't show.
9. **Mosaic aneuploidy check**: compute read depth and B-allele frequency per
   chromosome. Low-level whole-chromosome gains or losses in blood DNA would
   directly support MVA. This isn't scored, but it's strong evidence for the
   methods writeup.
10. **Non-coding review** of the top genes: promoter, UTR and deep-intronic
    variants (SpliceAI scores at all positions, not just near exons), plus
    any known hypomorphic haplotype.

### Phase C: ranking and submission

11. Merge the Exomiser ranks, our scores and the SV/CNV calls into one ranked
    list with calibrated confidence scores. F-max rewards a short list with
    sharp scores, so put the strong candidates (including both members of a
    compound-het pair) in the top 10, then add lower-confidence and
    secondary findings with low scores.
12. **Use the 6 submissions carefully**: start with a baseline submission
    once the format is confirmed, then refine. Keep at least 2 for the final
    days.
13. **Methods writeup**: QC, annotation sources and versions, scoring logic,
    phasing evidence, CNV/mosaicism evidence, the differential and why each
    alternative was ruled out, and limitations (singleton, no parental
    segregation).

## 4. Timeline (15 days left)

| Dates | Work |
|---|---|
| 9–10 Oct | Confirm the submission format; QC and normalisation; start alignment in the background; VEP and Exomiser |
| 11–13 Oct | Scoring, review the candidate genes, compound-het phasing; **first submission** |
| 14–17 Oct | SV/CNV and mosaicism analysis; non-coding review; revise the ranking |
| 18–21 Oct | Methods writeup; clean up the repo for the GitHub link |
| 22–24 Oct | Final submissions (deadline 24 Oct) |

## 5. Open questions

- **Submission file format**: which columns, which coordinate convention
  (GRCh38, `chr` prefix or not), HGVS or chrom/pos/ref/alt? Check the
  challenge Space before writing the output step.
- Is a confidence score column used for F-max, or only rank order?
- Do SV/CNV predictions count as scoreable variants?
- Will other probands (`Challenge_Clinical_Phenotype_2…`) be released? If so,
  keep the pipeline generic and avoid hardcoding to this case.
- Privacy: the phenotype document is confidential. Keep patient details out
  of the public repo and writeup beyond what the family has published.
