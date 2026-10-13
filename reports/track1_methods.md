# Track 1 methods: PROBAND01

> Draft. Rename to `<hf-username>_track1_report.md` before submitting.

## Summary

The proband carries **compound-heterozygous candidate variants in _BUB1B_**
(MANE Select transcript ENST00000287598.11), consistent with
**mosaic variegated aneuploidy syndrome 1 (MVA1, OMIM 257300)**:

| | Variant (GRCh38) | HGVS | Evidence |
|---|---|---|---|
| Allele 1 | chr15:40209701 T>G | c.2210T>G p.(Leu737\*) | Nonsense, upstream of the C-terminal kinase domain. ClinVar Pathogenic/Likely pathogenic for MVA1 (VCV 533901, multiple submitters, no conflicts). gnomAD v4 AC 115 / 1,461,846 exomes, 0 homozygotes. Exomiser ACMG: Pathogenic (PVS1, PP5_Strong, PM2_Supporting, PP4) |
| Allele 2 | chr15:40220612 T>G | c.3006T>G p.(Asn1002Lys) | Missense in the kinase domain (last exon, so no NMD). Not in ClinVar. gnomAD v4 AC 1 / 1,461,878 exomes, absent from genomes. AlphaMissense 0.92, MVP 0.85, REVEL 0.47 (borderline). Exomiser ACMG: VUS (PM2_Supporting, PP4) |

The genotype fits the established MVA1 architecture: one truncating allele plus
one hypomorphic (here missense, kinase-domain) allele. Biallelic complete loss of
BUB1B is not expected to be viable. The phenotype also fits MVA1:
rhabdomyosarcoma, prenatal growth restriction, short stature, failure to thrive,
and parental reproductive loss, which would be consistent with carrier-level
aneuploidy/premature chromatid separation in a parent.

## Data

- Proband-only WGS: Sentieon-called VCF (GRCh38 no_alt + hs38d1, GATK hard
  filters) and paired-end FASTQ (4 lanes, 2×149 bp).
- Phenotype: 8 HPO terms: HP:0002859, HP:0000121, HP:0004322, HP:0001508,
  HP:0003202, HP:0001622, HP:0001518, HP:0200067.
- No parental samples, so segregation and phase are not directly observable.

## Pipeline

All steps are in `scripts/pipeline/` (run in numeric order), with configs in
`config/` and Python in `src/rare_disease_mva/`. Tools were installed from
bioconda (micromamba). Versions: bcftools/samtools/htslib 1.21, minimap2 2.28,
Exomiser 14.0.0 (data 2406), Ensembl VEP 113, ClinVar 2026-10-04, GENCODE v46,
Manta, mosdepth 0.3.10.

1. **QC and normalisation** (`20_prep_vcf.sh`): renamed contigs to `chr`
   style, restricted to the primary assembly, split multiallelics, left-aligned
   against the reference, and removed `*` alleles. 4.73M PASS records,
   Ti/Tv 2.02.
   - **Sex:** XY (chrX non-PAR het/hom-alt ratio 0.013; 9,778 chrY calls).
   - **ROH** (`bcftools roh`): no autosomal ROH > 3 Mb, so no evidence of
     consanguinity or large isodisomy. A recessive cause is therefore most
     likely compound heterozygous.
2. **ClinVar annotation**, plus a genome-wide scan of ClinVar P/LP genotypes.
3. **Exomiser** (`30_exomiser.sh`, `config/exomiser/genome-analysis.yml`):
   hiPHIVE + OMIM prioritisation with all 8 HPO terms. gnomAD v4 / UK10K
   frequencies; REVEL, MVP, AlphaMissense and SpliceAI pathogenicity; all modes
   of inheritance, with compound-het modes weighted up.
4. **Candidate ranking** (`src/rare_disease_mva/rank.py`): one hypothesis per
   gene from Exomiser's contributing variants.
   - **Recessive genes** must have two hits.
   - **Cis penalty:** pairs of hets less than 100 bp apart are penalised, since
     they are likely a single complex allele in cis (this removed apparent
     FANCD2 and SERPINA1 "compound hets").
   - **HLA/KIR/paralog-artefact genes** are excluded.
   - **Panel boost:** genes on a phenotype-driven differential panel
     (`config/gene_panels.tsv`) get a boost. The panel covers MVA genes, other
     chromosome-instability disorders, sarcoma predisposition, and
     nephrocalcinosis.
5. **Realignment** (`10_align.sh`): minimap2 `-x sr`, then fixmate, sort,
   merge and markdup (1.08 billion reads, 13% duplicates). The BAM feeds
   read-level review, Manta SV calling and mosdepth depth (`40_sv_cnv.sh`):
   1 kb bins, per-base depth over the MVA genes, and MAPQ ≥ 20 100 kb bins for
   the aneuploidy check.

## Results

- **Exomiser:** BUB1B (AR) ranked 1 of all genes (combined score 0.982,
  p < 1e-4), with exactly the two variants above as its contributing pair.
  Next were HLA-DRB1 (an artefact-prone region) and FANCD2 (AR; two deletions
  2 bp apart at the c.1278+1 donor, which is likely a single cis event, and
  FANCD2 has a pseudogene).
- **ClinVar P/LP genotypes genome-wide:** BUB1B c.2210T>G (above);
  heterozygous carrier states in FLG, GNRHR and RBM8A; LZTR1 c.2244C>G
  p.(Tyr748\*) (see secondary findings); and PRSS1 at a VAF of 0.15 in the
  paralog-rich TRB/PRSS locus (interpreted as a mapping artefact).
- **MVA genes CEP57 and TRIP13:** no rare coding or splice variants.
- **Read-level review** (realigned BAM, MAPQ ≥ 20, BQ ≥ 20): both BUB1B
  variants are clean heterozygous calls. c.2210T>G has 25/46 alt reads in the
  VCF and alt support on both strands in the pileup. c.3006T>G has 13/28 alt
  reads, also on both strands.
- **Phase:** the two BUB1B variants are 10.9 kb apart, beyond read-pair
  distance. The 20 kb window holds only four heterozygous sites, 0.8 to 6.8 kb
  apart, so read-backed phasing can't bridge them either. Neither is in the 1000 Genomes high-coverage phased panel
  (3,202 genomes), so population-based phasing isn't possible. Trans
  configuration is inferred, not observed (see Limitations).
- **Mosaic aneuploidy screen from the VCF** (per-chromosome het-SNP B-allele
  fraction): chr20–22 show modestly wider allele balance and higher depth.
  This pattern is typical of segmental-duplication and GC artefacts on small
  chromosomes rather than mosaic trisomy: a depth increase of that size would
  imply near-complete trisomy, which is inconsistent with the median BAF of
  0.48. The BAM check below resolves this.
- **Mosaic aneuploidy check from the BAM** (mosdepth, MAPQ ≥ 20, median of
  covered 100 kb bins; autosomal median 45.6×): every autosome is within
  0.986–1.044× of the genome median. The highest values are on the GC-rich
  chromosomes 19, 17, 22 and 16 (1.03–1.04×), the pattern expected from GC
  bias. chrX and chrY are both 0.49×, consistent with a 46,XY karyotype. There
  is no whole-chromosome gain or loss above roughly 10% of cells. **This does
  not exclude MVA.** In MVA, different cells gain or lose different
  chromosomes, and aneuploid cells are often a small, variable fraction of
  blood. Bulk WGS averages over that, so karyotyping or interphase FISH would
  be needed to see it.
- **SV/CNV** (Manta 1.6.0, 9,577 PASS calls; mosdepth): no SV calls and no
  exon-level copy change in BUB1B, CEP57 or TRIP13.
  - **Per-exon depth** (MANE exons, relative to 45.6×): BUB1B 0.75–1.42×,
    TRIP13 0.67–1.07×.
  - **CEP57 exons 4 and 8** read at 0.59× and 0.57×. Both are short (122 bp
    and 78 bp) and AT-rich (34–38% GC). Depth slopes down gradually over
    about 300 bp instead of stepping down at breakpoints, and Manta found no
    candidate split or discordant reads there. We interpret these as
    coverage dips, not deletions.
  - **Large Manta calls are artefacts:** a 1.9 Mb duplication over FANCD2,
    a 97 Mb duplication on chr3 and a 48.8 Mb deletion on chr6 all have
    normal depth inside the call (0.98–1.02×), the same as their flanks.

- **Non-coding review** (VEP 113 over BUB1B, CEP57 and TRIP13 ±10 kb, all
  FILTER values): CEP57 and TRIP13 have no variants below 1% gnomAD AF other
  than common intronic SNPs. The BUB1B promoter and UTRs carry only common
  variants. The only rare non-coding BUB1B variant is c.2679-1026A>G (absent
  from gnomAD). It lies in a CA/TA-repeat tract and creates no new AG or GT
  splice dinucleotide, so it is a weak candidate for a hypomorphic second
  allele. No genome-wide SpliceAI score was available for it.

## Secondary findings

- **LZTR1** c.2244C>G p.(Tyr748\*), heterozygous: LZTR1 loss-of-function is
  associated with LZTR1-related schwannomatosis (AD, reduced penetrance) and,
  when biallelic, with AR Noonan syndrome. It doesn't explain the phenotype.
  Reported as a secondary finding for clinical consideration.

## Submission

The ranked CSV comes from `rank.py`. The BUB1B pair carries EPCR 0.85. The
next row pairs BUB1B p.(Leu737\*) with the next-rarest BUB1B het (intronic
c.2679-1026A>G), as a hedge against a non-coding hypomorphic second allele.
The remaining rows are lower-confidence Exomiser candidates.

## Limitations

- Singleton: no parental genotypes, so trans phase and inheritance are
  inferred. Parental Sanger testing (or long-read phasing) would confirm it.
- p.(Asn1002Lys) is novel and has no functional data. Its classification rests
  on rarity, kinase-domain location, in-silico prediction and phenotype fit.
  An assay of BUB1B protein level, or of mitotic checkpoint function and
  premature chromatid separation in patient cells, would be informative.
- Short-read WGS has limited sensitivity in segmental duplications and for
  low-level mosaicism.

## Use of AI assistance

Analysis was run with Claude Code (Anthropic, Claude Opus 5.5).
**[Fill in: plan/tier and training opt-out setting, e.g. "Claude Pro, training
opt-out enabled, no feedback provided".]**
