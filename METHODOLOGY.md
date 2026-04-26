# Methodology

## Project Overview
This project performs a comprehensive multi-omics analysis of breast cancer using publicly available datasets. The study integrates RNA-seq, ChIP-seq, GWAS, and ncRNA data to identify key regulatory mechanisms driving breast cancer progression.

---

## Task 1: RNA-seq Data Retrieval, Alignment and Quantification

### Objective
Retrieve raw RNA-seq data, align to reference genome, and generate gene expression count matrix.

### Dataset
- GEO Accession: GSE183947
- 12 samples: 6 tumor + 6 normal breast tissue
- Platform: Illumina HiSeq 2000 (paired-end, 150bp)

### Tools Used
- Galaxy Project (usegalaxy.org) — cloud-based bioinformatics platform
- fasterq-dump — SRA to FASTQ conversion
- RNA STAR v2.7.11b — splice-aware RNA-seq aligner
- featureCounts — read quantification
- Column Join — count matrix merging

### Steps
1. Downloaded 12 SRA files directly from NCBI SRA via AWS links
2. Converted SRA to paired-end FASTQ using fasterq-dump
3. Downloaded hg38 reference genome and ncbiRefSeq GTF annotation
4. Aligned all 12 samples to hg38 using RNA STAR with sjdbOverhang=149
5. Generated BAM files and per-gene read counts
6. Ran featureCounts on all 12 BAM files with GTF annotation
7. Merged 12 count tables into single matrix using Column Join
8. Output: FeatureCounts_Mod.txt (59,251 genes x 12 samples)

---

## Task 2: Differential Gene Expression Analysis

### Objective
Identify genes significantly differentially expressed between tumor and normal samples.

### Tool
- DESeq2 v1.42 in R v4.3.3

### Steps
1. Loaded count matrix and sample metadata into R
2. Constructed DESeqDataSet with design formula ~condition
3. Filtered genes with fewer than 10 total reads
4. Ran DESeq2 normalization and statistical testing
5. Extracted results comparing tumor vs normal
6. Applied FDR correction threshold padj < 0.01
7. Generated volcano plot and heatmap of top 50 DEGs

### Results
- Total genes tested: 31,840
- Significant DEGs: 942 (padj < 0.01)
- Upregulated in tumor: 411
- Downregulated in tumor: 531
- Top upregulated: SOX11, CXCL9, BIRC5, PLK1, CCNB1
- Top downregulated: DES, MYOCD, CIDEA

---

## Task 3: Pathway Enrichment Analysis

### Objective
Identify biological pathways and processes enriched in our DEG list.

### Tool
- clusterProfiler in R
- org.Hs.eg.db for gene ID conversion
- KEGG and GO databases

### Steps
1. Converted 942 DEG gene symbols to Entrez IDs using bitr()
2. Ran GO enrichment (Biological Process) using enrichGO()
3. Ran KEGG pathway enrichment using enrichKEGG()
4. Applied FDR correction (BH method, padj < 0.05)
5. Generated barplot for GO terms and dotplot for KEGG pathways

### Results
- Significant GO terms: 304
- Significant KEGG pathways: 9
- Top GO term: Mitotic cell cycle process
- Top KEGG pathway: Cell cycle (hsa04110)

---

## Task 4: ncRNA Regulatory Analysis

### Objective
Identify miRNA-gene regulatory interactions involving our DEGs.

### Tool
- multiMiR R package (queries miRTarBase, miRecords, TarBase)

### Steps
1. Used all 942 DEGs as query genes
2. Queried multiMiR for validated miRNA-target interactions
3. Filtered to keep only experimentally validated interactions
4. Ranked miRNAs by number of DEG targets
5. Built regulatory network of top 10 miRNAs vs top DEGs

### Results
- Total validated interactions: 90,452
- Unique miRNAs identified: 2,534
- Top miRNAs: hsa-miR-34a-5p, hsa-let-7a-5p, hsa-miR-15a-5p
- Key targets: CCNB1, E2F3, MKI67, SOX4, YAP1

---

## Task 5: ChIP-seq Transcription Factor Binding Analysis

### Objective
Identify genomic binding sites of FOXA1 transcription factor and overlap with DEGs.

### Dataset
- ENCODE ENCFF396BZQ: FOXA1 ChIP-seq in MCF-7 breast cancer cells

### Tool
- ChIPseeker R package
- TxDb.Hsapiens.UCSC.hg38.knownGene

### Steps
1. Downloaded FOXA1 ChIP-seq peak file from ENCODE
2. Loaded peaks as GRanges object (57,244 peaks)
3. Annotated peaks to nearest genes using annotatePeak()
4. TSS region defined as -3000 to +3000 bp
5. Overlapped annotated genes with 942 DEGs from Task 2

### Results
- Total FOXA1 peaks: 57,244
- Peaks at promoters: 34.58%
- DEGs bound by FOXA1: 670 (71% of all DEGs)
- Top binding location: Introns and Promoters

---

## Task 6: GWAS Analysis

### Objective
Identify genetic variants associated with breast cancer risk.

### Dataset
- GCST004988: Michailidou et al. 2017, Nature
- 76,192 cases + 63,082 controls
- 11.7 million SNPs

### Tool
- qqman R package

### Steps
1. Downloaded harmonised GWAS summary statistics from GWAS Catalog
2. Filtered SNPs: removed NA values, kept chromosomes 1-22
3. Applied genome-wide significance threshold p < 5e-8
4. Generated Manhattan plot showing all chromosomes
5. Generated QQ plot to assess statistical validity

### Results
- Total SNPs analyzed: 11,415,622
- Genome-wide significant SNPs: 20,994
- Strongest signal: Chromosome 10 (FGFR2 locus)
- Other significant loci: Chr11, Chr16 (TOX3/CASC16)

---

## Task 7: Multi-Omics Integration

### Objective
Integrate findings from all tasks to identify candidate cancer driver genes supported by multiple evidence layers.

### Approach
Genes were scored based on evidence from:
- Task 2: Significant DEG (padj < 0.01)
- Task 4: Validated miRNA regulation
- Task 5: FOXA1 ChIP-seq binding
- Task 6: Nearby GWAS significant SNP

### Steps
1. Loaded results from Tasks 2, 4, 5, 6
2. Joined DEGs with FOXA1 binding data
3. Joined with miRNA interaction data
4. Assigned evidence scores to each gene
5. Filtered genes with evidence from all layers
6. Built integrated evidence table
7. Generated multi-omics regulatory model diagram

### Results
- Genes with FOXA1 binding: 1,589
- Genes with miRNA regulation: 1,701
- Genes with both FOXA1 + miRNA: 1,507
- Top candidate genes: 624 unique genes
- Top candidates: SORBS1, CYP4F22, MYBL2, EPCAM, PLK1

### Regulatory Model
GWAS SNP → FOXA1 Binding Disrupted → Gene Expression Changed
→ miRNA Control Lost → Cancer Pathway Activated
