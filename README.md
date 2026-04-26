# 🧬 Breast Cancer Multi-Omics Analysis

> A comprehensive multi-omics bioinformatics pipeline integrating RNA-seq, ChIP-seq, GWAS, and ncRNA data to identify key regulatory mechanisms driving breast cancer progression.

---

## 📋 Table of Contents
- [Project Summary](#project-summary)
- [Key Findings](#key-findings)
- [Disease Background](#disease-background)
- [Pipeline Overview](#pipeline-overview)
- [Repository Structure](#repository-structure)
- [Datasets Used](#datasets-used)
- [Tools and Software](#tools-and-software)
- [How To Run](#how-to-run)
- [Results Summary](#results-summary)
- [Authors](#authors)
- [References](#references)

---

## 📌 Project Summary

| Step | Task | Description |
|------|------|-------------|
| 1 | RNA-seq Alignment and Quantification | Processing raw sequencing reads and generating gene count matrix |
| 2 | Differential Gene Expression | Identifying significantly upregulated and downregulated genes |
| 3 | Functional Enrichment Analysis | Discovering biological processes and pathways associated with DEGs |
| 4 | ncRNA-miRNA Regulatory Network | Mapping regulatory interactions between miRNAs and target genes |
| 5 | ChIP-seq Binding Analysis | Identifying FOXA1 transcription factor binding sites |
| 6 | GWAS SNP Analysis | Detecting significant genetic variants associated with breast cancer |
| 7 | Multi-Omics Integration | Combining all datasets to identify candidate cancer driver genes |

---

## 🔬 Key Findings

- **942 differentially expressed genes** identified between tumor and normal breast tissue
- **Cell cycle pathway** most significantly enriched (KEGG hsa04110)
- **2,534 miRNAs** regulate DEGs including tumor suppressors let-7 and miR-34a
- **FOXA1 transcription factor** binds near 670 (71%) of all DEGs
- **20,994 genome-wide significant SNPs** identified in breast cancer GWAS
- **624 candidate genes** supported by all 4 omics evidence layers

---

## 🎗️ Disease Background

**Disease:** Breast Cancer — ER-positive breast carcinoma

**Dataset:** GSE183947
- 30 tumor + 30 normal tissue pairs
- Source: Guangzhou Medical University, China
- Platform: Illumina HiSeq 2000 (paired-end, 150bp)

---

## 🔄 Pipeline Overview

```text
+--------------------------------------------------+
|       RNA-seq Raw Reads (GSE183947)              |
|       12 samples: 6 tumor + 6 normal             |
+--------------------------------------------------+
                        |
                        v
+--------------------------------------------------+
|  Task 1: Alignment and Quantification            |
|  Tools : STAR v2.7.11b + featureCounts           |
|  Output: Count Matrix (59,251 genes x 12 samples)|
+--------------------------------------------------+
                        |
                        v
+--------------------------------------------------+
|  Task 2: Differential Gene Expression            |
|  Tool  : DESeq2                                  |
|  Output: 942 significant DEGs (padj < 0.01)      |
+--------------------------------------------------+
                        |
                        v
+--------------------------------------------------+
|  Task 3: Pathway Enrichment Analysis             |
|  Tool  : clusterProfiler                         |
|  Output: 304 GO Terms + 9 KEGG Pathways          |
+--------------------------------------------------+
                        |
                        v
+--------------------------------------------------+
|  Task 4: ncRNA Regulatory Network                |
|  Tool  : multiMiR                                |
|  Output: 90,452 miRNA-Gene Interactions          |
+--------------------------------------------------+
                        |
                        v
+--------------------------------------------------+
|  Task 5: ChIP-seq Binding Analysis               |
|  Tool  : ChIPseeker (FOXA1 ENCODE dataset)       |
|  Output: 670 DEGs bound by FOXA1                 |
+--------------------------------------------------+
                        |
                        v
+--------------------------------------------------+
|  Task 6: GWAS SNP Analysis                       |
|  Tool  : qqman                                   |
|  Output: 20,994 Significant Breast Cancer SNPs   |
+--------------------------------------------------+
                        |
                        v
+--------------------------------------------------+
|  Task 7: Multi-Omics Integration                 |
|  Output: 624 Candidate Cancer Driver Genes       |
+--------------------------------------------------+
```


---

## 📁 Repository Structure

| Folder | Contents |
|--------|----------|
| `datasets/` | DATA_SOURCES.md with all download links |
| `task1_rnaseq_alignment_breast_cancer/` | Galaxy workflow + FeatureCounts_Mod.txt |
| `task2_differential_gene_expression/` | DESeq2 script + DEGs + volcano plot + heatmap |
| `task3_pathway_enrichment_analysis/` | Enrichment script + GO/KEGG results + plots |
| `task4_ncrna_mirna_regulatory_network/` | ncRNA script + interaction tables + network plot |
| `task5_chipseq_foxa1_binding_analysis/` | ChIPseeker script + peak annotations + overlap |
| `task6_gwas_breast_cancer_snp_analysis/` | GWAS script + significant SNPs + Manhattan plot |
| `task7_multiomics_integration/` | Integration script + candidate genes + eQTL diagram |
| `report/` | Final project report (PDF) |

---

## 🗄️ Datasets Used

| Task | Dataset ID | Source | Description |
|------|-----------|--------|-------------|
| Task 1 | GSE183947 | NCBI GEO | Breast cancer RNA-seq, 12 samples (6 tumor + 6 normal) |
| Task 5 | ENCFF396BZQ | ENCODE Project | FOXA1 ChIP-seq in MCF-7 breast cancer cell line |
| Task 6 | GCST004988 | GWAS Catalog | Breast cancer GWAS — Michailidou et al. 2017, ~139,274 samples |

> Full download instructions: [datasets/DATA_SOURCES.md](datasets/DATA_SOURCES.md)

---

## 🛠️ Tools and Software

| Tool | Version | Purpose |
|------|---------|---------|
| Galaxy Project | usegalaxy.org | Cloud-based RNA-seq alignment pipeline |
| RNA STAR | 2.7.11b | Splice-aware read alignment |
| featureCounts | 2.0.6 | Gene-level read quantification |
| R | 4.3.3 | Statistical computing environment |
| DESeq2 | 1.42 | Differential gene expression analysis |
| clusterProfiler | 4.10 | GO and KEGG pathway enrichment |
| multiMiR | 1.24 | miRNA-target interaction queries |
| ChIPseeker | 1.38 | ChIP-seq peak annotation |
| qqman | 0.1.9 | Manhattan and QQ plot generation |

---

## ▶️ How To Run

### Step 1 — Install R Packages

```r
install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2", "clusterProfiler", "org.Hs.eg.db",
  "enrichplot", "ChIPseeker", "multiMiR",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "ggplot2", "pheatmap", "ggrepel", "dplyr", "qqman"
))
```

### Step 2 — Download Datasets

See [datasets/DATA_SOURCES.md](datasets/DATA_SOURCES.md) for all download links and instructions.

### Step 3 — Run Analysis Scripts

```r
source("task2_differential_gene_expression/scripts/Task2_DESeq2.R")
source("task3_pathway_enrichment_analysis/scripts/Task3_enrichment.R")
source("task4_ncrna_mirna_regulatory_network/scripts/Task4_ncRNA.R")
source("task5_chipseq_foxa1_binding_analysis/scripts/Task5_ChIPseek.R")
source("task6_gwas_breast_cancer_snp_analysis/scripts/Task6_GWAS.R")
source("task7_multiomics_integration/scripts/Task7_integration.R")
```

> **Note:** Task 1 was performed on Galaxy Project (usegalaxy.org).
> See `task1_rnaseq_alignment_breast_cancer/scripts/Galaxy_workflow_Task1.ga`

---

## 📊 Results Summary

| Task | Key Output File | Key Result |
|------|----------------|------------|
| Task 1 | FeatureCounts_Mod.txt | 59,251 genes quantified across 12 samples |
| Task 2 | DEG_results.csv | 942 significant DEGs (411 up, 531 down) |
| Task 3 | GO_results.csv, KEGG_results.csv | Cell cycle most enriched pathway |
| Task 4 | ncRNA_interactions.csv | 90,452 validated miRNA-gene interactions |
| Task 5 | DEG_ChIP_overlap.csv | 670 DEGs (71%) bound by FOXA1 |
| Task 6 | significant_SNPs.csv | 20,994 genome-wide significant SNPs |
| Task 7 | candidate_genes.csv | 624 multi-evidence candidate driver genes |

---

## 👤 Authors

**Abid Hussain**
NUST University — Genomics End Semester Project 2026

---

## 📚 References

1. Zhang Y et al. (2021) *Front Genet.* GSE183947. PMID: 35046993
2. Michailidou K et al. (2017) *Nature.* GCST004988. PMID: 29059683
3. ENCODE Project Consortium. FOXA1 ChIP-seq MCF-7. ENCFF396BZQ
4. Love MI et al. (2014) *Genome Biology.* DESeq2
5. Yu G et al. (2012) *OMICS.* clusterProfiler
