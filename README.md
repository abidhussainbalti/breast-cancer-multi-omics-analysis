# Breast Cancer Multi-Omics Analysis

A comprehensive multi-omics analysis pipeline to investigate breast cancer by integrating RNA-seq, ChIP-seq, GWAS, and regulatory network data. This project identifies key genes, pathways, and regulatory mechanisms involved in breast cancer progression.

## Project Summary

This repository presents a step-by-step computational workflow covering:
```
RNA-seq alignment and quantification
Differential gene expression analysis
Functional enrichment (GO & KEGG)
ncRNA–miRNA regulatory network construction
ChIP-seq (FOXA1) binding analysis
GWAS-based SNP analysis
Multi-omics data integration
```
## Key Findings

- **942 differentially expressed genes** identified between tumor and normal breast tissue
- **Cell cycle pathway** most significantly enriched (KEGG hsa04110)
- **2,534 miRNAs** regulate DEGs including tumor suppressors let-7 and miR-34a
- **FOXA1 transcription factor** binds near 670 (71%) of all DEGs
- **20,994 genome-wide significant SNPs** identified in breast cancer GWAS
- **624 candidate genes** supported by all 4 omics evidence layers

## Disease

**Breast Cancer** — ER-positive breast carcinoma
Dataset: GSE183947 (30 tumor + 30 normal pairs, Guangzhou Medical University)

## Pipeline Overview
```
RNA-seq Raw Reads (GSE183947)
↓ Task 1: STAR Alignment + featureCounts
Count Matrix (59,251 genes x 12 samples)
↓ Task 2: DESeq2
942 Differentially Expressed Genes
↓ Task 3: clusterProfiler
304 GO Terms + 9 KEGG Pathways
↓ Task 4: multiMiR
90,452 miRNA-Gene Interactions
↓ Task 5: ChIPseeker (FOXA1 ChIP-seq)
670 DEGs bound by FOXA1
↓ Task 6: GWAS Analysis
20,994 Significant Breast Cancer SNPs
↓ Task 7: Multi-Omics Integration
624 Candidate Cancer Driver Genes
```
## Repository Structure
```breast-cancer-multi-omics-analysis/
├── README.md
├── METHODOLOGY.md
├── .gitignore
├── datasets/
│   └── DATA_SOURCES.md
├── task1_rnaseq_alignment_breast_cancer/
│   ├── scripts/
│   │   └── Galaxy_workflow_Task1.ga
│   └── results/
│       └── FeatureCounts_Mod.txt
├── task2_differential_gene_expression/
│   ├── scripts/
│   │   └── Task2_DESeq2.R
│   └── results/
│       ├── DEG_results.csv
│       ├── Top20_DEGs.csv
│       ├── volcano_plot.png
│       └── heatmap.png
├── task3_pathway_enrichment_analysis/
│   ├── scripts/
│   │   └── Task3_enrichment.R
│   └── results/
│       ├── GO_results.csv
│       ├── KEGG_results.csv
│       ├── GO_barplot.png
│       └── KEGG_dotplot.png
├── task4_ncrna_mirna_regulatory_network/
│   ├── scripts/
│   │   └── Task4_ncRNA.R
│   └── results/
│       ├── ncRNA_interactions.csv
│       ├── ncRNA_list.csv
│       └── network.png
├── task5_chipseq_foxa1_binding_analysis/
│   ├── scripts/
│   │   └── Task5_ChIPseq.R
│   └── results/
│       ├── annotated_peaks.csv
│       ├── DEG_ChIP_overlap.csv
│       ├── peak_annotation_pie.png
│       ├── peak_annotation_bar.png
│       └── overlap_annotation.png
├── task6_gwas_breast_cancer_snp_analysis/
│   ├── scripts/
│   │   └── Task6_GWAS.R
│   └── results/
│       ├── significant_SNPs.csv
│       ├── manhattan_plot.png
│       └── QQ_plot.png
├── task7_multiomics_integration/
│   ├── scripts/
│   │   └── Task7_integration.R
│   └── results/
│       ├── candidate_genes.csv
│       ├── integrated_evidence_table.csv
│       └── eQTL_diagram.png
└── report```

## Datasets Used
```

| Task | Dataset | Source | Description |
|------|---------|--------|-------------|
| Task 1 | GSE183947 | NCBI GEO | Breast cancer RNA-seq, 12 samples |
| Task 5 | ENCFF396BZQ | ENCODE | FOXA1 ChIP-seq in MCF-7 cells |
| Task 6 | GCST004988 | GWAS Catalog | Breast cancer GWAS, 139,274 samples |
```
Full download instructions: [datasets/DATA_SOURCES.md](datasets/DATA_SOURCES.md)

## Tools and Software
```

| Tool | Version | Purpose |
|------|---------|---------|
| Galaxy Project | usegalaxy.org | RNA-seq alignment pipeline |
| RNA STAR | 2.7.11b | Read alignment |
| featureCounts | 2.0.6 | Read quantification |
| R | 4.3.3 | Statistical analysis |
| DESeq2 | 1.42 | Differential expression |
| clusterProfiler | 4.10 | Pathway enrichment |
| multiMiR | 1.24 | miRNA interactions |
| ChIPseeker | 1.38 | ChIP-seq annotation |
| qqman | 0.1.9 | GWAS visualization |
```
## How To Run

### Install R Packages

```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db",
                       "enrichplot", "ChIPseeker", "multiMiR",
                       "TxDb.Hsapiens.UCSC.hg38.knownGene",
                       "ggplot2", "pheatmap", "ggrepel",
                       "dplyr", "qqman"))
```

### Run Analysis

```r
source("task2_differential_gene_expression/scripts/Task2_DESeq2.R")
source("task3_pathway_enrichment_analysis/scripts/Task3_enrichment.R")
source("task4_ncrna_mirna_regulatory_network/scripts/Task4_ncRNA.R")
source("task5_chipseq_foxa1_binding_analysis/scripts/Task5_ChIPseq.R")
source("task6_gwas_breast_cancer_snp_analysis/scripts/Task6_GWAS.R")
source("task7_multiomics_integration/scripts/Task7_integration.R")
```

## Results Summary
```

| Task | Output | Key Result |
|------|--------|------------|
| Task 1 | FeatureCounts_Mod.txt | 59,251 genes quantified |
| Task 2 | DEG_results.csv | 942 significant DEGs |
| Task 3 | GO_results.csv, KEGG_results.csv | Cell cycle most enriched |
| Task 4 | ncRNA_interactions.csv | 90,452 miRNA interactions |
| Task 5 | DEG_ChIP_overlap.csv | 670 FOXA1-bound DEGs |
| Task 6 | significant_SNPs.csv | 20,994 significant SNPs |
| Task 7 | candidate_genes.csv | 624 candidate genes |
```
## Authors

- Abid Hussain — NUST University, Genomics End Semester Project 2026

## References

1. Zhang Y et al. (2021) Front Genet. GSE183947. PMID: 35046993
2. Michailidou K et al. (2017) Nature. GCST004988. PMID: 29059683
3. ENCODE Project Consortium. ENCFF396BZQ
4. Love MI et al. (2014) DESeq2. Genome Biology
5. Yu G et al. (2012) clusterProfiler. OMICS
