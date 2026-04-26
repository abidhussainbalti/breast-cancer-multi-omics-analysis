# Dataset Sources and Download Instructions

## Task 1 — RNA-seq Dataset
**Dataset:** GSE183947  
**Title:** Identification of five cytotoxicity-related genes involved in the progression of breast cancer  
**Source:** NCBI GEO  
**Link:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183947  
**Samples:** 12 samples (6 tumor + 6 normal)  
**SRR Accessions:**
- Tumor: SRR15852399, SRR15852394, SRR15852395, SRR15852400, SRR15852407, SRR15852409
- Normal: SRR15852426, SRR15852427, SRR15852429, SRR15852432, SRR15852435, SRR15852438

**Download via SRA:**
```bash
# Download each sample using fasterq-dump
fasterq-dump SRR15852399 --split-files -O ./raw_data/
```

**Or download directly via AWS:**
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15852399/SRR15852399

**Reference Genome (hg38):**
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

**Gene Annotation (GTF):**
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz

---

## Task 5 — ChIP-seq Dataset
**Dataset:** ENCODE ENCFF396BZQ  
**Title:** FOXA1 ChIP-seq in MCF-7 breast cancer cells  
**Source:** ENCODE Project  
**Link:** https://www.encodeproject.org/files/ENCFF396BZQ/  
**Download:**
```bash
wget https://www.encodeproject.org/files/ENCFF396BZQ/@@download/ENCFF396BZQ.bed.gz \
     -O FOXA1_peaks.bed.gz
```

---

## Task 6 — GWAS Dataset
**Dataset:** GCST004988  
**Title:** Genome-wide association study of breast cancer (Michailidou et al. 2017, Nature)  
**Source:** GWAS Catalog  
**Link:** https://www.ebi.ac.uk/gwas/studies/GCST004988  
**Samples:** 76,192 cases + 63,082 controls  
**Download:**
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004988/harmonised/29059683-GCST004988-EFO_0000305-build37.f.tsv.gz

---

## Tasks 2, 3, 4 — Databases Used
| Database | Purpose | Link |
|---|---|---|
| DESeq2 | Differential expression | https://bioconductor.org/packages/DESeq2 |
| clusterProfiler | GO/KEGG enrichment | https://bioconductor.org/packages/clusterProfiler |
| multiMiR | miRNA interactions | https://bioconductor.org/packages/multiMiR |
| miRTarBase | Validated miRNA targets | https://mirtarbase.cuhk.edu.cn |
| KEGG | Pathway database | https://www.genome.jp/kegg |
| org.Hs.eg.db | Human gene annotation | https://bioconductor.org/packages/org.Hs.eg.db |
