
# Task 5: ChIP-seq Analysis
# FOXA1 ChIP-seq data: ENCODE ENCFF396BZQ (MCF-7 breast cancer cells)

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)

# Load data
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peaks <- readPeakFile("D:/Genomics_Project/Task5/FOXA1_peaks.bed.gz", as="GRanges")

# Annotate peaks
peakAnno <- annotatePeak(peaks, tssRegion=c(-3000,3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

# Save plots
png("D:/Genomics_Project/Task5/peak_annotation_pie.png", width=10, height=8, units="in", res=300)
plotAnnoPie(peakAnno, main="FOXA1 ChIP-seq Peak Annotation")
dev.off()

# Extract annotated genes
anno_df <- as.data.frame(peakAnno)
chip_genes <- anno_df %>% filter(!is.na(SYMBOL)) %>%
  dplyr::select(seqnames, start, end, annotation, distanceToTSS, SYMBOL, GENENAME) %>%
  distinct()
write.csv(chip_genes, "D:/Genomics_Project/Task5/annotated_peaks.csv", row.names=FALSE)

# Overlap with DEGs
sig_genes <- read.csv("D:/Genomics_Project/Task2/DEG_results.csv")
overlap <- chip_genes %>% filter(SYMBOL %in% sig_genes$gene) %>%
  dplyr::select(SYMBOL, annotation, distanceToTSS) %>% distinct()
write.csv(overlap, "D:/Genomics_Project/Task5/DEG_ChIP_overlap.csv", row.names=FALSE)

