
# Task 2: Differential Gene Expression Analysis
# Dataset: GSE183947 - Breast Cancer RNA-seq

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)

# Load data
counts <- read.table("D:/Genomics_Project/Task1/FeatureCounts_Mod.txt.tabular", 
                     header=TRUE, row.names=1, sep="	")
colnames(counts) <- gsub("_SRR.*", "", colnames(counts))

# Metadata
metadata <- read.csv("D:/Genomics_Project/Task2/metadata.csv")
metadata$condition <- factor(metadata$condition, levels=c("normal","tumor"))

# DESeq2
dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

# Results
results <- results(dds, contrast=c("condition","tumor","normal"))
results_df <- as.data.frame(results)
results_df$gene <- rownames(results_df)

# Filter DEGs
sig_genes <- results_df %>% filter(!is.na(padj) & padj < 0.01) %>% arrange(padj)
write.csv(sig_genes, "D:/Genomics_Project/Task2/DEG_results.csv", row.names=FALSE)

# Volcano plot
results_df$significance <- "Not Significant"
results_df$significance[results_df$padj < 0.01 & results_df$log2FoldChange > 1] <- "Upregulated"
results_df$significance[results_df$padj < 0.01 & results_df$log2FoldChange < -1] <- "Downregulated"
top_genes <- results_df %>% filter(!is.na(padj)) %>% arrange(padj) %>% head(10)
volcano <- ggplot(results_df, aes(x=log2FoldChange, y=-log10(padj), color=significance)) +
  geom_point(alpha=0.6, size=1.5) +
  scale_color_manual(values=c("Upregulated"="red","Downregulated"="blue","Not Significant"="grey")) +
  geom_text_repel(data=top_genes, aes(label=gene), size=3, color="black") +
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.01), linetype="dashed") +
  labs(title="Volcano Plot: Tumor vs Normal Breast Cancer",
       x="Log2 Fold Change", y="-Log10 Adjusted P-value") +
  theme_bw()
ggsave("D:/Genomics_Project/Task2/volcano_plot.png", volcano, width=10, height=8, dpi=300)

# Heatmap
norm_counts <- counts(dds, normalized=TRUE)
top50 <- sig_genes %>% arrange(padj) %>% head(50)
heat_data <- log2(norm_counts[top50$gene, ] + 1)
col_annotation <- data.frame(Condition=metadata$condition)
rownames(col_annotation) <- metadata$sample
pheatmap(heat_data, annotation_col=col_annotation, scale="row",
         main="Top 50 DEGs Heatmap: Tumor vs Normal",
         color=colorRampPalette(c("blue","white","red"))(100),
         fontsize_row=8,
         filename="D:/Genomics_Project/Task2/heatmap.png",
         width=12, height=14)

# Top 20 DEGs
top20_up <- sig_genes %>% filter(log2FoldChange > 0) %>% arrange(padj) %>% head(20)
top20_down <- sig_genes %>% filter(log2FoldChange < 0) %>% arrange(padj) %>% head(20)
write.csv(rbind(top20_up, top20_down), "D:/Genomics_Project/Task2/Top20_DEGs.csv", row.names=FALSE)

