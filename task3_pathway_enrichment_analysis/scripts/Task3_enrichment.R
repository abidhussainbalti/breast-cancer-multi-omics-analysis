install.packages("BiocManager")
BiocManager::install(c("DESeq2", "ggplot2", "pheatmap", "clusterProfiler",
"org.Hs.eg.db", "enrichplot", "ggrepel", "dplyr", "qqman", "biomaRt"))
packages <- c("DESeq2", "ggplot2", "pheatmap", "clusterProfiler",
"org.Hs.eg.db", "enrichplot", "ggrepel", "dplyr", "qqman", "biomaRt")
for(pkg in packages){
if(requireNamespace(pkg, quietly=TRUE)){
cat(pkg, "✓\n")
} else {
cat(pkg, "✗ NOT INSTALLED\n")
}
}
dir.create("C:/Users/abidhussain/Documents/Genomics_Project", recursive = TRUE)
dir.create("C:/Users/abidhussain/Documents/Genomics_Project/Task1", recursive = TRUE)
dir.create("C:/Users/abidhussain/Documents/Genomics_Project/Task2", recursive = TRUE)
dir.create("C:/Users/abidhussain/Documents/Genomics_Project/Task3", recursive = TRUE)
dir.create("C:/Users/abidhussain/Documents/Genomics_Project/Task4", recursive = TRUE)
dir.create("C:/Users/abidhussain/Documents/Genomics_Project/Task5", recursive = TRUE)
dir.create("C:/Users/abidhussain/Documents/Genomics_Project/Task6", recursive = TRUE)
dir.create("C:/Users/abidhussain/Documents/Genomics_Project/Task7", recursive = TRUE)
dir.create("D:/Genomics_Project/Task2", recursive = TRUE)
dir.create("D:/Genomics_Project/Task1", recursive = TRUE)
dir.create("D:/Genomics_Project/Task2", recursive = TRUE)
dir.create("D:/Genomics_Project/Task3", recursive = TRUE)
dir.create("D:/Genomics_Project/Task4", recursive = TRUE)
dir.create("D:/Genomics_Project/Task5", recursive = TRUE)
dir.create("D:/Genomics_Project/Task6", recursive = TRUE)
dir.create("D:/Genomics_Project/Task7", recursive = TRUE)
metadata <- data.frame(
sample = c("SRR15852399","SRR15852394","SRR15852395",
"SRR15852400","SRR15852407","SRR15852409",
"SRR15852426","SRR15852427","SRR15852429",
"SRR15852432","SRR15852435","SRR15852438"),
condition = c("tumor","tumor","tumor",
"tumor","tumor","tumor",
"normal","normal","normal",
"normal","normal","normal")
)
write.csv(metadata, "D:/Genomics_Project/Task2/metadata.csv", row.names=FALSE)
print(metadata)
# Load count matrix
counts <- read.table("D:/Genomics_Project/Task1/FeatureCounts_Mod.txt",
header=TRUE, row.names=1, sep="\t")
file.choose()
counts <- read.table("D:/Genomics_Project/Task1/FeatureCounts_Mod.txt.tabular",
header=TRUE, row.names=1, sep="\t")
dim(counts)
head(counts)
# Fix column names - keep only first SRR number
colnames(counts) <- gsub("_SRR.*", "", colnames(counts))
colnames(counts)
# Check if column names match metadata
all(colnames(counts) == metadata$sample)
library(DESeq2)
# Set reference level to normal
metadata$condition <- factor(metadata$condition, levels = c("normal", "tumor"))
# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts,
colData = metadata,
design = ~ condition)
# Filter low count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]
# Run DESeq2
dds <- DESeq(dds)
cat("DESeq2 analysis complete!\n")
cat("Number of genes after filtering:", nrow(dds), "\n")
# Extract results tumor vs normal
results <- results(dds, contrast = c("condition", "tumor", "normal"))
# Summary
summary(results)
library(dplyr)
# Convert to dataframe
results_df <- as.data.frame(results)
results_df$gene <- rownames(results_df)
# Filter significant genes padj < 0.01
sig_genes <- results_df %>%
filter(!is.na(padj) & padj < 0.01) %>%
arrange(padj)
cat("Total significant DEGs (padj < 0.01):", nrow(sig_genes), "\n")
cat("Upregulated:", sum(sig_genes$log2FoldChange > 0), "\n")
cat("Downregulated:", sum(sig_genes$log2FoldChange < 0), "\n")
# Save DEG results
write.csv(sig_genes, "D:/Genomics_Project/Task2/DEG_results.csv", row.names=FALSE)
library(ggplot2)
library(ggrepel)
# Prepare data
results_df$significance <- "Not Significant"
results_df$significance[results_df$padj < 0.01 & results_df$log2FoldChange > 1] <- "Upregulated"
results_df$significance[results_df$padj < 0.01 & results_df$log2FoldChange < -1] <- "Downregulated"
# Top 10 genes to label
top_genes <- results_df %>%
filter(!is.na(padj)) %>%
arrange(padj) %>% head(10)
# Plot
volcano <- ggplot(results_df, aes(x=log2FoldChange, y=-log10(padj), color=significance)) +
geom_point(alpha=0.6, size=1.5) +
scale_color_manual(values=c("Upregulated"="red", "Downregulated"="blue", "Not Significant"="grey")) +
geom_text_repel(data=top_genes, aes(label=gene), size=3, color="black") +
geom_vline(xintercept=c(-1,1), linetype="dashed") +
geom_hline(yintercept=-log10(0.01), linetype="dashed") +
labs(title="Volcano Plot: Tumor vs Normal Breast Cancer",
x="Log2 Fold Change", y="-Log10 Adjusted P-value") +
theme_bw()
ggsave("D:/Genomics_Project/Task2/volcano_plot.png", volcano, width=10, height=8, dpi=300)
print(volcano)
library(pheatmap)
# Get top 50 significant genes
top50 <- sig_genes %>% arrange(padj) %>% head(50)
# Extract normalized counts for top 50
norm_counts <- counts(dds, normalized=TRUE)
heat_data <- norm_counts[top50$gene, ]
# Log transform
heat_data <- log2(heat_data + 1)
# Annotation for columns
col_annotation <- data.frame(Condition = metadata$condition)
rownames(col_annotation) <- metadata$sample
# Plot heatmap
png("D:/Genomics_Project/Task2/heatmap.png", width=12, height=14, units="in", res=300)
pheatmap(heat_data,
annotation_col = col_annotation,
scale = "row",
show_rownames = TRUE,
show_colnames = TRUE,
main = "Top 50 DEGs Heatmap: Tumor vs Normal",
color = colorRampPalette(c("blue","white","red"))(100),
fontsize_row = 8)
dev.off()
cat("Heatmap saved!\n")
# Plot directly to RStudio viewer first
pheatmap(heat_data,
annotation_col = col_annotation,
scale = "row",
show_rownames = TRUE,
show_colnames = TRUE,
main = "Top 50 DEGs Heatmap: Tumor vs Normal",
color = colorRampPalette(c("blue","white","red"))(100),
fontsize_row = 8)
library(ggplot2)
library(reshape2)
# Melt data for ggplot
heat_melt <- melt(heat_data)
colnames(heat_melt) <- c("Gene", "Sample", "Expression")
# Add condition
heat_melt$Condition <- ifelse(heat_melt$Sample %in%
metadata$sample[metadata$condition=="tumor"],
"Tumor", "Normal")
# Plot
heatmap_gg <- ggplot(heat_melt, aes(x=Sample, y=Gene, fill=Expression)) +
geom_tile() +
scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
theme_bw() +
theme(axis.text.x = element_text(angle=90, hjust=1, size=7),
axis.text.y = element_text(size=7)) +
labs(title="Top 50 DEGs Heatmap: Tumor vs Normal")
ggsave("D:/Genomics_Project/Task2/heatmap.png", heatmap_gg,
width=14, height=12, dpi=300)
print(heatmap_gg)
cat("Heatmap saved!\n")
# Test if basic ggplot works
library(ggplot2)
test_plot <- ggplot(mtcars, aes(x=mpg, y=hp)) + geom_point()
print(test_plot)
# Add condition
# Check graphics device
dev.list()
dev.cur()
# Close all open devices
dev.off()
dev.off()
library(ggplot2)
test_plot <- ggplot(mtcars, aes(x=mpg, y=hp)) + geom_point()
print(test_plot)
# Now save heatmap properly
pheatmap(heat_data,
annotation_col = col_annotation,
scale = "row",
show_rownames = TRUE,
show_colnames = TRUE,
main = "Top 50 DEGs Heatmap: Tumor vs Normal",
color = colorRampPalette(c("blue","white","red"))(100),
fontsize_row = 8,
filename = "D:/Genomics_Project/Task2/heatmap.png",
width = 12,
height = 14)
cat("Heatmap saved!\n")
# Top 20 upregulated
top20_up <- sig_genes %>%
filter(log2FoldChange > 0) %>%
arrange(padj) %>%
head(20)
# Top 20 downregulated
top20_down <- sig_genes %>%
filter(log2FoldChange < 0) %>%
arrange(padj) %>%
head(20)
# Combine and save
top20_combined <- rbind(top20_up, top20_down)
write.csv(top20_combined, "D:/Genomics_Project/Task2/Top20_DEGs.csv", row.names=FALSE)
cat("Top 20 Upregulated:\n")
print(top20_up[,c("gene","log2FoldChange","padj")])
cat("\nTop 20 Downregulated:\n")
print(top20_down[,c("gene","log2FoldChange","padj")])
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
# Convert gene symbols to Entrez IDs (required by clusterProfiler)
gene_symbols <- sig_genes$gene
entrez_ids <- bitr(gene_symbols,
fromType = "SYMBOL",
toType = "ENTREZID",
OrgDb = org.Hs.eg.db)
cat("Genes successfully converted:", nrow(entrez_ids), "\n")
head(entrez_ids)
# GO Enrichment Analysis
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
OrgDb = org.Hs.eg.db,
ont = "BP",
pAdjustMethod = "BH",
pvalueCutoff = 0.05,
qvalueCutoff = 0.05,
readable = TRUE)
cat("Significant GO terms found:", nrow(go_results), "\n")
# Save results
go_df <- as.data.frame(go_results)
write.csv(go_df, "D:/Genomics_Project/Task3/GO_results.csv", row.names=FALSE)
head(go_df[,c("Description","GeneRatio","pvalue","p.adjust")], 10)
# GO Barplot
go_bar <- barplot(go_results,
showCategory=20,
title="Top 20 GO Biological Processes",
font.size=10)
ggsave("D:/Genomics_Project/Task3/GO_barplot.png",
go_bar, width=12, height=10, dpi=300)
cat("GO barplot saved!\n")
# KEGG Enrichment Analysis
kegg_results <- enrichKEGG(gene = entrez_ids$ENTREZID,
organism = "hsa",
pAdjustMethod = "BH",
pvalueCutoff = 0.05,
qvalueCutoff = 0.05)
cat("Significant KEGG pathways found:", nrow(kegg_results), "\n")
# Save results
kegg_df <- as.data.frame(kegg_results)
write.csv(kegg_df, "D:/Genomics_Project/Task3/KEGG_results.csv", row.names=FALSE)
head(kegg_df[,c("Description","GeneRatio","pvalue","p.adjust")], 10)
# KEGG Dotplot
kegg_dot <- dotplot(kegg_results,
showCategory=9,
title="KEGG Pathway Enrichment",
font.size=10)
ggsave("D:/Genomics_Project/Task3/KEGG_dotplot.png",
kegg_dot, width=10, height=8, dpi=300)
cat("KEGG dotplot saved!\n")
savehistory("D:/Genomics_Project/Task3/Task3_enrichment.R")
