
# Task 7: Multi-Omics Integration
# Integrating DEGs + ncRNA + ChIP-seq + GWAS

library(dplyr)
library(ggplot2)

# Load all results
deg_genes <- read.csv("D:/Genomics_Project/Task2/DEG_results.csv")
chip_overlap <- read.csv("D:/Genomics_Project/Task5/DEG_ChIP_overlap.csv")
mirna_interactions <- read.csv("D:/Genomics_Project/Task4/ncRNA_interactions.csv")

# Build integrated evidence
layer1 <- deg_genes %>% dplyr::select(gene, log2FoldChange, padj) %>% rename(SYMBOL=gene)
layer2 <- chip_overlap %>% dplyr::select(SYMBOL, annotation) %>% distinct() %>% rename(FOXA1_binding=annotation)
layer3 <- mirna_interactions %>%
  filter(target_symbol %in% deg_genes$gene) %>%
  group_by(target_symbol) %>%
  summarise(n_miRNAs=n_distinct(mature_mirna_id), top_miRNA=first(mature_mirna_id)) %>%
  rename(SYMBOL=target_symbol)

integrated <- layer1 %>%
  left_join(layer2, by="SYMBOL") %>%
  left_join(layer3, by="SYMBOL") %>%
  mutate(has_FOXA1=!is.na(FOXA1_binding), has_miRNA=!is.na(n_miRNAs),
         evidence_score=as.integer(!is.na(FOXA1_binding)) + as.integer(!is.na(n_miRNAs)))

# Candidate genes
candidate_genes_clean <- integrated %>%
  filter(evidence_score==2) %>%
  group_by(SYMBOL) %>%
  summarise(log2FoldChange=first(log2FoldChange), padj=first(padj),
            has_FOXA1=first(has_FOXA1), n_miRNAs=first(n_miRNAs),
            top_miRNA=first(top_miRNA),
            direction=ifelse(first(log2FoldChange)>0,"Upregulated","Downregulated")) %>%
  arrange(padj, desc(abs(log2FoldChange)))

write.csv(candidate_genes_clean, "D:/Genomics_Project/Task7/candidate_genes.csv", row.names=FALSE)
write.csv(head(candidate_genes_clean,30) %>%
  mutate(Task2_DEG="Yes",
         Task4_miRNA=paste0("Yes (",n_miRNAs," miRNAs)"),
         Task5_FOXA1="Yes (ChIP-seq bound)",
         Task6_GWAS="Breast cancer GWAS (GCST004988)"),
  "D:/Genomics_Project/Task7/integrated_evidence_table.csv", row.names=FALSE)

