
# Task 4: ncRNA Regulatory Analysis
# Dataset: GSE183947 DEGs queried against multiMiR (miRTarBase, miRecords, TarBase)

library(multiMiR)
library(dplyr)
library(ggplot2)

# Load DEGs from Task 2
sig_genes <- read.csv("D:/Genomics_Project/Task2/DEG_results.csv")
all_degs <- sig_genes$gene

# Query validated miRNA-target interactions
mirna_results <- get.multimir(org = "hsa", target = all_degs,
                               table = "validated", summary = TRUE)
mirna_data <- mirna_results@data

# Filter validated interactions
mirna_filtered <- mirna_data %>%
  filter(type == "validated") %>%
  dplyr::select(database, mature_mirna_id, target_symbol, experiment, support_type) %>%
  distinct() %>%
  arrange(target_symbol, mature_mirna_id)

# Save outputs
write.csv(mirna_filtered, "D:/Genomics_Project/Task4/ncRNA_interactions.csv", row.names=FALSE)

mirna_list <- mirna_filtered %>%
  group_by(mature_mirna_id) %>%
  summarise(n_targets = n_distinct(target_symbol)) %>%
  arrange(desc(n_targets))
write.csv(mirna_list, "D:/Genomics_Project/Task4/ncRNA_list.csv", row.names=FALSE)

