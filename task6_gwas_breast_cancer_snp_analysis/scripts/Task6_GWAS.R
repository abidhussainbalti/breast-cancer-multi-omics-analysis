
# Task 6: GWAS Analysis
# Dataset: GCST004988 - Michailidou et al. 2017, Nature
# Breast cancer GWAS - 76,192 cases, 63,082 controls

library(qqman)
library(dplyr)

# Load GWAS summary statistics
gwas_raw <- read.table(
  "D:/Genomics_Project/Task6/29059683-GCST004988-EFO_0000305-build37.f.tsv.gz",
  header=TRUE, sep="	", quote="", fill=TRUE)

# Clean data
gwas_clean <- gwas_raw %>%
  dplyr::select(variant_id, chromosome, base_pair_location, p_value) %>%
  rename(SNP=variant_id, CHR=chromosome, BP=base_pair_location, P=p_value) %>%
  filter(!is.na(P) & !is.na(CHR) & !is.na(BP)) %>%
  filter(P > 0 & P <= 1) %>%
  filter(CHR %in% 1:22) %>%
  mutate(CHR=as.integer(CHR), BP=as.integer(BP))

# Save significant SNPs
sig_snps <- gwas_clean %>% filter(P < 5e-8) %>% arrange(P)
write.csv(sig_snps, "D:/Genomics_Project/Task6/significant_SNPs.csv", row.names=FALSE)

# Subsample for plotting
set.seed(42)
gwas_plot <- gwas_clean %>%
  filter(P < 0.01) %>%
  bind_rows(gwas_clean %>% filter(P >= 0.01) %>% sample_n(500000))

# Manhattan plot
png("D:/Genomics_Project/Task6/manhattan_plot.png", width=14, height=7, units="in", res=300)
manhattan(gwas_plot, chr="CHR", bp="BP", snp="SNP", p="P",
          main="Manhattan Plot - Breast Cancer GWAS (Michailidou et al. 2017)",
          suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8),
          col=c("steelblue","orange"), cex=0.4)
dev.off()

# QQ plot
png("D:/Genomics_Project/Task6/QQ_plot.png", width=8, height=8, units="in", res=300)
qq(gwas_plot$P, main="QQ Plot - Breast Cancer GWAS (Michailidou et al. 2017)")
dev.off()

