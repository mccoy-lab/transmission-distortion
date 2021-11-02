library(data.table)
library(tidyverse)
library(cowplot)
library(ggpointdensity)
library(plyr)

#dt_p <- fread("~/Dropbox/projects/sperm_seq/nc17ab/csv_out/nc17ab_chr6_pval.csv")
dt_p <- fread("/Users/saracarioscia/mccoy-lab/transmission-distortion/raw_data_for_tests/nc17ab_chr6_pval.csv")

#dt <- fread("~/Downloads/nc17ab_goodcellsreplicatebcs_filteredhetsnps_6.cellsbyrow.txt")
dt <- fread("/Users/saracarioscia/mccoy-lab/transmission-distortion/raw_data_for_tests/nc17ab_goodcellsreplicatebcs_filteredhetsnps_6.cellsbyrow.txt") 

dt_mat <- dt %>%
  pivot_wider(., names_from = "cell", values_from = "gt") %>%
  arrange(., pos) %>%
  as.data.table()

# keep only the SNPs that are in both tables because some got filtered out before the p-values were calculated (127,030 vs. 122,461) 
dt_mat <- dt_mat[which(dt_mat$pos %in% dt_p$genomic_position),]
dt_mat_hold <- dt_mat

dt_p[, tot_obs := rowSums(!is.na(dt_mat))]

a <- ggplot(dt_p[genomic_position > 164014420], aes(x = genomic_position, y = -log10(pval))) +
  geom_point()

b <- ggplot(dt_p[genomic_position > 164014420], aes(x = genomic_position, y = tot_obs)) +
  geom_pointdensity() +
  theme(legend.position = "none")

ggplot(dt_p, aes(x = (pval < 1e-8), y = tot_obs)) +
  geom_boxplot() +
  ylab("Number of sperm containing SNP") + 
  xlab("SNP p-value < 1e-8") +
  theme(legend.position = "none")



dt_mat <- dt[pos > 168014420 & pos < 168177306] %>%
  pivot_wider(., names_from = "cell", values_from = "gt") %>%
  arrange(., pos) %>%
  as.data.table()

rowSums(dt_mat[, 2:ncol(dt_mat)] == 1, na.rm = TRUE)
rowSums(dt_mat[, 2:ncol(dt_mat)] == 0, na.rm = TRUE)

allele_counts <- data.table(x = rowSums(dt_mat[, 2:ncol(dt_mat)] == 1, na.rm = TRUE), 
                            y = rowSums(dt_mat[, 2:ncol(dt_mat)] == 0, na.rm = TRUE))

ggplot(data = allele_counts, aes(x = x, y = y)) +
  geom_point() +
  xlim(0, 100) +
  ylim(0, 100) +
  geom_abline(slope = 1) +
  geom_abline(slope = 2) +
  geom_abline(slope = 0.5)

dt_mat <- dt[pos > 158014420 & pos < 158177306] %>%
  pivot_wider(., names_from = "cell", values_from = "gt") %>%
  arrange(., pos) %>%
  as.data.table()

rowSums(dt_mat[, 2:ncol(dt_mat)] == 1, na.rm = TRUE)
rowSums(dt_mat[, 2:ncol(dt_mat)] == 0, na.rm = TRUE)

control_allele_counts <- data.table(x = rowSums(dt_mat[, 2:ncol(dt_mat)] == 1, na.rm = TRUE), 
                                    y = rowSums(dt_mat[, 2:ncol(dt_mat)] == 0, na.rm = TRUE))

ggplot(data = control_allele_counts, aes(x = x, y = y)) +
  geom_point() +
  xlim(0, 100) +
  ylim(0, 100) +
  geom_abline(slope = 1)