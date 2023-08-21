# Code to investigate potential duplications in the putative peak region 
# dt_p is the pvalues resulting from the TD test, to narrow down a region of the chromosome to investigate
# dt is the raw SNP and sperm data 

#

library(data.table)
library(tidyverse)
library(cowplot)
library(ggpointdensity)
library(plyr)

#dt_p <- fread("~/Dropbox/projects/sperm_seq/nc17ab/csv_out/nc17ab_chr6_pval.csv")
#dt_p <- fread("/Users/saracarioscia/mccoy-lab/transmission-distortion/raw_data_for_tests/nc17ab_chr10_pval.csv")
#dt_p <- fread("/Users/saracarioscia/mccoy-lab/transmission-distortion/raw_data_for_tests/nc12ab_chr10_pval.csv")
#dt_p <- fread("/Users/saracarioscia/mccoy-lab/transmission-distortion/raw_data_for_tests/nc27aboct17_chr22_pval.csv")
#dt_p <- fread("/Users/saracarioscia/mccoy-lab/transmission-distortion/raw_data_for_tests/nc17ab_chr6_pval.csv")

#dt <- fread("~/Downloads/nc17ab_goodcellsreplicatebcs_filteredhetsnps_6.cellsbyrow.txt")
#dt <- fread("/Users/saracarioscia/mccoy-lab/transmission-distortion/raw_data_for_tests/nc17ab_goodcellsreplicatebcs_filteredhetsnps_10.cellsbyrow.txt")
#dt <- fread("/Users/saracarioscia/mccoy-lab/transmission-distortion/raw_data_for_tests/nc12ab_goodcellsreplicatebcs_filteredhetsnps_10.cellsbyrow.txt")
#dt <- fread("/Users/saracarioscia/mccoy-lab/transmission-distortion/raw_data_for_tests/nc27aboct17_goodcellsreplicatebcs_filteredhetsnps_22.cellsbyrow.txt")
#dt <- fread("/Users/saracarioscia/mccoy-lab/transmission-distortion/raw_data_for_tests/nc17ab_goodcellsreplicatebcs_filteredhetsnps_6.cellsbyrow.txt") 

dt <- fread("/Users/saracarioscia/mccoy-lab/transmission-distortion/raw_data_for_tests/nc18ab_goodcellsreplicatebcs_filteredhetsnps_2.cellsbyrow.txt") 

dt_mat <- dt %>%
  pivot_wider(., names_from = "cell", values_from = "gt") %>%
  arrange(., pos) %>%
  as.data.table()
# keep only the SNPs that are in both tables 
dt_mat <- dt_mat[which(dt_mat$pos %in% dt_p$genomic_position),]
dt_mat_hold <- dt_mat

dt_p[, tot_obs := rowSums(!is.na(dt_mat))]

# a <- ggplot(dt_p[genomic_position > 164014420], aes(x = genomic_position, y = -log10(pval))) +
#   geom_point()
# 
# b <- ggplot(dt_p[genomic_position > 164014420], aes(x = genomic_position, y = tot_obs)) +
#   geom_pointdensity() +
#   theme(legend.position = "none")
a <- ggplot(dt_p[genomic_position > 164014420], aes(x = genomic_position, y = -log10(pval))) +
  geom_point()

b <- ggplot(dt_p[genomic_position > 164014420], aes(x = genomic_position, y = tot_obs)) +
  geom_pointdensity() +
  theme(legend.position = "none")


ggplot(dt_p, aes(x = (pval < 1e-8), y = tot_obs)) +
  geom_boxplot() +
  theme(legend.position = "none")


# dt_mat <- dt[pos > 168014420 & pos < 168177306] %>%
#   pivot_wider(., names_from = "cell", values_from = "gt") %>%
#   arrange(., pos) %>%
#   as.data.table()
# dt_mat <- dt[pos > 78284602 & pos < 78484602] %>%
#   pivot_wider(., names_from = "cell", values_from = "gt") %>%
#   arrange(., pos) %>%
#   as.data.table()
# 167900550 to 168329339
dt_mat <- dt[pos > 167900550 & pos < 168329339] %>%
  pivot_wider(., names_from = "cell", values_from = "gt") %>%
  arrange(., pos) %>%
  as.data.table()
# dt_mat <- dt[pos > 25311350 & pos < 25508290] %>%
#   pivot_wider(., names_from = "cell", values_from = "gt") %>%
#   arrange(., pos) %>%
#   as.data.table()

rowSums(dt_mat[, 2:ncol(dt_mat)] == 1, na.rm = TRUE)
rowSums(dt_mat[, 2:ncol(dt_mat)] == 0, na.rm = TRUE)

allele_counts <- data.table(x = rowSums(dt_mat[, 2:ncol(dt_mat)] == 1, na.rm = TRUE), 
                            y = rowSums(dt_mat[, 2:ncol(dt_mat)] == 0, na.rm = TRUE))
#allele_counts_nc17_chr10 <- allele_counts
allele_counts_nc12_chr10 <- allele_counts
#allele_counts_nc27_chr22 <- allele_counts

# check whether any alleles are closer to the line slope = 1 (expected), slope = 2 (potential over-representation), slope = 0.5 (potential under-representation)
# some variation is expected, but a tight clustering to the 0.5 and 2 lines suggests a segmental duplication 
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