library(data.table)
library(tidyverse)

acc_5k <- fread("~/Downloads/acc_5k_snp.csv", header = TRUE) %>%
  pivot_longer(-V1, names_to = "coverage") %>%
  as.data.table() %>%
  setnames(., c("n_gametes", "coverage", "value")) %>%
  .[, n_snps := 5000] %>%
  .[, metric := "accuracy"]

acc_30k <- fread("~/Downloads/acc_30k_snp.csv", header = TRUE) %>%
  pivot_longer(-V1, names_to = "coverage") %>%
  as.data.table() %>%
  setnames(., c("n_gametes", "coverage", "value")) %>%
  .[, n_snps := 30000] %>%
  .[, metric := "accuracy"]

acc_100k <- fread("~/Downloads/acc_100k_snp.csv", header = TRUE) %>%
  pivot_longer(-V1, names_to = "coverage") %>%
  as.data.table() %>%
  setnames(., c("n_gametes", "coverage", "value")) %>%
  .[, n_snps := 100000] %>%
  .[, metric := "accuracy"]

com_5k <- fread("~/Downloads/com_5k_snp.csv", header = TRUE) %>%
  pivot_longer(-V1, names_to = "coverage") %>%
  as.data.table() %>%
  setnames(., c("n_gametes", "coverage", "value")) %>%
  .[, n_snps := 5000] %>%
  .[, metric := "completeness"]

com_30k <- fread("~/Downloads/com_30k_snp.csv", header = TRUE) %>%
  pivot_longer(-V1, names_to = "coverage") %>%
  as.data.table() %>%
  setnames(., c("n_gametes", "coverage", "value")) %>%
  .[, n_snps := 30000] %>%
  .[, metric := "completeness"]

com_100k <- fread("~/Downloads/com_100k_snp.csv", header = TRUE) %>%
  pivot_longer(-V1, names_to = "coverage") %>%
  as.data.table() %>%
  setnames(., c("n_gametes", "coverage", "value")) %>%
  .[, n_snps := 100000] %>%
  .[, metric := "completeness"]

snp_names <- as_labeller(c(`5000` = "5000 SNPs", `30000` = "30000 SNPs",`1e+05` = "100000 SNPs"))
metric_names <- as_labeller(c("accuracy" = "Phasing\nAccuracy (%)", "completeness" = "Phasing\nCompleteness (%)"))

ggplot(data = rbind(acc_5k, acc_30k, acc_100k, com_5k, com_30k, com_100k),
       aes(x = factor(coverage), y = factor(n_gametes), fill = value)) +
  geom_tile() + ggtitle('Donor Haplotype Phasing') +
  facet_grid(metric ~ n_snps, labeller = labeller(metric = metric_names, n_snps = snp_names)) +
  theme(panel.background = element_blank(), panel.grid = element_blank()) +
  scale_fill_viridis_c() + xlab("Coverage (x)") + ylab ("Number of gametes") + theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  theme(text = element_text(size=25))