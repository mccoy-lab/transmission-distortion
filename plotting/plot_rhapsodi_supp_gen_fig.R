library(itertools)
library(tidyverse)
library(data.table)
library(ggplot2)
library(gridExtra)

nsnps <- c(5000, 30000, 100000)
nsnps_iter <- as.list(itertools::enumerate(nsnps))
ngams <- c(3, 15, 50, 150, 500, 1000, 2500, 5000)
ngams_iter <- as.list(itertools::enumerate(ngams))
covs <- c(0.001, 0.01, 0.1, 0.223, 0.357, 0.511, 0.693, 1.204, 2.303)
covs_iter <- as.list(itertools::enumerate(covs))
num_vals <- length(nsnps) * length(ngams)  * length(covs) * 7
rsds <- c(42, 357, 1848)
gen_seqerror <- 0.005
mod_seqerror <- 0.005
gen_avgrecomb <- 1
mod_avgrecomb <- 1

data_mat <- array(rep(NA, num_vals), dim = c(length(nsnps), length(ngams), length(covs), 7))
last_dim_metric <- list("phasing ser" = 1, 
                        "phasing lhs" = 2, 
                        "gam_imputation ser" = 3,
                        "recomb precision" = 4, 
                        "recomb accuracy" = 5, 
                        "recomb specificity" = 6, 
                        "recomb fpr" = 7)
dir_base <- "/home/kweave23/gamete_data/gen_model_results_noDNM/"
true_snp_file <- "true_nsnp_gen_model_noDNM.Rdata"

load(paste0(dir_base, true_snp_file))

for (i in last_dim_metric){
  metricOI <- names(last_dim_metric[i])
  metricOI_top <- strsplit(metricOI, " ")[[1]][1]
  metricOI_bot <- strsplit(metricOI, " ")[[1]][2]
  for (nsnp in nsnps_iter){
    for (ngam in ngams_iter){
      for (cov in covs_iter){
        dir_data <- paste0("g", ngam$value, "_s", as.character(as.integer(nsnp$value)), "_c", cov$value, "_se", gen_seqerror, "_r", gen_avgrecomb, "/")
        data_vec <- c()
        for (rsd in rsds){
          data_file <- paste0("assess_out_rs_", rsd, ".Rdata")
          fileOI <- paste0(dir_base, dir_data, data_file)
          if (file_test("-f", fileOI)){
            load(fileOI)
            data_val <- assess_out[[metricOI_top]][[metricOI_bot]]
            if (metricOI_bot == "acc" | metricOI_bot == "accuracy"){
              data_val <- data_val / 100
            }
            if (metricOI_bot == "lhs"){
              true_snp_number <- true_nsnp_dict[[format(nsnp$value, scientific = FALSE, trim = TRUE)]][[as.character(ngam$value)]][[as.character(cov$value)]][[as.character(gen_seqerror)]][[as.character(gen_avgrecomb)]][[as.character(rsd)]]
              if (!is.na(true_snp_number)){
                data_val <- data_val / true_snp_number
              }
            }
            if (metricOI_top == "gam_imputation"){
              data_val <- mean(data_val, na.rm = TRUE)
            }
          } else { data_val <- NA}
        }
        data_vec <- c(data_vec, data_val)
        to_add <- mean(data_vec, na.rm = TRUE)
        data_mat[nsnp$index, ngam$index, cov$index, i] <- to_add
      }
    }
  }
}

snp_names <- as_labeller(c(`5000` = "5000 SNPs", `30000` = "30000 SNPs",`1e+05` = "100000 SNPs"))
metric_names <- as_labeller(c("phasing ser" = "Phasing\nSER", "phasing lhs" = "Phasing\nLHS",
                              "gam_imputation ser" = "Inference\nSER", "recomb precision" = "Discovery\nPrecision", 
                              "recomb specificity" = "Discovery\nTNR", "recomb accuracy" = "Discovery\nAccuracy", 
                              "recomb fpr" = "Discovery\nFPR"))
dt <- data.table(c())
for (i in c(1,2)){
  metricOI <- names(last_dim_metric[i])
  for (nsnp in nsnps_iter){
    nsnp_mat <- data.table(data_mat[nsnp$index,,,i]) %>%
      cbind(ngams, .) %>%
      `colnames<-`(c("ngams", covs)) %>%
      pivot_longer(-ngams,  names_to = "coverage") %>%
      as.data.table() %>%
      setnames(., c("n_gametes", "coverage", "value")) %>%
      .[, n_snps := nsnp$value] %>%
      .[, metric := metricOI]
    dt <- rbind(dt, nsnp_mat)
  }
}
g1 <- ggplot(data = dt,
             aes(x = factor(coverage), y = factor(n_gametes), fill = value)) +
  geom_tile() + ggtitle("Donor Haplotype") +
  facet_grid(metric ~ n_snps, labeller = labeller(metric = metric_names, n_snps = snp_names)) +
  theme(panel.background = element_blank(), panel.grid = element_blank()) +
  scale_fill_viridis_c(limits=c(0,1)) + xlab("Coverage (x)") + ylab("Number of gametes") + theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1)) +
  theme(text = element_text(size=14)) + coord_fixed(clip=FALSE)

dt <- data.table(c())
i <- 3
metricOI <- names(last_dim_metric[i])
for (nsnp in nsnps_iter){
  nsnp_mat <- data.table(data_mat[nsnp$index,,,i]) %>%
    cbind(ngams, .) %>%
    `colnames<-`(c("ngams", covs)) %>%
    pivot_longer(-ngams,  names_to = "coverage") %>%
    as.data.table() %>%
    setnames(., c("n_gametes", "coverage", "value")) %>%
    .[, n_snps := nsnp$value] %>%
    .[, metric := metricOI]
  dt <- rbind(dt, nsnp_mat)
}
g2 <- ggplot(data = dt,
             aes(x = factor(coverage), y = factor(n_gametes), fill = value)) +
  geom_tile() + ggtitle("Gamete Genotype") +
  facet_grid(metric ~ n_snps, labeller = labeller(metric = metric_names, n_snps = snp_names)) +
  theme(panel.background = element_blank(), panel.grid = element_blank()) +
  scale_fill_viridis_c(limits=c(0,1)) + xlab("Coverage (x)") + ylab("Number of gametes") + theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1)) +
  theme(text = element_text(size=14)) + coord_fixed(clip=FALSE)

dt <- data.table(c())
for (i in c(4,5,6,7)){
  metricOI <- names(last_dim_metric[i])
  for (nsnp in nsnps_iter){
    nsnp_mat <- data.table(data_mat[nsnp$index,,,i]) %>% 
      cbind(ngams, .) %>% 
      `colnames<-`(c("ngams", covs)) %>% 
      pivot_longer(-ngams,  names_to = "coverage") %>% 
      as.data.table() %>% 
      setnames(., c("n_gametes", "coverage", "value")) %>% 
      .[, n_snps := nsnp$value] %>% 
      .[, metric := metricOI]
    dt <- rbind(dt, nsnp_mat)
  }
}
g3 <- ggplot(data = dt,
             aes(x = factor(coverage), y = factor(n_gametes), fill = value)) +
  geom_tile() + ggtitle("Meiotic Recombination") +
  facet_grid(metric ~ n_snps, labeller = labeller(metric = metric_names, n_snps = snp_names)) +
  theme(panel.background = element_blank(), panel.grid = element_blank()) +
  scale_fill_viridis_c(limits=c(0,1)) + xlab("Coverage (x)") + ylab("Number of gametes") + theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1)) +
  theme(text = element_text(size=14)) + coord_fixed(clip=FALSE)

combo_plot <- grid.arrange(g1, g2, g3, ncol = 1, nrow=3, heights = c(7, 4.4, 12.25))
ggsave("supp_metrics_suppfig2.png", combo_plot, width = 9, height = 18)
