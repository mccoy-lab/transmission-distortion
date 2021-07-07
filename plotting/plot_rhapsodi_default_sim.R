library(itertools)
library(tidyverse)
library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

nsnps <- c(5000, 30000, 100000)
nsnps_iter <- as.list(itertools::enumerate(nsnps))
ngams <- c(3, 15, 50, 150, 500, 1000, 2500, 5000)
ngams_iter <- as.list(itertools::enumerate(ngams))
covs <- c(0.001, 0.01, 0.1, 0.223, 0.357, 0.511, 0.693, 1.204, 2.303)
covs_iter <- as.list(itertools::enumerate(covs))
num_vals <- length(nsnps) * length(ngams)  * length(covs)
rsds <- c(42, 357, 1848)
gen_seqerror <- args[1]
mod_seqerror <- 0.005
gen_avgrecomb <- args[2]
mod_avgrecomb <- 1
top_three <- args[3]
metricOI <- args[4]
fig_title <- args[5]

stopifnot(top_three == "phasing" | top_three == "gam_imputation" | top_three == "recomb")
if (top_three == "phasing" | top_three == "gam_imputation" ){
  stopifnot(metricOI == "acc" | metricOI == "com" | metricOI == "lhs" | metricOI == "ser")
} else if (top_three == "recomb"){
  stopifnot(metricOI == "precision" | metricOI == "recall" | metricOI == "accuracy" | metricOI == "specificity" | metricOI == "fdr" | metricOI == "fpr" | metricOI == "f1" | metricOI == "true_n" | metricOI == "pred_n" | metricOI == "tn" | metricOI == "fn" | metricOI == "tp" | metricOI == "fp")
}

data_mat <- array(rep(NA, num_vals), dim = c(length(nsnps), length(ngams), length(covs)))

dir_base <- "/home/kweave23/gamete_data/gen_model_results_noDNM/"
for (nsnp in nsnps_iter){
  for (ngam in ngams_iter){
    for (cov in covs_iter){
      dir_data <- paste0("g", ngam$value, "_s", as.character(as.integer(nsnp$value)), "_c", cov$value, "_se", gen_seqerror, "_r", gen_avgrecomb, "/")
      data_vec <- c()
      outfig  <- paste0(dir_base, "heatmap_", top_three, "_", metricOI, "_", gen_seqerror, "_", gen_avgrecomb, ".png")
      for (rsd in rsds){
        data_file <- paste0("assess_out_rs_", rsd, ".Rdata")
        fileOI <- paste0(dir_base, dir_data, data_file)
        if (file_test("-f", fileOI)){
          load(fileOI)
          data_val <- assess_out[[top_three]][[metricOI]]
          if (top_three == "gam_imputation"){
            data_val <- mean(data_val, na.rm = TRUE)
          }
          
        } else { data_val <- NA}
      }
      data_vec <- c(data_vec, data_val)
      to_add <- mean(data_vec, na.rm = TRUE)
      data_mat[nsnp$index, ngam$index, cov$index] <- to_add
    }
  }
}

dt <- data.table(c())
for (nsnp in nsnps_iter){
  nsnp_mat <- data.table(data_mat[nsnp$index,,]) %>% 
    cbind(ngams, .) %>% 
    `colnames<-`(c("ngams", covs)) %>% 
    pivot_longer(-ngams,  names_to = "coverage") %>% 
    as.data.table() %>% 
    setnames(., c("n_gametes", "coverage", "value")) %>% 
    .[, n_snps := nsnp$value] %>% 
    .[, metric := metricOI]
  dt <- rbind(dt, nsnp_mat)
}

snp_names <- as_labeller(c(`5000` = "5000 SNPs", `30000` = "30000 SNPs",`1e+05` = "100000 SNPs"))
metric_names <- as_labeller(c("acc" = "Accuracy (%)", "com" = "Completeness (%)", "ser" = "Switch Error Rate", "lhs" = "Longest Haplotype Segment",
                              "accuracy" = "Accuracy (%)", "precision" = "Precision", "recall" = "Recall", "f1" = "F1 Score",  "specificity" = "Specificity",
                              "fdr" = "False Discovery Rate", "fpr" = "False Positive Rate", "true_n" = "# of true breakpoints", "pred_n" = "# of predicted breakpoints", 
                              "tp" = "True Positives", "fp" = "False Positives", "tn" =  "True Negatives", "fn" = "False Negatives"))
g <- ggplot(data = dt,
       aes(x = factor(coverage), y = factor(n_gametes), fill = value)) +
  geom_tile() + ggtitle(fig_title) +
  facet_grid(metric ~ n_snps, labeller = labeller(metric = metric_names, n_snps = snp_names)) +
  theme(panel.background = element_blank(), panel.grid = element_blank()) +
  scale_fill_viridis_c() + xlab("Coverage (x)") + ylab("Number of gametes") + theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1)) +
  theme(text = element_text(size=14)) + coord_fixed(clip=FALSE)
ggsave(outfig)
