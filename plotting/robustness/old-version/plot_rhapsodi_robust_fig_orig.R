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
rsds <- c(42, 357, 1848)
gen_seqerror <- c(0.001, 0.05)
gse_iter <- as.list(itertools::enumerate(gen_seqerror))
mod_seqerror <- 0.005
gen_avgrecomb <- c(0.6, 3)
gar_iter <- as.list(itertools::enumerate(gen_avgrecomb))
mod_avgrecomb <- 1
num_vals <- length(nsnps) * length(ngams)  * length(covs) * 2 * 7

data_mat_gse <- array(rep(NA, num_vals), dim=c(length(nsnps), length(ngams), length(covs), length(gen_seqerror), 7))
data_mat_gar <- array(rep(NA, num_vals), dim=c(length(nsnps), length(ngams), length(covs), length(gen_avgrecomb), 7))
last_dim_metric <- list("phasing acc" = 1, 
                        "phasing com" = 2, 
                        "gam_imputation acc" = 3, 
                        "gam_imputation com" = 4, 
                        "recomb recall" = 5, 
                        "recomb fdr" = 6, 
                        "recomb f1" = 7)
dir_base <- "/home/kweave23/gamete_data/gen_model_results_noDNM/"

for (i in last_dim_metric){
  metricOI <- names(last_dim_metric[i])
  metricOI_top <- strsplit(metricOI, " ")[[1]][1]
  metricOI_bot <- strsplit(metricOI, " ")[[1]][2]
  for (nsnp in nsnps_iter){
    for (ngam in ngams_iter){
      for (cov in covs_iter){
        for (gse in gse_iter){
          dir_data_match <- paste0("g", ngam$value, "_s", as.character(as.integer(nsnp$value)), "_c", cov$value, "_se", mod_seqerror, "_r", mod_avgrecomb, "/")
          dir_data_mismatch <- paste0("g", ngam$value, "_s", as.character(as.integer(nsnp$value)), "_c", cov$value, "_se", gse$value, "_r", mod_avgrecomb, "/")
          data_vec_match <- c()
          data_vec_mismatch <- c()
          for (rsd in rsds){
            data_file <- paste0("assess_out_rs_", rsd, ".Rdata")
            fileOI_match <- paste0(dir_base, dir_data_match, data_file)
            fileOI_mismatch <- paste0(dir_base, dir_data_mismatch, data_file)
            if (file_test("-f", fileOI_match)){
              load(fileOI_match)
              data_val_match <- assess_out[[metricOI_top]][[metricOI_bot]]
              if (metricOI_bot == "acc" | metricOI_bot == "accuracy"){
                data_val_match <- data_val_match / 100
              }
              if (metricOI_top == "gam_imputation"){
                data_val_match <- mean(data_val_match, na.rm = TRUE)
              }
            } else { data_val_match <- NA}
            if (file_test("-f", fileOI_mismatch)){
              load(fileOI_mismatch)
              data_val_mismatch <- assess_out[[metricOI_top]][[metricOI_bot]]
              if (metricOI_bot == "acc" | metricOI_bot == "accuracy"){
                data_val_mismatch <- data_val_mismatch / 100
              }
              if (metricOI_top == "gam_imputation"){
                data_val_mismatch <- mean(data_val_mismatch, na.rm = TRUE)
              }
            } else { data_val_mismatch <- NA}
          }
          data_vec_match <- c(data_vec_match, data_val_match)
          data_vec_mismatch <- c(data_vec_mismatch, data_val_mismatch)
          to_add <- mean(data_vec_match, na.rm = TRUE) - mean(data_vec_mismatch, na.rm = TRUE)
          data_mat_gse[nsnp$index, ngam$index, cov$index, gse$index,i] <- to_add
        }
        for (gar in gar_iter){
          dir_data_match <- paste0("g", ngam$value, "_s", as.character(as.integer(nsnp$value)), "_c", cov$value, "_se", mod_seqerror, "_r", mod_avgrecomb, "/")
          dir_data_mismatch <- paste0("g", ngam$value, "_s", as.character(as.integer(nsnp$value)), "_c", cov$value, "_se", mod_seqerror, "_r", gar$value, "/")
          data_vec_match <- c()
          data_vec_mismatch <- c()
          for (rsd in rsds){
            data_file <- paste0("assess_out_rs_", rsd, ".Rdata")
            fileOI_match <- paste0(dir_base, dir_data_match, data_file)
            fileOI_mismatch <- paste0(dir_base, dir_data_mismatch, data_file)
            if (file_test("-f", fileOI_match)){
              load(fileOI_match)
              data_val_match <- assess_out[[metricOI_top]][[metricOI_bot]]
              if (metricOI_bot == "acc" | metricOI_bot == "accuracy"){
                data_val_match <- data_val_match / 100
              }
              if (metricOI_top == "gam_imputation"){
                data_val_match <- mean(data_val_match, na.rm = TRUE)
              }
            } else { data_val_match <- NA}
            if (file_test("-f", fileOI_mismatch)){
              load(fileOI_mismatch)
              data_val_mismatch <- assess_out[[metricOI_top]][[metricOI_bot]]
              if (metricOI_bot == "acc" | metricOI_bot == "accuracy"){
                data_val_mismatch <- data_val_mismatch / 100
              }
              if (metricOI_top == "gam_imputation"){
                data_val_mismatch <- mean(data_val_mismatch, na.rm = TRUE)
              }
            } else { data_val_mismatch <- NA}
          }
          data_vec_match <- c(data_vec_match, data_val_match)
          data_vec_mismatch <- c(data_vec_mismatch, data_val_mismatch)
          to_add <- mean(data_vec_match, na.rm = TRUE) - mean(data_vec_mismatch, na.rm = TRUE)
          data_mat_gar[nsnp$index, ngam$index, cov$index, gar$index,i] <- to_add
        }
      }
    }
  }
}

snp_names <- as_labeller(c(`5000` = "5000 SNPs", `30000` = "30000 SNPs",`1e+05` = "100000 SNPs"))
metric_names <- as_labeller(c("phasing acc" = "Phasing\nAccuracy", "phasing com" = "Phasing\nCompleteness",
                              "gam_imputation acc" = "Inference\nAccuracy", "gam_imputation com" = "Inference\nCompleteness", "recomb recall" = "Discovery\nTPR", "recomb f1" = "Discovery\nF1 Score",  
                              "recomb fdr" = "Discovery\nFDR"))

#Over estimated (generative paramters < model parameters), j = 1
#Under estimated (generative parameters > model parameters), j = 2
out_base = c("mparam_gt_gparam", "mparam_lt_gparam")
out_title = c("Overestimated", "Underestimated")
for (j in c(1,2)){
  dtp_gse <- data.table(c())
  dtp_gar <- data.table(c())
  for (i in c(1,2)){
    metricOI <- names(last_dim_metric[i])
    for (nsnp in nsnps_iter){
      nsnp_mat_gse_p <- data.table(data_mat_gse[nsnp$index,,,j,i]) %>%
        cbind(ngams, .) %>%
        `colnames<-`(c("ngams", covs)) %>%
        pivot_longer(-ngams,  names_to = "coverage") %>%
        as.data.table() %>%
        setnames(., c("n_gametes", "coverage", "difference")) %>%
        .[, n_snps := nsnp$value] %>%
        .[, metric := metricOI]
      dtp_gse <- rbind(dtp_gse, nsnp_mat_gse_p)
      nsnp_mat_gar_p <- data.table(data_mat_gar[nsnp$index,,,j,i]) %>%
        cbind(ngams, .) %>%
        `colnames<-`(c("ngams", covs)) %>%
        pivot_longer(-ngams,  names_to = "coverage") %>%
        as.data.table() %>%
        setnames(., c("n_gametes", "coverage", "difference")) %>%
        .[, n_snps := nsnp$value] %>%
        .[, metric := metricOI]
      dtp_gar <- rbind(dtp_gar, nsnp_mat_gar_p)
    }
  }
  g1_gse <- ggplot(data = dtp_gse,
                   aes(x = factor(coverage), y = factor(n_gametes), fill = difference)) +
    geom_tile() + ggtitle(paste0(out_title[j], " Sequencing Error Rate")) +
    facet_grid(metric ~ n_snps, labeller = labeller(metric = metric_names, n_snps = snp_names)) +
    theme(panel.background = element_blank(), panel.grid = element_blank()) +
    scale_fill_viridis_c(limits=c(-1,1)) + xlab("Coverage (x)") + ylab("Number of gametes") + theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1)) +
    theme(text = element_text(size=14)) + coord_fixed(clip=FALSE)
  
  g1_gar <- ggplot(data=dtp_gar,  
                   aes(x = factor(coverage), y = factor(n_gametes), fill = difference)) +
    geom_tile() + ggtitle(paste0(out_title[j], " Average Recombination Rate")) +
    facet_grid(metric ~ n_snps, labeller = labeller(metric = metric_names, n_snps = snp_names)) +
    theme(panel.background = element_blank(), panel.grid = element_blank()) +
    scale_fill_viridis_c(limits=c(-1,1)) + xlab("Coverage (x)") + ylab("Number of gametes") + theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1)) +
    theme(text = element_text(size=14)) + coord_fixed(clip=FALSE)
  
  combo_plot_ph <- grid.arrange(g1_gse, g1_gar, ncol = 1, nrow=2, top = "Donor Haplotype Phasing Robustness", heights=c(7,7))
  ggsave(paste0("phasing_robustness_", out_base[j], ".png"), combo_plot_ph, width=9, height=15)
  
  
  dti_gse <- data.table(c())
  dti_gar <- data.table(c())
  for (i in c(3,4)){
    metricOI <- names(last_dim_metric[i])
    for (nsnp in nsnps_iter){
      nsnp_mat_gse_i <- data.table(data_mat_gse[nsnp$index,,,j,i]) %>%
        cbind(ngams, .) %>%
        `colnames<-`(c("ngams", covs)) %>%
        pivot_longer(-ngams,  names_to = "coverage") %>%
        as.data.table() %>%
        setnames(., c("n_gametes", "coverage", "difference")) %>%
        .[, n_snps := nsnp$value] %>%
        .[, metric := metricOI]
      dti_gse <- rbind(dti_gse, nsnp_mat_gse_i)
      nsnp_mat_gar_i <- data.table(data_mat_gar[nsnp$index,,,j,i]) %>%
        cbind(ngams, .) %>%
        `colnames<-`(c("ngams", covs)) %>%
        pivot_longer(-ngams,  names_to = "coverage") %>%
        as.data.table() %>%
        setnames(., c("n_gametes", "coverage", "difference")) %>%
        .[, n_snps := nsnp$value] %>%
        .[, metric := metricOI]
      dti_gar <- rbind(dti_gar, nsnp_mat_gar_i)
    }
  }
  g2_gse <- ggplot(data = dti_gse,
                   aes(x = factor(coverage), y = factor(n_gametes), fill = difference)) +
    geom_tile() + ggtitle(paste0(out_title[j], " Sequencing Error Rate")) +
    facet_grid(metric ~ n_snps, labeller = labeller(metric = metric_names, n_snps = snp_names)) +
    theme(panel.background = element_blank(), panel.grid = element_blank()) +
    scale_fill_viridis_c(limits=c(-1,1)) + xlab("Coverage (x)") + ylab("Number of gametes") + theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1)) +
    theme(text = element_text(size=14)) + coord_fixed(clip=FALSE)
  
  g2_gar <- ggplot(data = dti_gar,
                   aes(x = factor(coverage), y = factor(n_gametes), fill = difference)) +
    geom_tile() + ggtitle(paste0(out_title[j], " Average Recombination Rate")) +
    facet_grid(metric ~ n_snps, labeller = labeller(metric = metric_names, n_snps = snp_names)) +
    theme(panel.background = element_blank(), panel.grid = element_blank()) +
    scale_fill_viridis_c(limits=c(-1,1)) + xlab("Coverage (x)") + ylab("Number of gametes") + theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1)) +
    theme(text = element_text(size=14)) + coord_fixed(clip=FALSE)
  
  combo_plot_imp <- grid.arrange(g2_gse, g2_gar, ncol = 1, nrow=2, top = "Gamete Genotype Inference Robustness", heights = c(7, 7))
  ggsave(paste0("imputation_robustness_", out_base[j], ".png"), combo_plot_imp, width = 9, height = 15)
  
  dtd_gse <- data.table(c())
  dtd_gar <- data.table(c())
  for (i in c(5,6,7)){
    metricOI <- names(last_dim_metric[i])
    for (nsnp in nsnps_iter){
      nsnp_mat_gse_d <- data.table(data_mat_gse[nsnp$index,,,j,i]) %>%
        cbind(ngams, .) %>%
        `colnames<-`(c("ngams", covs)) %>%
        pivot_longer(-ngams,  names_to = "coverage") %>%
        as.data.table() %>%
        setnames(., c("n_gametes", "coverage", "difference")) %>%
        .[, n_snps := nsnp$value] %>%
        .[, metric := metricOI]
      dtd_gse <- rbind(dtd_gse, nsnp_mat_gse_d)
      nsnp_mat_gar_d <- data.table(data_mat_gar[nsnp$index,,,j,i]) %>%
        cbind(ngams, .) %>%
        `colnames<-`(c("ngams", covs)) %>%
        pivot_longer(-ngams,  names_to = "coverage") %>%
        as.data.table() %>%
        setnames(., c("n_gametes", "coverage", "difference")) %>%
        .[, n_snps := nsnp$value] %>%
        .[, metric := metricOI]
      dtd_gar <- rbind(dtd_gar, nsnp_mat_gar_d)
    }
  }
  g3_gse <- ggplot(data = dtd_gse,
                   aes(x = factor(coverage), y = factor(n_gametes), fill = difference)) +
    geom_tile() + ggtitle(paste0(out_title[j] , " Sequencing Error Rate")) +
    facet_grid(metric ~ n_snps, labeller = labeller(metric = metric_names, n_snps = snp_names)) +
    theme(panel.background = element_blank(), panel.grid = element_blank()) +
    scale_fill_viridis_c(limits=c(-1,1)) + xlab("Coverage (x)") + ylab("Number of gametes") + theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1)) +
    theme(text = element_text(size=14)) + coord_fixed(clip=FALSE)
  
  g3_gar <- ggplot(data = dtd_gar,
                   aes(x = factor(coverage), y = factor(n_gametes), fill = difference)) +
    geom_tile() + ggtitle(paste0(out_title[j] ," Average Recombination Rate")) +
    facet_grid(metric ~ n_snps, labeller = labeller(metric = metric_names, n_snps = snp_names)) +
    theme(panel.background = element_blank(), panel.grid = element_blank()) +
    scale_fill_viridis_c(limits=c(-1,1)) + xlab("Coverage (x)") + ylab("Number of gametes") + theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1)) +
    theme(text = element_text(size=14)) + coord_fixed(clip=FALSE)
  
  combo_plot_re <- grid.arrange(g3_gse, g3_gar, ncol = 1, nrow=2, top = "Meiotic Recombination Discovery Robustness", heights=c(9.5, 9.5))
  ggsave(paste0("recombination_robustness_", out_base[j], ".png"), combo_plot_re, width=9, height=18)
}  

