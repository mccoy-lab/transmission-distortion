library(itertools)
library(tidyverse)
library(data.table)
library(ggplot2)
#library(gridExtra)
library(patchwork)

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

# for (i in last_dim_metric){
#   metricOI <- names(last_dim_metric[i])
#   metricOI_top <- strsplit(metricOI, " ")[[1]][1]
#   metricOI_bot <- strsplit(metricOI, " ")[[1]][2]
#   for (nsnp in nsnps_iter){
#     for (ngam in ngams_iter){
#       for (cov in covs_iter){
#         for (gse in gse_iter){
#           dir_data_match <- paste0("g", ngam$value, "_s", as.character(as.integer(nsnp$value)), "_c", cov$value, "_se", mod_seqerror, "_r", mod_avgrecomb, "/")
#           dir_data_mismatch <- paste0("g", ngam$value, "_s", as.character(as.integer(nsnp$value)), "_c", cov$value, "_se", gse$value, "_r", mod_avgrecomb, "/")
#           data_vec_match <- c()
#           data_vec_mismatch <- c()
#           for (rsd in rsds){
#             data_file <- paste0("assess_out_rs_", rsd, ".Rdata")
#             fileOI_match <- paste0(dir_base, dir_data_match, data_file)
#             fileOI_mismatch <- paste0(dir_base, dir_data_mismatch, data_file)
#             if (file_test("-f", fileOI_match)){
#               load(fileOI_match)
#               data_val_match <- assess_out[[metricOI_top]][[metricOI_bot]]
#               if (metricOI_bot == "acc" | metricOI_bot == "accuracy"){
#                 data_val_match <- data_val_match / 100
#               }
#               if (metricOI_top == "gam_imputation"){
#                 data_val_match <- mean(data_val_match, na.rm = TRUE)
#               }
#             } else { data_val_match <- NA}
#             if (file_test("-f", fileOI_mismatch)){
#               load(fileOI_mismatch)
#               data_val_mismatch <- assess_out[[metricOI_top]][[metricOI_bot]]
#               if (metricOI_bot == "acc" | metricOI_bot == "accuracy"){
#                 data_val_mismatch <- data_val_mismatch / 100
#               }
#               if (metricOI_top == "gam_imputation"){
#                 data_val_mismatch <- mean(data_val_mismatch, na.rm = TRUE)
#               }
#             } else { data_val_mismatch <- NA}
#           }
#           data_vec_match <- c(data_vec_match, data_val_match)
#           data_vec_mismatch <- c(data_vec_mismatch, data_val_mismatch)
#           to_add <- mean(data_vec_match, na.rm = TRUE) - mean(data_vec_mismatch, na.rm = TRUE)
#           data_mat_gse[nsnp$index, ngam$index, cov$index, gse$index,i] <- to_add
#         }
#         for (gar in gar_iter){
#           dir_data_match <- paste0("g", ngam$value, "_s", as.character(as.integer(nsnp$value)), "_c", cov$value, "_se", mod_seqerror, "_r", mod_avgrecomb, "/")
#           dir_data_mismatch <- paste0("g", ngam$value, "_s", as.character(as.integer(nsnp$value)), "_c", cov$value, "_se", mod_seqerror, "_r", gar$value, "/")
#           data_vec_match <- c()
#           data_vec_mismatch <- c()
#           for (rsd in rsds){
#             data_file <- paste0("assess_out_rs_", rsd, ".Rdata")
#             fileOI_match <- paste0(dir_base, dir_data_match, data_file)
#             fileOI_mismatch <- paste0(dir_base, dir_data_mismatch, data_file)
#             if (file_test("-f", fileOI_match)){
#               load(fileOI_match)
#               data_val_match <- assess_out[[metricOI_top]][[metricOI_bot]]
#               if (metricOI_bot == "acc" | metricOI_bot == "accuracy"){
#                 data_val_match <- data_val_match / 100
#               }
#               if (metricOI_top == "gam_imputation"){
#                 data_val_match <- mean(data_val_match, na.rm = TRUE)
#               }
#             } else { data_val_match <- NA}
#             if (file_test("-f", fileOI_mismatch)){
#               load(fileOI_mismatch)
#               data_val_mismatch <- assess_out[[metricOI_top]][[metricOI_bot]]
#               if (metricOI_bot == "acc" | metricOI_bot == "accuracy"){
#                 data_val_mismatch <- data_val_mismatch / 100
#               }
#               if (metricOI_top == "gam_imputation"){
#                 data_val_mismatch <- mean(data_val_mismatch, na.rm = TRUE)
#               }
#             } else { data_val_mismatch <- NA}
#           }
#           data_vec_match <- c(data_vec_match, data_val_match)
#           data_vec_mismatch <- c(data_vec_mismatch, data_val_mismatch)
#           to_add <- mean(data_vec_match, na.rm = TRUE) - mean(data_vec_mismatch, na.rm = TRUE)
#           data_mat_gar[nsnp$index, ngam$index, cov$index, gar$index,i] <- to_add
#         }
#       }
#     }
#   }
# }

load("data_mat_gar.Rdata")
load("data_mat_gse.Rdata")

make_ggobj <- function(nsnps_iter, plots_list, j, out_base, out_title, to_it_vec, last_dim_metric, data_mat_gse, data_mat_gar){
  snp_names <- as_labeller(c(`5000` = "5000 SNPs", `30000` = "30000 SNPs",`1e+05` = "100000 SNPs"))
  metric_names <- as_labeller(c("phasing acc" = "Phasing\nAccuracy", "phasing com" = "Phasing\nCompleteness",
                                "gam_imputation acc" = "Imputation\nAccuracy", "gam_imputation com" = "Imputation\nCompleteness", "recomb recall" = "Discovery\nTPR", "recomb f1" = "Discovery\nF1 Score",
                                "recomb fdr" = "Discovery\nFDR"))
    se_dt <- data.table(c())
    ar_dt <- data.table(c())
    for (i in to_it_vec){
      metricOI <- names(last_dim_metric[i])
      for (nsnp in nsnps_iter){
        nsnp_mat_se <- data.table(data_mat_gse[nsnp$index,,,j,i]) %>%
          cbind(ngams, .) %>%
          `colnames<-`(c("ngams", covs)) %>%
          pivot_longer(-ngams,  names_to = "coverage") %>%
          as.data.table() %>%
          setnames(., c("n_gametes", "coverage", "difference")) %>%
          .[, n_snps := nsnp$value] %>%
          .[, metric := metricOI]
        se_dt <- rbind(se_dt, nsnp_mat_se)
        nsnp_mat_ar <- data.table(data_mat_gar[nsnp$index,,,j,i]) %>%
          cbind(ngams, .) %>%
          `colnames<-`(c("ngams", covs)) %>%
          pivot_longer(-ngams,  names_to = "coverage") %>%
          as.data.table() %>%
          setnames(., c("n_gametes", "coverage", "difference")) %>%
          .[, n_snps := nsnp$value] %>%
          .[, metric := metricOI]
        ar_dt <- rbind(ar_dt, nsnp_mat_ar)
      }
    }
    if (max(to_it_vec) < 3){
      if (j == 1){
        g_se <- ggplot(data = se_dt,
                         aes(x = factor(coverage), y = factor(n_gametes), fill = difference)) +
          geom_tile() + ggtitle(paste0(toupper(letters[1]), ": ", out_title[j], " Sequencing Error Rate")) +
          facet_grid(metric ~ n_snps, labeller = labeller(metric = metric_names, n_snps = snp_names)) +
          theme(panel.background = element_blank(), panel.grid = element_blank()) +
          scale_fill_viridis_c(limits=c(-1,1)) +
          xlab("Coverage (x)") + ylab("Number of gametes") +
          theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1), text = element_text(size=14)) +
          theme(plot.title = element_text(size = 24, face = "bold")) +
          coord_fixed(clip=FALSE)

        g_ar <- ggplot(data = ar_dt,
                         aes(x = factor(coverage), y = factor(n_gametes), fill = difference)) +
          geom_tile() + ggtitle(paste0(toupper(letters[2]), ": ", out_title[j], " Average Recombination Rate")) +
          facet_grid(metric ~ n_snps, labeller = labeller(metric = metric_names, n_snps = snp_names)) +
          theme(panel.background = element_blank(), panel.grid = element_blank()) +
          scale_fill_viridis_c(limits=c(-1,1)) +
          xlab("Coverage (x)") + ylab("Number of gametes") +
          theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1), text = element_text(size=14)) +
          theme(plot.title = element_text(size = 24, face = "bold")) +
          coord_fixed(clip=FALSE)
        } else {
        g_se <- ggplot(data = se_dt,
                         aes(x = factor(coverage), y = factor(n_gametes), fill = difference)) +
          geom_tile() + ggtitle(paste0(toupper(letters[3]), ": ", out_title[j], " Sequencing Error Rate")) +
          facet_grid(metric ~ n_snps, labeller = labeller(metric = metric_names, n_snps = snp_names)) +
          theme(panel.background = element_blank(), panel.grid = element_blank()) +
          scale_fill_viridis_c(limits=c(-1,1)) +
          xlab("Coverage (x)") + ylab("Number of gametes") +
          theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1), text = element_text(size=14)) +
          theme(plot.title = element_text(size = 24, face = "bold")) +
          coord_fixed(clip=FALSE)

        g_ar <- ggplot(data = ar_dt,
                         aes(x = factor(coverage), y = factor(n_gametes), fill = difference)) +
          geom_tile() + ggtitle(paste0(toupper(letters[4]), ": ", out_title[j], " Average Recombination Rate")) +
          facet_grid(metric ~ n_snps, labeller = labeller(metric = metric_names, n_snps = snp_names)) +
          theme(panel.background = element_blank(), panel.grid = element_blank()) +
          scale_fill_viridis_c(limits=c(-1,1)) +
          xlab("Coverage (x)") + ylab("Number of gametes") +
          theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1), text = element_text(size=14)) +
          theme(plot.title = element_text(size = 24, face = "bold")) +
          coord_fixed(clip=FALSE)
        }
      } else {
          g_se <- ggplot(data = se_dt,
                           aes(x = factor(coverage), y = factor(n_gametes), fill = difference)) +
            geom_tile() +
            facet_grid(metric ~ n_snps, labeller = labeller(metric = metric_names, n_snps = snp_names)) +
            theme(panel.background = element_blank(), panel.grid = element_blank()) +
            scale_fill_viridis_c(limits=c(-1,1)) +
            xlab("Coverage (x)") + ylab("Number of gametes") +
            theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1), text = element_text(size=14)) +
            coord_fixed(clip=FALSE)


          g_ar <- ggplot(data = ar_dt,
                           aes(x = factor(coverage), y = factor(n_gametes), fill = difference)) +
            geom_tile() +
            facet_grid(metric ~ n_snps, labeller = labeller(metric = metric_names, n_snps = snp_names)) +
            theme(panel.background = element_blank(), panel.grid = element_blank()) +
            scale_fill_viridis_c(limits=c(-1,1)) + xlab("Coverage (x)") + ylab("Number of gametes") +
            theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1), text = element_text(size=14)) +
            coord_fixed(clip=FALSE)

        }
        plots_list[[paste0(out_base[j], "_", unlist(strsplit(names(last_dim_metric[i]), " "))[1], "_se")]] <- g_se
        plots_list[[paste0(out_base[j], "_", unlist(strsplit(names(last_dim_metric[i]), " "))[1], "_ar")]] <- g_ar
        return(plots_list)

      }

#Over estimated (generative paramters < model parameters), j = 1
#Under estimated (generative parameters > model parameters), j = 2
list_of_figs <- list()
out_base = c("mparam_gt_gparam", "mparam_lt_gparam")
out_title = c("Overestimated", "Underestimated")

for (j in c(1,2)){
  for (i_vec in list(c(1,2), c(3,4), c(5,6,7))){
    list_of_figs <- make_ggobj(nsnps_iter, list_of_figs, j, out_base, out_title, i_vec, last_dim_metric, data_mat_gse, data_mat_gar)
  }
}


#Want a plot with 12 panels
# phasing heatmaps in the first row
# imputation heatmaps in the second row
# discovery heatmaps in the third row
# First column overestimated seqerror (0.005 rhapsodi and 0.001 gen)
# second column overestimated recomb (1 rhapsodi and 0.6 gen)
# third column underestimated seqerror (0.005 rhapsodi and 0.05 gen)
# fourth column underestimated recomb (1 rhapsodi and 3 gen)
# Try to keep the colorbar in the 4th column and second row, but lose all the rest of them.

combo_plot <- list_of_figs[["mparam_gt_gparam_phasing_se"]] + list_of_figs[["mparam_gt_gparam_phasing_ar"]] + list_of_figs[["mparam_lt_gparam_phasing_se"]] + list_of_figs[["mparam_lt_gparam_phasing_ar"]] +
              list_of_figs[["mparam_gt_gparam_gam_imputation_se"]] + list_of_figs[["mparam_gt_gparam_gam_imputation_ar"]] + list_of_figs[["mparam_lt_gparam_gam_imputation_se"]] + list_of_figs[["mparam_lt_gparam_gam_imputation_ar"]] +
              list_of_figs[["mparam_gt_gparam_recomb_se"]] + list_of_figs[["mparam_gt_gparam_recomb_ar"]] + list_of_figs[["mparam_lt_gparam_recomb_se"]] + list_of_figs[["mparam_lt_gparam_recomb_ar"]] +
              plot_layout(ncol = 4, nrow = 3, widths = c(9, 9, 9, 9), heights = c(7, 7, 9.5), guides="collect")

ggsave("robustness_overall.png", combo_plot, width=42, height=23)
