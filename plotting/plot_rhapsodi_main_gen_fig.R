library(itertools)
library(tidyverse)
library(data.table)
library(ggplot2)
#library(gridExtra)
library(patchwork)

load_data <- TRUE

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

if (!load_data){
  data_mat <- array(rep(NA, num_vals), dim = c(length(nsnps), length(ngams), length(covs), 7))
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
}

snp_names <- as_labeller(c(`5000` = "5000 SNPs", `30000` = "30000 SNPs",`1e+05` = "100000 SNPs"))
#metric_names <- as_labeller(c("phasing acc" = "Phasing\nAccuracy", "phasing com" = "Phasing\nCompleteness",
#                              "gam_imputation acc" = "Imputation\nAccuracy", "gam_imputation com" = "Imputation\nCompleteness", "recomb recall" = "Discovery\nTPR", "recomb f1" = "Discovery\nF1 Score",
#                              "recomb fdr" = "Discovery\nFDR"))
metric_names <- as_labeller(c("phasing acc" = "Accuracy", "phasing com" = "Completeness",
                              "gam_imputation acc" = "Accuracy", "gam_imputation com" = "Completeness", "recomb recall" = "TPR", "recomb f1" = "F1 Score",
                              "recomb fdr" = "FDR"))

if (!load_data){
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
  save(dt, file='main_phasing_dt.Rdata')
} else {
  load('main_phasing_dt.Rdata')
}

dt$outline <- ifelse((dt$n_gametes == 1000) & (dt$coverage == 0.01) & (dt$n_snps == 30000), TRUE, NA)
g1 <- ggplot(data = dt,
     aes(x = factor(coverage), y = factor(n_gametes), fill = value)) +
  geom_tile(show.legend = c(value=TRUE, outline=FALSE)) +
  geom_tile(data = dt[!is.na(dt$outline),] , aes(colour=outline), size = 0.75, show.legend = FALSE) + 
  ggtitle("A: Donor Haplotype Phasing") +
  facet_grid(metric ~ n_snps, labeller = labeller(metric = metric_names, n_snps = snp_names)) +
  theme(panel.background = element_blank(), panel.grid = element_blank()) +
  scale_fill_viridis_c(limits=c(0,1)) +
  xlab("Coverage (x)") + ylab("Number of gametes") + 
  theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1)) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(text = element_text(size=14)) + coord_fixed(clip=FALSE)

if(!load_data){
  dt <- data.table(c())
  for (i in c(3,4)){
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
  save(dt, file='main_imputation_dt.Rdata')
} else {
  load('main_imputation_dt.Rdata')
}

dt$outline <- ifelse((dt$n_gametes == 1000) & (dt$coverage == 0.01) & (dt$n_snps == 30000), TRUE, NA)
g2 <- ggplot(data = dt,
            aes(x = factor(coverage), y = factor(n_gametes), fill = value)) +
  geom_tile(show.legend = c(value=TRUE, outline=FALSE)) +
  geom_tile(data = dt[!is.na(dt$outline),] , aes(colour=outline), size = 0.75, show.legend = FALSE) +
  ggtitle("B: Gamete Genotype Imputation") +
  facet_grid(metric ~ n_snps, labeller = labeller(metric = metric_names, n_snps = snp_names)) +
  theme(panel.background = element_blank(), panel.grid = element_blank()) +
  scale_fill_viridis_c(limits=c(0,1)) + 
  xlab("Coverage (x)") + ylab("Number of gametes") + 
  theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1)) +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(text = element_text(size=14)) + coord_fixed(clip=FALSE)

if (!load_data){
  dt <- data.table(c())
  for (i in c(5,6,7)){
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
  save(dt, file='main_discovery_dt.Rdata')
} else {
  load('main_discovery_dt.Rdata')
}

dt$outline <- ifelse((dt$n_gametes == 1000) & (dt$coverage == 0.01) & (dt$n_snps == 30000), TRUE, NA)
g3 <- ggplot(data = dt,
             aes(x = factor(coverage), y = factor(n_gametes), fill = value)) + 
    geom_tile(show.legend = c(value=TRUE, outline=FALSE)) +
    geom_tile(data = dt[!is.na(dt$outline),] , aes(colour=outline), size = 0.75, show.legend = FALSE) +
    ggtitle("C: Meiotic Recombination Discovery") +
    facet_grid(metric ~ n_snps, labeller = labeller(metric = metric_names, n_snps = snp_names)) +
    theme(panel.background = element_blank(), panel.grid = element_blank()) +
    scale_fill_viridis_c(limits=c(0,1)) + 
    xlab("Coverage (x)") + ylab("Number of gametes") + 
    theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1)) +
    theme(plot.title = element_text(size = 15, face = "bold")) +
    theme(text = element_text(size=14)) + coord_fixed(clip=FALSE)

combo_plot <- g1 + g2 + g3 + plot_layout(ncol = 1, nrow=3, heights = c(7, 7, 10.5), widths=c(9), guides="collect") & theme(legend.position="bottom", legend.text = element_text(angle=45, vjust=0.5), legend.key.width = unit(2.5, "cm"))
ggsave("main_metrics.pdf", combo_plot, width = 9, height = 18)
