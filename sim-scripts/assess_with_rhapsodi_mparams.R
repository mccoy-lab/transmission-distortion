library(rhapsodi)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
num_gametes <- args[1]
num_snps <- args[2]
cov <- args[3]
seq_error <- args[4]
mseq_error <- as.numeric(args[4])
avg_recomb <- args[5]
mavg_recomb <- as.numeric(args[5])
random_seed <- args[6]
threads <- as.integer(args[7])

sim_base <- paste0("/home/kweave23/gamete_data/gen_model_results_noDNM/g",num_gametes,"_s",num_snps,"_c",cov,"_se",seq_error,"_r",avg_recomb,"/runGen_gam_",num_gametes,"_snp_", num_snps, "_cov_", cov, "_seqerr_", seq_error, "_avgr_", avg_recomb, "_rs_", random_seed)
sim_out_dir <- paste0("/home/kweave23/gamete_data/assess_changing_modelparams/g",num_gametes,"_s",num_snps,"_c",cov,"_se",seq_error,"_r",avg_recomb,"_mse", mseq_error, "_mr", mavg_recomb)
if (!dir.exists(sim_out_dir)){
  dir.create(sim_out_dir, showWarnings = FALSE, recursive = TRUE)
}
sim_out <- paste0(sim_out_dir, "/assess_out_rs_", random_seed, ".Rdata")
pred_out <- paste0(sim_out_dir,"/rhapsodi_out_rs_", random_seed, ".Rdata")

survivor_df <- data.frame(rep(0, 7))
survivor_csv_out <- paste0(sim_out_dir, "/survivor_arr_rs_", random_seed, ".csv")

input_file <- paste0(sim_base,"_gametedf_na_truth_afseqednm.csv")
dt <- read.delim(input_file, sep=",", na.strings = c("NA"))
if (nrow(dt) > 0){
  survivor_df[1,1] <- 1
} else {
  write_csv(survivor_df, survivor_csv_out)
}

#Running rhapsodi
#rhapsodi_out <- rhapsodi_autorun(NULL, use_dt = TRUE, input_dt = dt, mcstop=FALSE, threads = threads, seqError_model = mseq_error, avg_recomb_model = mavg_recomb)
standard_input_out <- standard_input(NULL, use_dt = TRUE, input_dt = dt)
complete_haplotypes <- tryCatch(impute_donor_haplotypes(standard_input_out$dt, standard_input_out$positions, threads = threads, mcstop = FALSE),
                                error= function(e) {return(2)})
if (length(complete_haplotypes)==1){
  if (complete_haplotypes == 2){
    write_csv(survivor_df, survivor_csv_out)
    q()
  }
} 
survivor_df[2,1] <- 1

filled_gametes <- tryCatch(fill_gametes(standard_input_out$dt, complete_haplotypes, threads = threads, sequencing_error = mseq_error, avg_recomb = mavg_recomb),
                           error= function(e) {return(2)})
if (length(filled_gametes) == 1){
  if (filled_gametes == 2){
    write_csv(survivor_df, survivor_csv_out)
    q()
  }
}
survivor_df[3,1] <- 1

rhapsodi_out <- tryCatch(report_gametes(TRUE, FALSE, complete_haplotypes, standard_input_out$dt, filled_gametes, standard_input_out$positions, "sampleT", "chrT", threads=threads),
                         error= function(e) {return(2)})
if (length(rhapsodi_out)==1){
  if (rhapsodi_out == 2){
    write_csv(survivor_df, survivor_csv_out)
    q()
  }
} 
survivor_df[4,1] <- 1
save(rhapsodi_out, file = pred_out)

#Assessing rhapsodi
true_ci_file <- paste0(sim_base, "_crossoverIndices_truth_ptfseqednm.csv")
true_ci <- read.delim(true_ci_file, sep = ",", na.strings = c("NA"))

true_donor_haps_file <- paste0(sim_base, "_donorHaps_truth_afseqednm.csv")
true_dh <- read.delim(true_donor_haps_file, sep=",", na.strings = c("NA"))

true_gamete_full_file <- paste0(sim_base, "_gametedf_full_truth_afseqednm.csv")
true_gf <- read.delim(true_gamete_full_file, sep = ",", na.strings = c("NA"))

#assess_out <- sim_assess_it(true_dh, rhapsodi_out$donor_haps, true_ci, rhapsodi_out$recomb_breaks, true_gf, rhapsodi_out$gamete_genotypes)
assess_out <- list()
assess_phasing_out <- tryCatch(sim_assess_phasing(true_dh, rhapsodi_out$donor_haps, nrow(true_dh)),
                               error= function(e) {return(2)})
if (length(assess_phasing_out) == 1){
  if (assess_phasing_out == 2){
    write_csv(survivor_df, survivor_csv_out)
    q()
  }
}
survivor_df[5,1] <- 1
assess_out$phasing <- assess_phasing_out

assess_gam_out <- tryCatch(sim_assess_gam_imputation(true_gf, rhapsodi_out$gamete_genotypes, nrow(true_gf), ncol(true_gf[,-1])),
                           error= function(e) {return(2)})
if (length(assess_gam_out) == 1){
  if (assess_gam_out == 2){
    write_csv(survivor_df, survivor_csv_out)
    save(assess_out, file=sim_out)
    q()
  }
}
survivor_df[6,1] <- 1
assess_out$gam_imputation <- assess_gam_out

assess_recomb_out <- tryCatch(sim_assess_recomb(true_ci, rhapsodi_out$recomb_breaks),
                              error= function(e) {return(2)})
if (length(assess_recomb_out) == 1){
  if (assess_recomb_out == 2){
    write_csv(survivor_df, survivor_csv_out)
    save(assess_out, file=sim_out)
    q()
  }
}
survivor_df[7,1] <- 1
assess_out$recomb <- assess_recomb_out

write_csv(survivor_df, survivor_csv_out)
save(assess_out, file=sim_out)
