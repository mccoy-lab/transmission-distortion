library(rhapsodi)
library(here)

setwd("/home/kweave23/gamete_data/gen_model_results_noDNM/")
args <- commandArgs(trailingOnly = TRUE)
num_gametes <- args[1]
num_snps <- args[2]
cov <- args[3]
seq_error <- args[4]
avg_recomb <- args[5]
random_seed <- args[6]
threads <- as.integer(args[7])

sim_base <- paste0("g",num_gametes,"_s",num_snps,"_c",cov,"_se",seq_error,"_r",avg_recomb,"/runGen_gam_",num_gametes,"_snp_", num_snps, "_cov_", cov, "_seqerr_", seq_error, "_avgr_", avg_recomb, "_rs_", random_seed)
sim_out <- paste0("g",num_gametes,"_s",num_snps,"_c",cov,"_se",seq_error,"_r",avg_recomb,"/assess_out_rs_", random_seed, ".Rdata")

input_file <- paste0(sim_base,"_gametedf_na_truth_afseqednm.csv")
dt <- read.delim(input_file, sep=",", na.strings = c("NA"))

#Running rhapsodi
rhapsodi_out <- rhapsodi_autorun(NULL, use_dt = TRUE, input_dt = dt, , mcstop=FALSE, threads = threads)

#Assessing rhapsodi
true_ci_file <- paste0(sim_base, "_crossoverIndices_truth_ptfseqednm.csv")
true_ci <- read.delim(true_ci_file, sep = ",", na.strings = c("NA"))

true_donor_haps_file <- paste0(sim_base, "_donorHaps_truth_afseqednm.csv")
true_dh <- read.delim(true_donor_haps_file, sep=",", na.strings = c("NA"))

true_gamete_full_file <- paste0(sim_base, "_gametedf_full_truth_afseqednm.csv")
true_gf <- read.delim(true_gamete_full_file, sep = ",", na.strings = c("NA"))

assess_out <- sim_assess_it(true_dh, rhapsodi_out$donor_haps, true_ci, rhapsodi_out$recomb_breaks, true_gf, rhapsodi_out$gamete_genotypes)
save(assess_out, file=sim_out)
