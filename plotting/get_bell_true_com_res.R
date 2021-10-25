library(rhapsodi)
library(tidyverse)
library(pbmcapply)

args <- commandArgs(trailingOnly = TRUE)
threads <- as.integer(args[1])
output_file <- args[2]
samples <- c("nc1abnov17", "nc2absept17", "nc3aboct17", "nc4abnov17", "nc6abcd", "nc8ab", "nc9ab", "nc10oldoil", "nc11ab", "nc12ab", "nc13ab", "nc14ab", "nc15ab", "nc16ab", "nc17ab", "nc18ab", "nc22abcd", "nc25abcd", "nc26abcd", "nc27aboct17", "ff3a", "ff4a", "pb2a", "pb3a", "pb4a")
samples_desc <- c(rep("ninf", 20), rep("inf", 5))
chroms <- 1:22
dir_base1 <- "/home/kweave23/gamete_data/sperm_seq_rhapsodi/final_filter_files/"
dir_base2 <- "/csv_out/"

re_recode_gametes <- function(dt, complete_haplotypes) {
  to_return <- data.frame(matrix(NA_real_, nrow=nrow(dt), ncol=ncol(dt)))
  for (i in 1:ncol(dt)) {
    locs_h1 <- dt[,i] == "haplotype1"
    locs_h1[which(is.na(locs_h1))] <- FALSE
    locs_h2 <- dt[,i] == "haplotype2"
    locs_h2[which(is.na(locs_h2))] <- FALSE
    to_return[locs_h1, i] <- complete_haplotypes$h1[locs_h1]
    to_return[locs_h2, i] <- complete_haplotypes$h2[locs_h2]
    colnames(to_return) <- colnames(dt)
  }
  return (to_return)
}

assess_it <- function(dir_base1, sample, sample_desc, dir_base2, chr){
  rhapsodi_haps <- read.delim(paste0(dir_base1, sample, dir_base2, sample, "_", chr, "_parental_hap.csv"), sep=",", na.strings = c("NA"))
  rhapsodi_gametes <- read.delim(paste0(dir_base1, sample, dir_base2, sample, "_", chr, "_filled_sperm_unsmoothed.csv"), sep=",", na.strings = c("NA"))
  rhapsodi_recode_gametes <- re_recode_gametes(rhapsodi_gametes, rhapsodi_haps)
  rhapsodi_recomb <- read.delim(paste0(dir_base1, sample, dir_base2, sample, "_", chr, "_recombination_locs.csv"), sep=",", na.strings = c("NA"))
  full_filter_input <- read.delim(paste0(dir_base1, sample, dir_base2, sample, "_", chr, "_full_filtered_dt.csv"), sep=",", na.strings = c("NA"))
  
  
  num_snps <- nrow(rhapsodi_haps)
  num_gametes <- ncol(rhapsodi_gametes)
  fake_true_phasing <- matrix(rep(0, num_snps*2), ncol = 2, nrow=num_snps) %>% as.data.frame() %>% `colnames<-`(c("donor1", "donor2"))
  fake_true_gametes <- cbind(1:num_snps, matrix(rep(0, num_snps * num_gametes), ncol=num_gametes, nrow = num_snps)) %>% as.data.frame() %>% `colnames<-`(c("positions", colnames(rhapsodi_gametes)))

  phasing_assess <- rhapsodi::sim_assess_phasing(fake_true_phasing, rhapsodi_haps, num_snps)
  phasing_completeness <- phasing_assess$com
  
  imputation_assess <- rhapsodi::sim_assess_gam_imputation(fake_true_gametes, rhapsodi_recode_gametes, num_snps, num_gametes)
  imputation_completeness_vec <- imputation_assess$com
  imputation_completeness_str <- paste0(imputation_completeness_vec, collapse="_")
  
  input_assess <- rhapsodi::sim_assess_gam_imputation(fake_true_gametes, full_filter_input[,-1], num_snps, num_gametes)
  input_completeness_vec <- input_assess$com
  input_completeness_str <- paste0(input_completeness_vec, collapse="_")
  
  res_vec <- rhapsodi_recomb$Genomic_end[!is.na(rhapsodi_recomb$Genomic_end)] - rhapsodi_recomb$Genomic_start[!is.na(rhapsodi_recomb$Genomic_start)]
  res_vec_str <- paste0(res_vec, collapse = "_")
  
  out <- data.frame(sample = sample, sample_desc = sample_desc,  chr = chrom, num_snps = num_snps, num_gametes = num_gametes, phase_com = phasing_completeness, imp_com = imputation_completeness_str, input_com = input_completeness_str, res = res_vec_str)
  return(out)
}

full_df <- data.frame(c())
for (chrom in chroms){
  chrom_df <- do.call(rbind, pbmcapply::pbmclapply(1:length(samples), function(x) assess_it(dir_base1, samples[x], samples_desc[x], dir_base2, chrom),
                                                   mc.cores = getOption("mc.cores", threads)))
  full_df <- rbind(full_df, chrom_df)
}
save(full_df, file=output_file)
