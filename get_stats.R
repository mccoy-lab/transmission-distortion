library(tidyverse)
library(data.table)
library(pbapply)

read_files <- function(file_name){
  fread(file_name) %>% .[, file_name := basename(file_name)] %>% return()
}


#SAMPLE SPECIFIC
flist_sample_meta <- list.files("/scratch/groups/rmccoy22/kweave23/sc_transmission_distortion/runworkflow_20211105/",
                                "_sample_meta.csv", 
                                full.names = TRUE, recursive = TRUE)

dt_sample_meta <- rbindlist(pblapply(1:length(flist_sample_meta), function(x) read_files(flist_sample_meta[x])))

##number of SNPs
total_snps <- sum(dt_sample_meta$nsnp_in)
total_snps_af <- sum(dt_sample_meta$nsnp_segdup)
snps_filt <- total_snps - total_snps_af #15138461
snps_filt_per <- snps_filt / total_snps * 100 #30.315

##number of gametes
donor_spec_gam <- lapply(1:length(unique(dt_sample_meta$donor)), function(x) dt_sample_meta$ngam[which(unique(dt_sample_meta$donor)[x] == dt_sample_meta$donor)]) %>% `names<-`(unique(dt_sample_meta$donor))
max_donor <- unlist(lapply(1:length(donor_spec_gam), function(x) max(donor_spec_gam[[x]])))
total_gam <- sum(max_donor) #41189
min_gam <- min(dt_sample_meta$ngam) #969
max_gam <- max(dt_sample_meta$nngam) #3377

##phasing COM
mean_phase_com <- mean(dt_sample_meta$phasing_com) * 100 #99.90375
sd_phase_com <- sd(dt_sample_meta$phasing_com) * 100 #0.1573844


#GAMETE SPECIFIC
flist_gamete_meta <- list.files("/scratch/groups/rmccoy22/kweave23/sc_transmission_distortion/runworkflow_20211105/",
                                "_gamete_meta.csv", 
                                full.names = TRUE, recursive = TRUE)

dt_gamete_meta <- rbindlist(pblapply(1:length(flist_gamete_meta), function(x) read_files(flist_gamete_meta[x])))

##input COM and cov
mean_input_com <- mean(dt_gamete_meta$input_com) * 100 #1.146296
sd_input_com <- sd(dt_gamete_meta$input_com) * 100 #0.76332
mean_input_cov <- mean(dt_gamete_meta$input_approx_cov) #0.01155937
sd_input_cov <- sd(dt_gamete_meta$input_approx_cov) #0.007797991

##imputation COM
mean_imp_com <- mean(dt_gamete_meta$output_com) * 100 #99.29736
sd_imp_com <- sd(dt_gamete_meta$output_com) * 100 #1.231243

##recomb events
mean_recomb_per <- mean(dt_gamete_meta$n_crossovers) #1.171595
median_recomb_per <- median(dt_gamete_meta$n_crossovers) #1

recomb_summary_by_gamete <- group_by(dt_gamete_meta, donor, gamete) %>%
  summarize(., tot_crossovers = sum(n_crossovers)) %>%
  as.data.table() %>%
  setorder(tot_crossovers)

recomb_summary_across_gametes <- table(recomb_summary_by_gamete$tot_crossovers) %>% 
  as.data.table() %>%
  setnames(., c("genome_wide_crossovers", "n_gametes")) %>%
  setorder(., -n_gametes)
recomb_summary_across_gametes[, genome_wide_crossovers := as.numeric(genome_wide_crossovers)]

mode_genome_wide <- recomb_summary_across_gametes$genome_wide_crossovers[1] #24

flist_recomb_meta <- list.files("/scratch/groups/rmccoy22/kweave23/sc_transmission_distortion/runworkflow_20211105/",
                                "_recomb_res.csv", 
                                full.names = TRUE, recursive = TRUE)

dt_recomb_meta <- rbindlist(pblapply(1:length(flist_recomb_meta), function(x) read_files(flist_recomb_meta[x])))

##breakpoint resolution
mean_res <- mean(dt_recomb_meta$res) #663582.7
sd_res <- sd(dt_recomb_meta$res) #1254802
median_res <- median(dt_recomb_meta$res) #357381