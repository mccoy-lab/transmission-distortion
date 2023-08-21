library(rhapsodi)
library(tidyverse)
library(data.table)
library(pbapply)
library(pbmcapply)

# args1 is iteration, args2 is the outdir
args <- commandArgs(trailingOnly = TRUE)
iteration <- as.numeric(args[1])

set.seed(42)
seeds <- rev(sample(1:1e6, 1000, replace = FALSE))

seed <- seeds[iteration]

random_seed <- seed
threads <- 24

setwd("/scratch/groups/rmccoy22/kweave23/sc_transmission_distortion/runworkflow_20211105/")

donor_chrom_meta <- read.csv("donor_chrom_meta_dt.csv")
#donor, chrom, ngam, and nsnp_segdup columns are important here

prune_meta <- read.csv("ld_prune_snps.csv")
#donor, chrom, wc are important here

donor_chrom_meta <- merge(donor_chrom_meta, prune_meta, by=c("donor", "chrom"))

gamete_meta <- read.csv("gamete_meta_dt.csv")
#donor, chrom, n_crossovers columns are important here

donor_chrom_combo <- paste0(donor_chrom_meta$donor, "_", donor_chrom_meta$chrom)

find_n_crossovers <- function(donor, chrom, gamete_meta){
  roi <- which((gamete_meta$donor == donor) & (gamete_meta$chrom == chrom))
  return(gamete_meta$n_crossovers[roi])
}

find_ngam <- function(donor, chrom, donor_chrom_meta){
  roi_meta <- which(donor_chrom_meta$donor == donor & donor_chrom_meta$chrom == chrom)
  return(donor_chrom_meta$ngam[roi_meta])
}

find_num_snps <- function(donor, chrom, donor_chrom_meta){
  roi_meta <- which(donor_chrom_meta$donor == donor & donor_chrom_meta$chrom == chrom)
  return(donor_chrom_meta$nsnp_segdup[roi_meta])
}

find_num_snps_prune <- function(donor, chrom, donor_chrom_meta){
  roi_meta <- which(donor_chrom_meta$donor == donor & donor_chrom_meta$chrom == chrom)
  return(donor_chrom_meta$wc[roi_meta])
}

find_chrom <- function(donor_chrom_combo_element){
  return(as.integer(unlist(strsplit(donor_chrom_combo_element, "_"))[2]))
}

find_donor <- function(donor_chrom_combo_element){
  return(unlist(strsplit(donor_chrom_combo_element, "_"))[1])
}

gam_full_df_list <- pbmcapply::pbmclapply(1:length(donor_chrom_combo), 
                                          function(x) rhapsodi::sim_generate_global(find_ngam(find_donor(donor_chrom_combo[x]), find_chrom(donor_chrom_combo[x]), donor_chrom_meta),
                                                                                    find_num_snps_prune(find_donor(donor_chrom_combo[x]), find_chrom(donor_chrom_combo[x]), donor_chrom_meta),
                                                                                    find_n_crossovers(find_donor(donor_chrom_combo[x]), find_chrom(donor_chrom_combo[x]), gamete_meta),
                                                                                    random_seed = random_seed),
                                          mc.cores = getOption("mc.cores", threads)) %>% `names<-`(donor_chrom_combo)


calc_dist <- function(dt, chrom, donor) {
  
  # remove positions
  dt <- dt[, -1]
  # calculate distance matrix 
  my_dists <- as.matrix(dist(t(dt))^2)
  my_dists[upper.tri(my_dists, diag = TRUE)] <- NA
  
  dist_long <- my_dists %>%
    as.data.table(keep.rownames = "s1") %>%
    pivot_longer(!s1, names_to = "s2", values_to = "n_mismatch") %>%
    as.data.table() %>%
    .[, tot_snps := nrow(dt)] %>%
    .[, n_match := tot_snps - n_mismatch] %>%
    .[complete.cases(.)] %>%
    .[, donor := donor] %>%
    .[, chrom := chrom]
  
  return(dist_long)
}

calc_dist_wrapper <- function(sim_gametes, index) {
  
  df_name <- names(sim_gametes)[index]
  donor <- strsplit(df_name, "_")[[1]][1]
  chrom <- strsplit(df_name, "_")[[1]][2]
  results <- calc_dist(sim_gametes[[index]], chrom, donor)
  return(results)
  
}

donor_wrapper <- function(sim_gametes, donor_index) {
  donor_ids <- unique(unlist(map(strsplit(names(sim_gametes), "_"), 1)))
  donor_indices <- which(grepl(donor_ids[donor_index], names(sim_gametes)))
  
  combined_results <- rbindlist(pbmclapply(donor_indices, 
                                           function(x) calc_dist_wrapper(sim_gametes, x), 
                                           mc.cores = getOption("mc.cores", threads)))
  
  dist_long_donor <- combined_results %>% 
    group_by(s1, s2) %>%
    summarise(tot_n_match = sum(n_match), tot_snps = sum(tot_snps)) %>%
    as.data.table()
  
  return(dist_long_donor)
}

results_all_donors <- rbindlist(pblapply(1:25, function(x) donor_wrapper(gam_full_df_list, x)))

setwd(args[2])

fwrite(data.table(seed = seed, mean_match = mean(results_all_donors$tot_n_match / results_all_donors$tot_snps)),
       file = paste0("sim_", iteration, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)