library(rhapsodi)
library(tidyverse)
library(pbmcapply)

random_seed <- **
threads <- **

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




