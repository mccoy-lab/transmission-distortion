#wc of .prune.in files

library(data.table)
library(pbapply)
library(tidyverse)

read_files <- function(file_name){
  fread(file_name) %>% .[, file_name := basename(file_name)] %>% return()
}

flist <- list.files("/scratch/groups/rmccoy22/scarios1/transmission-distortion/run_plink/output_files_plink_110521/", "*.prune.in", full.names = TRUE, recursive = TRUE)
dt_numrows <- rbindlist(pblapply(1:length(flist), function(x) data.frame("wc"=table(read_files(flist[x])$file_name)))) %>% `names<-`(c("file", "wc"))
dt_numrows$donor_chrom <- str_replace(dt_numrows$file, ".prune.in", "")
dt_numrows[, c("donor", "chrom") := tstrsplit(donor_chrom, "_", fixed=TRUE)]
write.csv(dt_numrows, file="/scratch/groups/rmccoy22/kweave23/sc_transmission_distortion/ld_prune_snps.csv")