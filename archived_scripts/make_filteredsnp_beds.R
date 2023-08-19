library(tidyverse)
library(rhapsodi)
library(data.table)

filter_files <- list.files("/home/kweave23/data_rmccoy22/sperm_seq_rhapsodi/bell_data/final_filter_files", "full_filtered_dt.csv", full.names = TRUE, recursive = TRUE)
for (filename in filter_files){
  input_dt <- read.delim(filename, sep=",", na.strings = c("NA"))
  input_data <- rhapsodi::read_data(NULL, use_dt = TRUE, input_dt = input_dt)
  nonnapos <- rbindlist(lapply(1:ncol(input_data$dt), function(x) data.frame(input_data$positions[which(!is.na(input_data$dt[,x]))])))
  donor <- unlist(strsplit(unlist(strsplit(filename, '/'))[9], "_"))[1]
  chr <- paste0("chr", unlist(strsplit(unlist(strsplit(filename, '/'))[9], "_"))[2])
  beddt <- data.frame(chrom = chr, pos1 = nonnapos - 1, pos2 = nonnapos + 1) %>% `colnames<-`(c("chrom", "chromStart", "chromEnd"))
  write.table(beddt, paste0("/home/kweave23/scr16_rmccoy22/kweave23/sc_transmission_distortion/runworkflow_20211105/", donor, "/csv_out/assess_filtered/", donor, "_", chr, "_nonnapos.bed"), row.names = FALSE, quote = FALSE, col.names = FALSE, sep="\t")
}