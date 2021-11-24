library(readr)
library(tidyverse)
library(data.table)
library(pbapply)
library(ggplot2)

window_size <- 1e6

donor_chrom_meta <- read.csv("/scratch/groups/rmccoy22/kweave23/sc_transmission_distortion/runworkflow_20211105/donor_chrom_meta_dt.csv")

donors <- unique(donor_chrom_meta$donor)
chroms <- sort(unique(donor_chrom_meta$chrom))

chrom_sizes <- read_delim("~/resources/hg38.chrom.sizes_main.txt", delim="\t", col_names = c("chrom", "size"))
chrom_sizes$num_windows <- ceiling(chrom_sizes$size/1e6)

loc_to_window <- function(loc, window_size){
  window <- ceiling(loc/window_size)
  return(window)
}

rown <- max(chrom_sizes$num_windows)
coln <- length(chroms)
storage_mat <- matrix(rep(0, rown*coln), nrow = rown , ncol=coln) %>% `colnames<-`(paste0("chr", 1:coln))

read_files <- function(file_name){
  fread(file_name) %>% .[, file_name := basename(file_name)] %>% return()
}

flist_recomb <- list.files("/scratch/groups/rmccoy22/kweave23/sc_transmission_distortion/runworkflow_20211105/", "_recombination_locs.csv", full.names=TRUE, recursive = TRUE)

dt_recomb <- rbindlist(pblapply(1:length(flist_recomb), function(x) read_files(flist_recomb[x])))

dt_recomb$midpoint <- (dt_recomb$Genomic_start + dt_recomb$Genomic_end)/2
dt_recomb[, c("donor", "chrom", "gamete") := tstrsplit(Ident, "_", fixed=TRUE)]
dt_recomb <- dt_recomb[!is.na(Genomic_start)]

for (chrom in chroms){
  roi <- which(dt_recomb$chrom == chrom)
  window_mid <- loc_to_window(dt_recomb$midpoint[roi], window_size)
  table_wm <- table(window_mid)
  storage_mat[as.integer(names(table_wm)), paste0("chr", chrom)] <- table_wm
}

for (i in 1:length(chrom_sizes$chrom)){
  chrom <- chrom_sizes$chrom[i]
  max_chrom_window <- chrom_sizes$num_windows[i]
  if ((max_chrom_window < rown) & (chrom %in% paste0("chr",chroms))){
    storage_mat[(max_chrom_window+1):rown, chrom] <- NA
  }
}

melted_df <- pivot_longer(data.frame(cbind(bin=1:nrow(storage_mat), storage_mat)), starts_with("chr"), names_to = "chr", values_to = "count")
melted_df$gen_pos <- melted_df$bin * window_size
melted_df$chr_num <- as.integer(str_replace(melted_df$chr, "chr", ""))
melted_df <- setorder(melted_df, chr_num)
melted_df$chr_num <- factor(melted_df$chr_num, labels=c(paste0("chr", 1:length(chroms))))
melted_df <- drop_na(melted_df)
p <- ggplot(melted_df, aes(x=gen_pos, y=count)) + geom_bar(stat='identity') + 
     facet_wrap(~chr_num, nrow=floor(length(chroms)/4), scales="free_x", shrink = FALSE, labeller=label_parsed) + 
     theme_bw() + theme(panel.background = element_blank(), panel.grid = element_blank()) + 
     xlab("Genomic position") + ylab("Inferred number of crossovers") + 
     theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1)) + theme(text = element_text(size=10))

ggsave("/scratch/groups/rmccoy22/kweave23/sc_transmission_distortion/plot_rdata/plots/recombination_map.pdf", plot=p, height=7.48, width=6.78, units="in")