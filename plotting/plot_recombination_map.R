library(readr)
library(tidyverse)
library(data.table)
library(pbapply)
library(ggplot2)
library(ggpubr)

load_data <- TRUE

if (!load_data){
window_size <- 1e6

donor_chrom_meta <- read.csv("/home/kweave23/scr16_rmccoy22/kweave23/sc_transmission_distortion/runworkflow_20211105/donor_chrom_meta_dt.csv", row.names =1)

donors <- unique(donor_chrom_meta$donor)
chroms <- sort(unique(donor_chrom_meta$chrom))

chrom_sizes <- read_delim("/home/kweave23/scr4_rmccoy22/kweave23/resources/hg38.chrom.sizes_main.txt", delim="\t", col_names = c("chrom", "size"))
chrom_sizes$num_windows <- ceiling(chrom_sizes$size/window_size)

loc_to_window <- function(loc, window_size){
  window <- ceiling(loc/window_size)
  return(window)
}

rown <- max(chrom_sizes$num_windows)
coln <- length(chroms)
donorn <- length(donors)
storage_mat_donor <- array(rep(0, rown*coln*donorn), c(rown, coln, donorn)) %>% `dimnames<-`(list(c(), paste0("chr", 1:coln), donors))

read_files <- function(file_name){
  fread(file_name) %>% .[, file_name := basename(file_name)] %>% return()
}

flist_recomb <- list.files("/home/kweave23/scr16_rmccoy22/kweave23/sc_transmission_distortion/runworkflow_20211105", "_recombination_locs.csv", full.names=TRUE, recursive = TRUE)

dt_recomb <- rbindlist(pblapply(1:length(flist_recomb), function(x) read_files(flist_recomb[x])))

dt_recomb$midpoint <- (dt_recomb$Genomic_start + dt_recomb$Genomic_end)/2
dt_recomb[, c("donor", "chrom", "gamete") := tstrsplit(Ident, "_", fixed=TRUE)]
dt_recomb <- dt_recomb[!is.na(Genomic_start)]

for (chrom in chroms){
  for (donor in donors){
    roi <- which(dt_recomb$chrom == chrom & dt_recomb$donor == donor)
    window_mid <- loc_to_window(dt_recomb$midpoint[roi], window_size)
    table_wm <- table(window_mid)
    storage_mat_donor[as.integer(names(table_wm)), paste0("chr", chrom), donor] <- table_wm
  }
}

for (i in 1:length(chrom_sizes$chrom)){
  chrom <- chrom_sizes$chrom[i]
  max_chrom_window <- chrom_sizes$num_windows[i]
  if ((max_chrom_window < rown) & (chrom %in% paste0("chr",chroms))){
    storage_mat_donor[(max_chrom_window+1):rown, chrom,] <- NA
  }
}

melted_df_donor <- pivot_longer(cbind(bin = 1:rown, as.data.frame(storage_mat_donor)), cols = starts_with("chr"), values_to = "count", names_to = "meta") %>%
  separate(meta, c("chr", "donor"), sep ="\\.", remove=TRUE)

melted_df <- melted_df_donor %>% group_by(bin, chr) %>% summarize(count = sum(count))

edit_df <- function(dfoi, window_size = 1e6, chroms = paste0("chr", 1:22), merge_bool = FALSE, to_merge = NULL){
  if (!merge_bool){
    dfoi$gen_pos <- dfoi$bin * window_size
    dfoi$chr_num <- as.integer(str_replace(dfoi$chr, "chr", ""))
    dfoi <- setorder(dfoi, chr_num)
    dfoi$chr_num <- factor(dfoi$chr_num, labels=c(paste0("chr", 1:length(chroms))))
    dfoi <- drop_na(dfoi)
    return(dfoi)
  } else{
    dfoi_merged <- merge(dfoi, to_merge, by=c("chr", "bin"), all = TRUE) %>%
      replace_na(list(rate=0))
    return(dfoi_merged)
  }
}

melted_df <- edit_df(melted_df, window_size =window_size, chroms = chroms)
melted_df_donor <- edit_df(melted_df_donor, window_size =window_size, chroms = chroms)

save(melted_df, file="supp_recomb/recomb_map_meltdf.Rdata")
#save(melted_df_donor, file = "supp_recomb/recomb_map_meltdonordf.Rdata")

paternal_map <- read.delim("/home/kweave23/scr4_rmccoy22/kweave23/resources/aau1043_datas1", skip=7)
paternal_map$interval_size <- paternal_map$End - paternal_map$Begin
paternal_map$midpoint <- (paternal_map$Begin + paternal_map$End) / 2
paternal_map$bin <- loc_to_window(paternal_map$midpoint, window_size)

paternal_map_wa <- paternal_map %>%
  group_by(Chr, bin) %>%
  summarize(rate = weighted.mean(cMperMb, interval_size/window_size)) %>%
  `colnames<-`(c("chr", "bin", "rate"))

merged_df <- edit_df(melted_df, merge_bool = TRUE, to_merge = paternal_map_wa)
save(merged_df, file="supp_recomb/recomb_map_decodecomp_merged_df.Rdata")
merged_df_donor <- edit_df(melted_df_donor, merge_bool = TRUE, to_merge = paternal_map_wa)
save(merged_df_donor, file="supp_recomb/recomb_map_decodecomp_merged_df.Rdata")

} else {
load("supp_recomb/recomb_map_meltdf.Rdata")
load("supp_recomb/recomb_map_decodecomp_merged_df.Rdata")
load("supp_recomb/recomb_map_decodecomp_merged_df.Rdata")
}

p <- ggplot(melted_df, aes(x=gen_pos, y=count)) + geom_bar(stat='identity') +
     facet_wrap(~chr_num, nrow=floor(22/4), scales="free_x", shrink = FALSE, labeller=label_parsed) +
     theme_bw() + theme(panel.background = element_blank(), panel.grid = element_blank()) +
     xlab("Genomic position") + ylab("Inferred number of crossovers") +
     theme(axis.text.x = element_text(angle=50, vjust=1, hjust = 1)) + theme(text = element_text(size=10))

ggsave("plots/recombination_map.pdf", plot=p, height=7.48, width=6.78, units="in")

p2 <- ggplot(merged_df, aes(x=rate, y = count)) + geom_point(alpha = 0.5) +
  theme_bw() + theme(panel.background = element_blank(), panel.grid = element_blank()) +
  xlab("Published male-specific recombination rate (cM/Mb)") + ylab("Inferred number of crossovers per genomic bin from Sperm-seq") +
  stat_cor()

ggsave("plots/decode_comparison.pdf", plot=p2)

p3 <- ggplot(merged_df_donor, aes(x=rate, y = count)) + geom_point(alpha=0.5) +
  theme_bw() + theme(panel.background = element_blank(), panel.grid = element_blank()) +
  xlab("Published male-specific recombination rate (cM/Mb)") + ylab("Inferred number of crossovers per genomic bin from Sperm-seq") +
  stat_cor() +
  facet_wrap(~donor, nrow=5, scales="free_x", shrink = FALSE)

ggsave("plots/decode_comparison_facetdonor.pdf", plot=p3)

bins_to_color <- read.csv("/home/kweave23/scr16_rmccoy22/kweave23/sc_transmission_distortion/runworkflow_20211105/filtered_data_specific_first_last_bins.csv", row.names =1) %>% `colnames<-`(c("bin", "freq", "color", "chr"))
merged_df <- merged_df %>% mutate(across(chr, gsub, pattern="chr", replacement = ""))
merged_df <- merge(merged_df, bins_to_color, by =c("chr", "bin"), all = TRUE) %>%
  replace_na(list(color="black", freq=0))

p4 <- ggplot(merged_df, aes(x=rate, y=count, color=color)) + geom_point(alpha=0.5) +
  theme_bw() + theme(panel.background = element_blank(), panel.grid = element_blank()) +
  xlab("Published male-specific recombination rate (cM/Mb)") + ylab("Inferred number of crossovers per genomic bin from Sperm-seq") +
  scale_color_manual(name = 'chromosome  bin', values = c("black", "purple", "blue"), labels = c('not first or last bin', 'first bin', 'last bin')) +
  stat_cor()

ggsave("plots/decode_comparison_binloccolor.pdf", plot=p4)

p5 <- ggplot(merged_df, aes(x=rate, y=count, color=color, size=freq)) + geom_point(alpha=0.5) +
  theme_bw() + theme(panel.background = element_blank(), panel.grid = element_blank()) +
  xlab("Published male-specific recombination rate (cM/Mb)") + ylab("Inferred number of crossovers per genomic bin from Sperm-seq") +
  scale_color_manual(name = 'chromosome  bin', values = c("black", "purple", "blue"), labels = c('not first or last bin', 'first bin', 'last bin'))

ggsave("plots/decode_comparison_binloccolorfreqsize.pdf", plot=p5)

correlations_donorspeconly <- unlist(lapply(1:length(donors), function(x) cor.test(merged_df_donor[which(merged_df_donor$donor == donors[x]), "rate"], merged_df_donor[which(merged_df_donor$donor == donors[x]), "count"])$estimate))
cor_dspec_df <- data.frame(donor = donors, pearson_cor = correlations_donorspeconly)

foranovadonorspec <- full_join(donor_meta_df %>% group_by(donor, fertility) %>% summarize(avgcov = mean(cov)) %>% select(donor, fertility, avgcov), cor_dspec_df, by = "donor") %>%
  mutate(pearson_cor = as.numeric(pearson_cor)) %>%
  mutate(fertility = as.factor(fertility)) %>%
  mutate(avgcov = as.numeric(avgcov))

resdspec <- lm(pearson_cor ~ fertility + avgcov, data = foranovadonorspec)
summary(resdspec)
