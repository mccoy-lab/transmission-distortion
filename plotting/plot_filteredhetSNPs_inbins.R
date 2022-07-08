library(tidyverse)
library(ggplot2)
library(data.table)

donors <- c("ff3a", "ff4a", "nc10oldoil", "nc11ab", "nc12ab", "nc13ab", "nc14ab", "nc15ab", "nc16ab", "nc17ab", "nc18ab", "nc1abnov17", "nc22abcd",
            "nc25abcd", "nc26abcd", "nc27aboct17", "nc2absept17", "nc3aboct17", "nc4abnov17", "nc6abcd", "nc8ab", "nc9ab", "pb2a", "pb3a", "pb4a")

binnedfilelist <- list.files("/home/kweave23/scr16_rmccoy22/kweave23/sc_transmission_distortion/runworkflow_20211105", "0pt5Mbp_coverage.txt", full.names = TRUE, recursive = TRUE)

read_files <- function(filename){
  dt <- read_delim(filename, delim = "\t", col_names = FALSE) %>% `colnames<-`(c("Chrom", "WindowStart", "WindowEnd", "Count"))
  donor <- unlist(strsplit(filename, '/'))[8]
  dt$Donor <- donor
  dt$Window <- 1:nrow(dt)
  return(dt)
}

full_dt <- rbindlist(lapply(1:length(binnedfilelist), function(x) read_files(binnedfilelist[x])))

for (chrom in paste0("chr", 1:22)){
  p <- ggplot(full_dt[full_dt$Chrom == chrom,], aes(x=Window, y=Count)) + geom_bar(stat = "identity") + facet_wrap(~Donor) + 
    theme_bw() + theme(panel.background = element_blank(), panel.grid = element_blank()) +
    xlab("Bin (0.5 Mbps)") + ylab("filtered hetSNP Count") + ggtitle(chrom)
  ggsave(paste0("/home/kweave23/scr16_rmccoy22/kweave23/sc_transmission_distortion/plot_rdata/plots/filteredhetSNP_0pt5Mbp_coverage", chrom, ".pdf"), plot = p, height=8.31, width=8.31, units="in")
}

full_dt <- full_dt %>% unite(chromwindow, c("Chrom", "Window"), remove = FALSE, sep="_")
roizerocount <- which(full_dt$Count == 0)
full_dt_subset <- full_dt[roizerocount]
full_dt_zerocountsummarize <- full_dt_subset %>% group_by(chromwindow) %>%
  summarize(n = n())
roinot25 <- which(full_dt_zerocountsummarize$n < 25)
ndonorpresent <- c()
for (roi in roinot25){
  donors_present <- full_dt_subset[which(full_dt_subset$chromwindow == full_dt_zerocountsummarize[roi,]$chromwindow), "Donor"]
  ndonorpresent <- c(nrow(donors_present), ndonorpresent)
}
roi1 <- which(full_dt_zerocountsummarize$n == 1)
windowszerocount_1donor <- rbindlist(lapply(1:length(roi1), function(x) full_dt_subset[which(full_dt_subset$chromwindow == full_dt_zerocountsummarize[roi1[x],]$chromwindow),]))

windowszerocount_1donor <- windowszerocount_1donor %>% setorder(Chrom, Window)

lappend <- function (a, lst){
  lst <- c(lst, list(a))
  return(lst)
}


donororder <- c()
ndonorcontig <- c()
donorwindows <- c()
donorwindowstogether <- list()
previous_num <- 1
for (roi in 2:nrow(windowszerocount_1donor)){
  previousChrom <- windowszerocount_1donor$Chrom[roi-1]
  currentChrom <- windowszerocount_1donor$Chrom[roi]
  previousWindow <- windowszerocount_1donor$Window[roi-1]
  donorwindows <- c(windowszerocount_1donor$chromwindow[roi-1], donorwindows)
  currentWindow <- windowszerocount_1donor$Window[roi]
  previousDonor <- windowszerocount_1donor$Donor[roi-1]
  currentDonor <- windowszerocount_1donor$Donor[roi]
  if (previousChrom == currentChrom & (previousWindow+1) == currentWindow & previousDonor == currentDonor){
    previous_num <- previous_num + 1
    donorwindows <- c(windowszerocount_1donor$chromwindow[roi], donorwindows)
  } else{
    donororder <- c(donororder, previousDonor)
    ndonorcontig <- c(ndonorcontig, previous_num)
    donorwindowstogether <- lappend(donorwindows, donorwindowstogether)
    donorwindows <- c()
    previous_num <- 1
  }
}

donorwt <- unlist(lapply(1:length(donorwindowstogether), function(x) str_c(unique(unlist(donorwindowstogether[x])), collapse=',')))
df <- data.frame(donors = donororder, ncontig = ndonorcontig, windows = donorwt)

write.csv(df, "~/scr16_rmccoy22/kweave23/sc_transmission_distortion/contig_donor_windows_withoutHetSNPs_df.csv")

full_dt_windowmap <- full_dt %>% group_by(chromwindow) %>% select(chromwindow, Chrom, WindowStart, WindowEnd) %>% distinct()
write.csv(full_dt_windowmap, "~/scr16_rmccoy22/kweave23/sc_transmission_distortion/chromwindow_map_dt.csv")



