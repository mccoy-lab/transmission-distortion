# Investigate long stretches on each chr where one donor has no SNPs (but the rest do) 
library(dplyr)
library(tidyverse)
"
`contig_donor_windows_withoutHetSNPs_df.csv` is the file that lists which 0.5 Mbp windows have only one donor without hetSNPs in that window 
(for the fully filtered data). Col1 is donor. Col2 is number of contiguous 0.5Mbp windows that meet that criteria. 
Col3 is the chromwindow ID for those contiguous windows.

`chromwindow_map_dt.csv` is a dataframe that has the chromwindow ID in the first column, the chromosome in the second, 
the window bp start in the third and the window bp end in the fourth
"

# read in info with windows where only one donor has a stretch of homozygosity
contig_donor_windows_withoutHetSNPs_df <- read.csv("/Users/saracarioscia/mccoy-lab/transmission-distortion/investigate_stretches_homozygosity/contig_donor_windows_withoutHetSNPs_df.csv")
# read in info with windows genomic coordinates 
chromwindow_map_dt <- read.csv("/Users/saracarioscia/mccoy-lab/transmission-distortion/investigate_stretches_homozygosity/chromwindow_map_dt.csv")

# characterize the bins that stretch 7 or more contigs 
contig_donor_windows_withoutHetSNPs_df_long <- contig_donor_windows_withoutHetSNPs_df[contig_donor_windows_withoutHetSNPs_df$ncontig > 6,]
# add column for affected chromosome 
contig_donor_windows_withoutHetSNPs_df_long$chrom <- c(11, 4, 4, 4, 6)
# add column for min window 
contig_donor_windows_withoutHetSNPs_df_long$min_window <- sapply(strsplit(contig_donor_windows_withoutHetSNPs_df_long$windows, ",", fixed=TRUE), tail, 1)
# add column for max window 
contig_donor_windows_withoutHetSNPs_df_long$max_window <- sapply(strsplit(contig_donor_windows_withoutHetSNPs_df_long$windows, ",", fixed=TRUE), head, 1)


# make columns with the same name as the other datatable to allow merges 
chromwindow_map_dt$min_window <- chromwindow_map_dt$chromwindow
chromwindow_map_dt$max_window <- chromwindow_map_dt$chromwindow


# add column for WindowStart and WindowEnd based on window 
contig_donor_windows_withoutHetSNPs_df_long_full <- full_join(contig_donor_windows_withoutHetSNPs_df_long, chromwindow_map_dt, by = "min_window")
# make column for start of large ncontig windows 
contig_donor_windows_withoutHetSNPs_df_long_full$contig_start <- contig_donor_windows_withoutHetSNPs_df_long_full$WindowStart
# get WindowEnd to be the end of the window (i.e., the end of the largest sub window)
contig_donor_windows_withoutHetSNPs_df_long_full_max <- full_join(contig_donor_windows_withoutHetSNPs_df_long, chromwindow_map_dt, by = "max_window")
# make column for end of large ncontig windows 
contig_donor_windows_withoutHetSNPs_df_long_full$contig_end <- contig_donor_windows_withoutHetSNPs_df_long_full_max$WindowEnd
# keep only the 5 rows of interest 
contig_donor_windows_withoutHetSNPs_df_long_full_interest <- contig_donor_windows_withoutHetSNPs_df_long_full[1:5,]
# get column of contig length, to confirm 
contig_donor_windows_withoutHetSNPs_df_long_full_interest$contig_length <- contig_donor_windows_withoutHetSNPs_df_long_full_interest$contig_end - contig_donor_windows_withoutHetSNPs_df_long_full_interest$contig_start
