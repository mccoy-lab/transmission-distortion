# Script to run rhapsodi on infertile data from Leung et al. (transferred from Avery to MARCC: `/scratch/groups/rmccoy22/scarios1/bell_infertile_data`)

# Only major change is the read-in of dt. For the original data, we processed it to remove cells of low-quality, according to Avery's study.
# This infertile study has no such companion file, so we use the raw data directly but it's in a different format. 

library(data.table)
library(tidyverse)
library(pbapply)
library(pbmcapply)
library(rhapsodi)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
blacklist_file <- args[2]
giab_file <- args[3]
grch38_file <- args[4]
snps_1kgp_file <- args[5]
sampleName <- args[6]
chrom <- args[7]
outDir <- args[8]
threads <- as.integer(args[9])
seqError <- 0.005
avgRecomb <- 1
window_length <- 3000
smooth_imputed_genotypes <- FALSE
smooth_crossovers <- TRUE

full_data <- read_delim(input_file, delim = "\t") %>%
  pivot_wider(., names_from = "cell", values_from = "gt") %>%
  arrange(., pos) %>%
  as.data.frame()

dt <- full_data

dt <- dt[((rowSums(dt[, 2:ncol(dt)] == 0, na.rm = TRUE) > 0) & (rowSums(dt[, 2:ncol(dt)] == 1, na.rm = TRUE) > 0)),]

# remove the first column (positions)
positions <- dt[, 1]

# Define "not in" function
'%!in%' <- function(x,y)!('%in%'(x,y))

# Get chromosome of interest
col_chr_name <- paste0("chr", chrom)

# 1. Remove problem regions from grch38 T2T
# Make positions into a bed-like file for comparison
positions_dt <- dt[,1] %>% as.data.table()
positions_dt[, chr := col_chr_name]
positions_dt[, start := positions]
positions_dt[, end := positions]
# Read in T2T 100 bp segments data
grch38_poor_map_regions <- read.table(grch38_file, sep = "\t", header = FALSE) %>% as.data.table()
colnames(grch38_poor_map_regions) <- c("chr", "start", "end")
# Get just chr of interest
grch38_poor_map_regions_chr <- grch38_poor_map_regions[grch38_poor_map_regions$chr == col_chr_name]

# Find regions in the raw data that are in the grch38 difficult-to-map sections
setkey(grch38_poor_map_regions_chr, start, end)
raw_data_in_grch38_poor_map <- data.table::foverlaps(positions_dt, grch38_poor_map_regions_chr, type="any", nomatch=NULL)

# Remove those regions from dt
dt_filtered_1 <- positions_dt[positions_dt$start %!in% raw_data_in_grch38_poor_map$start]


# 2. Remove problem regions from ENCODE
# Read in ENCODE blacklist data
blacklist <- read_delim(blacklist_file, col_names = FALSE, delim = "\t") %>%
  as.data.table()
colnames(blacklist) <- c("chr", "start", "end")
# Get just chr of interest
blacklist_chr <- blacklist[blacklist$chr == col_chr_name]
# Get positions in raw data that are in the blacklist
setkey(blacklist_chr, start, end)
raw_data_in_blacklist <- data.table::foverlaps(dt_filtered_1, blacklist_chr, type="any", nomatch=NULL)
# Remove those regions from dt
dt_filtered_2 <- dt_filtered_1[dt_filtered_1$start %!in% raw_data_in_blacklist$start]

# 3. Keep only good regions from GIAB
# Read in GIAB union data
giab_union <- read.table(giab_file, sep = "\t", header = FALSE) %>% as.data.table()
colnames(giab_union) <- c("chr", "start", "end")
# Get just chr of interest
giab_union_chr <- giab_union[giab_union$chr == col_chr_name]
# Get positions in the raw data that are in the giab union
setkey(giab_union_chr, start, end)
raw_data_in_giab_union <- data.table::foverlaps(dt_filtered_2, giab_union_chr, type="any", nomatch=NULL)
colnames(raw_data_in_giab_union) <- c("chr", "start", "end", "positions", "i.chr", "i.start", "i.end")
# In raw data, keep only the positions that are in the union
dt_filtered_3 <- dt_filtered_2[dt_filtered_2$start %in% raw_data_in_giab_union$positions] # 107741

# 4. Remove uncommon variants, per 1kGP
# Read in 1kgp data
snps_1kgp_all <- read.table(snps_1kgp_file, header = FALSE)
colnames(snps_1kgp_all) <- c("chr", "positions")
# get just chr of interest
snps_1kgp_chr <- snps_1kgp_all[snps_1kgp_all$chr == col_chr_name,]
# Make bed-like file for 1kgp snps
snps_1kgp_chr_bed <- snps_1kgp_chr$positions %>% as.data.table()
snps_1kgp_chr_bed[, chr := col_chr_name]
snps_1kgp_chr_bed[, start := snps_1kgp_chr$positions]
snps_1kgp_chr_bed[, end := snps_1kgp_chr$positions]
# Get positions in raw data that are in the 1kgp SNP list
setkey(snps_1kgp_chr_bed, start, end)
raw_data_in_1kgp_snps <- data.table::foverlaps(dt_filtered_3, snps_1kgp_chr_bed, type="any", nomatch=NULL)
colnames(raw_data_in_1kgp_snps) <- c("position", "chr", "start", "end", "positions", "i.chr", "i.start", "i.end")
# In raw data, keep only the positions that are in the 1kgp list
dt_filtered_4 <- dt_filtered_3[dt_filtered_3$start %in% raw_data_in_1kgp_snps$positions] # 107333

# Use the newly filtered set of SNPs to subset the raw data
dt <- dt[dt$positions %in% dt_filtered_4$start,] # 107333
# Remove any SNPs that are not heterozygous
dt <- dt[((rowSums(dt[, 2:ncol(dt)] == 0, na.rm = TRUE) > 0) & (rowSums(dt[, 2:ncol(dt)] == 1, na.rm = TRUE) > 0)),] # 104324

message("good before step5")
#5. Remove SNPs that have more reads than expected
allele_counts <- data.table(x = rowSums(dt[, 2:ncol(dt)] == 1, na.rm = TRUE),
                            y = rowSums(dt[, 2:ncol(dt)] == 0, na.rm = TRUE))
allele_counts$positions <- dt[,1]
allele_counts$sum <- allele_counts$x + allele_counts$y
mean_read_count <- mean(allele_counts$x + allele_counts$y)
std_dev_read_count <- sqrt(var(allele_counts$x + allele_counts$y))
# Remove any SNPs that have a read count more than a standard deviation above the mean
cutoff_read_count <- mean_read_count + std_dev_read_count
allele_counts_high <- allele_counts[allele_counts$sum > cutoff_read_count]# this gets rid of 12.9% --> 13,490
dt_filtered_5 <- dt[dt$positions %!in% allele_counts_high$positions,]#90,834

# Use dt for downstream analyses
dt <- dt_filtered_5

filename_dt <- paste0(outDir, sampleName, "_", chrom, "_full_filtered_dt.csv")
write_csv(dt, filename_dt)

#Run rhapsodi
standard_input_out <- rhapsodi::standard_input(NULL, use_dt = TRUE, input_dt = dt)
complete_haplotypes <- rhapsodi::impute_donor_haplotypes(standard_input_out$dt, standard_input_out$positions, threads = threads, window_length = window_length)
filled_gametes <- rhapsodi::fill_gametes(standard_input_out$dt, complete_haplotypes, threads = threads, sequencing_error = seqError, avg_recomb = avgRecomb)
rhapsodi_out <- rhapsodi::report_gametes(smooth_crossovers, smooth_imputed_genotypes, complete_haplotypes, standard_input_out$dt, filled_gametes, standard_input_out$positions, sampleName, chrom, threads=threads)

#Report Donor Haps
filename_hap <- paste0(outDir, sampleName, "_", chrom, "_parental_hap.csv")
write_csv(rhapsodi_out$donor_haps, filename_hap)

#Report filled sperm
filename_fs <- paste0(outDir, sampleName, "_", chrom, "_filled_sperm_smoothed.csv")
write_csv(filled_gametes, filename_fs)
filename_fs <- paste0(outDir, sampleName, "_", chrom, "_filled_sperm_unsmoothed.csv")
write_csv(rhapsodi_out$gamete_haps, filename_fs)

# Report recombination
filename_rs <- paste0(outDir, sampleName, "_", chrom, "_recombination_locs.csv")
write_csv(rhapsodi_out$recomb_breaks, filename_rs)

td_test <- function(sperm_matrix, row_index) {
  test_row <- sperm_matrix[row_index,]
  gt_vector <- unlist(test_row)[-1]
  one_count <- sum(gt_vector == "haplotype1", na.rm = TRUE)
  two_count <- sum(gt_vector == "haplotype2", na.rm = TRUE)
  p_value <- binom.test(c(one_count, two_count))$p.value
  return(c(p_value, one_count, two_count))
}

df_counts_pvals <- do.call(rbind, pbmclapply(1:nrow(rhapsodi_out$gamete_haps),
                                             function(x) td_test(rhapsodi_out$gamete_haps, x),
                                             mc.cores=getOption("mc.cores", threads))) %>%
  as_tibble() %>%
  add_column(standard_input_out$positions) #bind the positions vector to df_counts_pvals
colnames(df_counts_pvals) <- c("pval", "h1_count", "h2_count", "genomic_position")