library(data.table)
library(tidyverse)
library(pbapply)
library(pbmcapply)
library(rhapsodi)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
blacklist_file <- args[2]
giab_file <- args[3]
snps_1kgp_file <- args[4] 
sampleName <- args[5]
chrom <- args[6]
outDir <- args[7]
threads <- as.integer(args[8])

dt <- read_delim(input_file, delim = "\t", col_types = cols(.default = "d")) %>%
  as.data.frame()

orig_snps <- nrow(dt)
num_gams <- ncol(dt)-1

#0. Select only heterozygous SNPs
dt <- dt[((rowSums(dt[, 2:ncol(dt)] == 0, na.rm = TRUE) > 0) & (rowSums(dt[, 2:ncol(dt)] == 1, na.rm = TRUE) > 0)),]
snps0 <- nrow(dt)

# Define "not in" function
'%!in%' <- function(x,y)!('%in%'(x,y))

# Get chromosome of interest
col_chr_name <- paste0("chr", chrom)

# Make positions into a bed-like file for comparison
positions <- dt[,1]
positions_dt <- dt[,1] %>% as.data.table()
positions_dt[, chr := col_chr_name]
positions_dt[, start := positions]
positions_dt[, end := positions]

# 1. Remove problem regions from ENCODE
# Read in ENCODE blacklist data
blacklist <- read_delim(blacklist_file, col_names = FALSE, delim = "\t") %>%
  as.data.table()
colnames(blacklist) <- c("chr", "start", "end")
# Get just chr of interest
blacklist_chr <- blacklist[blacklist$chr == col_chr_name]
# Get positions in raw data that are in the blacklist
setkey(blacklist_chr, start, end)
raw_data_in_blacklist <- data.table::foverlaps(positions_dt, blacklist_chr, type="any", nomatch=NULL)
# Remove those regions from dt
dt_filtered_1 <- positions_dt[positions_dt$start %!in% raw_data_in_blacklist$start]
snps1 <- nrow(dt_filtered_1)

# 2. Keep only good regions from GIAB
# Read in GIAB union data
giab_union <- read.table(giab_file, sep = "\t", header = FALSE) %>% as.data.table()
colnames(giab_union) <- c("chr", "start", "end")
# Get just chr of interest
giab_union_chr <- giab_union[giab_union$chr == col_chr_name]
# Get positions in the raw data that are in the giab union
setkey(giab_union_chr, start, end)
raw_data_in_giab_union <- data.table::foverlaps(dt_filtered_1, giab_union_chr, type="any", nomatch=NULL)
colnames(raw_data_in_giab_union) <- c("chr", "start", "end", "positions", "i.chr", "i.start", "i.end")
# In raw data, keep only the positions that are in the union
dt_filtered_2 <- dt_filtered_1[dt_filtered_1$start %in% raw_data_in_giab_union$positions]
snps2 <- nrow(dt_filtered_2)


# 3. Remove uncommon variants, per 1kGP
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
raw_data_in_1kgp_snps <- data.table::foverlaps(dt_filtered_2, snps_1kgp_chr_bed, type="any", nomatch=NULL)
colnames(raw_data_in_1kgp_snps) <- c("position", "chr", "start", "end", "positions", "i.chr", "i.start", "i.end")
# In raw data, keep only the positions that are in the 1kgp list
dt_filtered_3 <- dt_filtered_2[dt_filtered_2$start %in% raw_data_in_1kgp_snps$positions]
snps3 <- nrow(dt_filtered_3)

# Use the newly filtered set of SNPs to subset the raw data
dt <- dt[dt$positions %in% dt_filtered_3$start,]

#4. Remove SNPs that have more reads than expected
allele_counts <- data.table(x = rowSums(dt[, 2:ncol(dt)] == 1, na.rm = TRUE),
                            y = rowSums(dt[, 2:ncol(dt)] == 0, na.rm = TRUE))
allele_counts$positions <- dt[,1]
allele_counts$sum <- allele_counts$x + allele_counts$y
mean_read_count <- mean(allele_counts$x + allele_counts$y)
std_dev_read_count <- sqrt(var(allele_counts$x + allele_counts$y))
# Remove any SNPs that have a read count more than a standard deviation above the mean
cutoff_read_count <- mean_read_count + std_dev_read_count
allele_counts_high <- allele_counts[allele_counts$sum > cutoff_read_count]
dt_filtered_4 <- dt[dt$positions %!in% allele_counts_high$positions,]
snps4 <- nrow(dt_filtered_4)

# Use dt_filtered_4 for downstream analyses
filename_dt <- paste0(outDir, sampleName, "_", chrom, "_full_filtered_dt.csv")
write_csv(dt_filtered_4, filename_dt)

#5. Run rhapsodi
rhapsodi_out <- rhapsodi::rhapsodi_autorun(NULL, use_dt = TRUE, input_dt = dt_filtered_4, threads = threads, mcstop=FALSE, verbose = TRUE, sampleName = sampleName, chrom = chrom)

##Report Donor Haps
filename_hap <- paste0(outDir, sampleName, "_", chrom, "_parental_hap.csv")
write_csv(rhapsodi_out$donor_haps, filename_hap)

##Report filled sperm
filename_fs <- paste0(outDir, sampleName, "_", chrom, "_filled_sperm_smoothed.csv")
write_csv(rhapsodi_out$gamete_haps, filename_fs)

filename_fs <- paste0(outDir, sampleName, "_", chrom, "_filled_sperm_unsmoothed.csv")
write_csv(rhapsodi_out$unsmoothed_gamete_haps, filename_fs)

##Report recombination
filename_rs <- paste0(outDir, sampleName, "_", chrom, "_recombination_locs.csv")
write_csv(rhapsodi_out$recomb_breaks, filename_rs)

#6. TD scan
td_test <- function(sperm_matrix, row_index) {
  test_row <- sperm_matrix[row_index,]
  gt_vector <- unlist(test_row)[-1]
  one_count <- sum(gt_vector == "h1", na.rm = TRUE)
  two_count <- sum(gt_vector == "h2", na.rm = TRUE)
  p_value <- binom.test(c(one_count, two_count))$p.value
  return(c(p_value, one_count, two_count))
}

df_counts_pvals <- do.call(rbind, pbmclapply(1:nrow(rhapsodi_out$unsmoothed_gamete_haps),
                                             function(x) td_test(rhapsodi_out$unsmoothed_gamete_haps[,-c(1,2)], x),
                                             mc.cores=getOption("mc.cores", threads))) %>%
  as_tibble() %>%
  add_column(rhapsodi_out$unsmoothed_gamete_haps[,2]) #bind the positions vector to df_counts_pvals
colnames(df_counts_pvals) <- c("pval", "h1_count", "h2_count", "genomic_position")
filename_df <- paste0(outDir, sampleName, "_", chrom, "_pval.csv")
write_csv(df_counts_pvals, filename_df)

#7. Look at some metadata for the data and output

fake_true_phasing <- matrix(rep(0, snps4*2), ncol = 2, nrow=snps4) %>% as.data.frame() %>% `colnames<-`(c("donor1", "donor2"))
fake_true_gametes <- cbind(1:snps4, matrix(rep(0, snps4 * num_gams), ncol=num_gams, nrow=snps4)) %>% as.data.frame() %>% `colnames<-`(c("positions", colnames(rhapsodi_out$gamete_genotypes[,-c(1,2)])))

#Input
## Input COM by gamete
input_assess <- rhapsodi::sim_assess_gam_imputation(fake_true_gametes, dt_filtered_4[,-1], snps4, num_gams)

## Input approx COV by gamete

input_cov <- -log(1 - input_assess$com)

## First non-NA SNP by gamete

first_genomic_loc_by_gam <- unlist(lapply(1:ncol(dt_filtered_4[,-1]), function(x) dt_filtered_4[,1][which(!is.na(dt_filtered_4[,x+1]))[1]])) #%>% `names<-`(colnames(dt_filtered_4[,-1]))
first_index_loc_by_gam <- unlist(lapply(1:ncol(dt_filtered_4[,-1]), function(x) which(!is.na(dt_filtered_4[,x+1]))[1])) #%>% `names<-`(colnames(dt_filtered_4[,-1]))

## Last non-NA SNP by gamete

last_genomic_loc_by_gam <- unlist(lapply(1:ncol(dt_filtered_4[,-1]), function(x) dt_filtered_4[,1][rev(which(!is.na(dt_filtered_4[,x+1])))[1]])) #%>% `names<-`(colnames(dt_filtered_4[,-1]))
last_index_loc_by_gam <- unlist(lapply(1:ncol(dt_filtered_4[,-1]), function(x) rev(which(!is.na(dt_filtered_4[,x+1])))[1])) #%>% `names<-`(colnames(dt_filtered_4[,-1]))

gamete_input_metadata <- data.frame(gamete = colnames(dt_filtered_4[,-1]),
                                    donor = sampleName, chrom = chrom,
                                    input_com = input_assess$com,
                                    input_approx_cov = input_cov,
                                    first_nonNA_pos = first_genomic_loc_by_gam ,
                                    first_nonNA_index = first_index_loc_by_gam,
                                    last_nonNA_pos = last_genomic_loc_by_gam, 
                                    last_nonNA_index = last_index_loc_by_gam)

#Phasing
## Phasing COM

phasing_assess <- rhapsodi::sim_assess_phasing(fake_true_phasing, rhapsodi_out$donor_haps, snps4)

#Imputation
## Imputation COM by gamete

imputation_assess <- rhapsodi::sim_assess_gam_imputation(fake_true_gametes, rhapsodi_out$gamete_genotypes[,-c(1,2)], snps4, num_gams)

## Differences between smoothed and unsmoothed by gamete

smooth_unsmooth_dif_by_gam <- colSums(rhapsodi_out$gamete_genotypes[,-c(1,2)] != rhapsodi_out$unsmoothed_gamete_genotypes[,-c(1,2)], na.rm=TRUE)

gamete_genotype_metadata <- data.frame(gamete = colnames(rhapsodi_out$gamete_genotypes[,-c(1,2)]),
                                       output_com = imputation_assess$com,
                                       unsmoothed_ndif = smooth_unsmooth_dif_by_gam) 

#Discovery
nonrecomb_chroms <- rhapsodi_out$recomb_breaks[is.na(rhapsodi_out$recomb_breaks$Genomic_start),]$Ident
num_nonrecomb <- length(nonrecomb_chroms)
recomb_breaks <- rhapsodi_out$recomb_breaks[!is.na(rhapsodi_out$recomb_breaks$Genomic_start),]
recomb_out <- data.frame(Ident = recomb_breaks$Ident, res = recomb_breaks$Genomic_end - recomb_breaks$Genomic_start)


## Number of crossovers by gamete
recomb_summary <- table(recomb_breaks$Ident) %>%
  as.data.table() %>%
  setnames(., c("donor_chrom_gamete", "n_crossovers")) %>%
  rbind(data.table(donor_chrom_gamete = nonrecomb_chroms, n_crossovers = 0)) %>%
  setorder(n_crossovers)

recomb_summary[, c("donor", "chrom", "gamete") := tstrsplit(donor_chrom_gamete, "_", fixed = TRUE)]



#8. Report the metadata stats
sample_metadata_out <- data.frame(donor=sampleName, chrom = chrom, ngam = num_gams, 
                           nsnp_in = orig_snps, nsnp_het = snps0,
                           nsnp_encode = snps1, nsnp_giab = snps2, 
                           nsnp_1kgp = snps3, nsnp_segdup = snps4,
                           phasing_com = phasing_assess$com, 
                           num_recombinant_gam = num_gams - num_nonrecomb,
                           num_non_recombinant_gam = num_nonrecomb)

sample_metadata_df_file <- paste0(outDir, sampleName, "_", chrom, "_sample_meta.csv")
write_csv(sample_metadata_out, sample_metadata_df_file)

gamete_metadata_out <- merge(gamete_input_metadata, gamete_genotype_metadata, by = "gamete") #join gamete_input_metadata and gamete_genotype_metadata by gamete
gamete_metadata_out <- merge(gamete_metadata_out, recomb_summary, by=c("gamete", "donor", "chrom")) #then join this with recomb_summary's n_crossovers by gamete

gamete_metadata_df_file <- paste0(outDir, sampleName, "_", chrom, "_gamete_meta.csv")
write_csv(gamete_metadata_out, gamete_metadata_df_file)

recomb_out_file <- paste0(outDir, sampleName, "_", chrom, "_recomb_res.csv")
write_csv(recomb_out, recomb_out_file)

                          

