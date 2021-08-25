library(data.table)
library(tidyverse)
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
hapProb <- 1 - seqError
window_length <- 3000
smooth_imputed_genotypes <- FALSE
smooth_crossovers <- TRUE

dt <- read_delim(input_file, delim = "\t", col_types = cols(.default = "d")) %>% 
  as.data.frame()

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

dt <- dt[, -1]


getmode <- function(x) { #from https://stackoverflow.com/questions/56552709/r-no-mode-and-exclude-na?noredirect=1#comment99692066_56552709
  ux <- unique(na.omit(x))
  tx <- tabulate(match(x, ux))
  if(length(ux) != 1 & sum(max(tx) == tx) > 1) {
    if (is.character(ux)) return(NA_character_) else return(NA_real_)
  }
  max_tx <- tx == max(tx)
  return(ux[max_tx])
}

invertBits <- function(df) {
  return(abs(df-1))
}

# overlapping window function from https://stackoverflow.com/questions/8872376/split-vector-with-overlapping-samples-in-r
splitWithOverlap <- function(vec, seg.length, overlap) {
  starts = seq(1, length(vec), by=seg.length-overlap)
  ends   = starts + seg.length - 1
  ends[ends > length(vec)] = length(vec)
  lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
}

# use overlaps of window length/2
windows <- splitWithOverlap(rank(positions), window_length, overlap = window_length / 2)

# merge the last two windows to avoid edge effect
if (length(windows) > 1){
  combined <- unique(c(windows[[length(windows) - 1]], windows[[length(windows)]]))
  combined <- combined[order(combined)]
  total_combined <- windows[-c((length(windows) - 1), length(windows))]
  total_combined[[length(total_combined) + 1]] <- combined
  windows <- total_combined
}

num_snps <- nrow(dt)
message(paste0("Number of windows with overlap of ", window_length / 2 , " and ", num_snps, " number of SNPs following filtering: ", length(windows)))


# function to reconstruct parental haplotypes
reconstruct_hap <- function(input_dt, input_positions, window_indices) {
  window_start <- min(window_indices)
  window_end <- max(window_indices)
  positions_for_window <- input_positions[window_start:window_end]
  # compute a distance matrix
  d <- dist(t(as.matrix(input_dt)[window_start:window_end,]), method = "binary")
  # plut in 0.5 for any NA entries of the distance matrix
  d[is.na(d)] <- 0.5
  # cluster the distance matrix
  tree <- hclust(d, method = "ward.D2")
  # plot(tree, cex = 0.1) # uncomment to plot
  # cut the tree generated by clustering into two groups (haplotypes)
  haplotypes <- cutree(tree, k=2)
  # get the names of the sperm cells falling into the two groups
  h1_sperm <- names(haplotypes[haplotypes == 1])
  h2_sperm <- names(haplotypes[haplotypes == 2])
  # reconstruct the original haplotypes by majority vote after inverting the opposite haplotype
  h1_inferred <- unname(apply(cbind(input_dt[window_start:window_end, h1_sperm],
                                    invertBits(input_dt[window_start:window_end, h2_sperm])),
                              1, function(x) getmode(x)))
  # h2_inferred <- unname(apply(cbind(input_dt[window_start:window_end, h2_sperm],
  #                                   invertBits(input_dt[window_start:window_end, h1_sperm])),
  #                             1, function(x) getmode(x)))
  return(tibble(index = window_indices, pos = positions_for_window, h1 = h1_inferred))
}

# infer the haplotypes within the overlapping windows
before_time <- Sys.time()
inferred_haplotypes <- pbmclapply(1:length(windows), 
                                  function(x) reconstruct_hap(dt, positions, windows[[x]]),
                                  mc.cores = getOption("mc.cores", threads))

# stitch together the haplotypes
initial_haplotype <- inferred_haplotypes[[1]]
for (hap_window in 1:length(windows)) {
  olap_haps <- merge(initial_haplotype, inferred_haplotypes[[hap_window]], by = "index")
  olap_haps_complete <- merge(initial_haplotype, inferred_haplotypes[[hap_window]], by = "index", all = TRUE)
  #mean_concordance <- mean(olap_haps$h1.x == olap_haps$h1.y)
  mean_concordance2 <- mean(olap_haps$h1.x == olap_haps$h1.y, na.rm = TRUE)
  message(paste("Window:", hap_window, mean_concordance2, sep = " "))
  if (mean_concordance2 < 0.1) {
    olap_haps_complete$h1.y <- invertBits(olap_haps_complete$h1.y)
  } else if (mean_concordance2 < 0.9) {
    stop(paste0("Haplotypes within overlapping windows are too discordant to merge. Mean: ", mean_concordance2, " and Window: ", hap_window))
  }
  initial_haplotype <- tibble(index = olap_haps_complete$index,
                              pos = c(olap_haps_complete[is.na(olap_haps_complete$pos.y),]$pos.x,
                                      olap_haps_complete[!is.na(olap_haps_complete$pos.x) &
                                                         !is.na(olap_haps_complete$pos.y),]$pos.x,
                                      olap_haps_complete[is.na(olap_haps_complete$pos.x),]$pos.y),
                              h1 = c(olap_haps_complete[is.na(olap_haps_complete$pos.y),]$h1.x,
                                     olap_haps_complete[!is.na(olap_haps_complete$pos.x) &
                                                        !is.na(olap_haps_complete$pos.y),]$h1.x,
                                     olap_haps_complete[is.na(olap_haps_complete$pos.x),]$h1.y))
}
complete_haplotypes <- initial_haplotype %>%
  mutate(h2 = invertBits(h1))
after_time <- Sys.time()
time_to_impute <- difftime(after_time, before_time, units="mins")
message(paste0("Time used for inferring parental haplotypes: ", time_to_impute))

filename_hap <- paste0(outDir, sampleName, "_", chrom, "_parental_hap.csv")
write_csv(complete_haplotypes, filename_hap)

# Going through each sperm, if an allele (0 or 1) in a sperm matches the allele (0 or 1)
# in h1 at that position, replace the allele with "h1". Do the same for h2.
for (i in 1:ncol(dt)) {
  dt[i][dt[i] == complete_haplotypes$h1] <- "h1"
  dt[i][dt[i] == complete_haplotypes$h2] <- "h2"
  dt[c(which(dt[,i] == 0 | dt[,i] == 1)),i] <- NA
}

# Scan sperm by sperm to interpret state given emission
# First, we initialize our HMM
# set denominator for transition probability - one recombination event per chromosome
num_snps <- nrow(complete_haplotypes)
# two states
states <- c("haplotype1", "haplotype2")
# probability of state at position x+1 given state at position x
hap1Prob <- c(1-(1/num_snps), 1/num_snps)
hap2Prob <- c(1/num_snps, 1-(1/num_snps))
transProb <- matrix(c(hap1Prob, hap2Prob), 2)
# Two emissions (observations): an allele from h1 or an allele from h2
emissions <- c("h1","h2")
# Prob of emitting an h1 allele, prob of emitting an h2 allele in state `haplotype1`
h1ProbEmiss <- c(hapProb, seqError)
# Prob of emitting an h1 allele, prob of emitting an h2 allele in state `haplotype2`
h2ProbEmiss <- c(seqError, hapProb)
emissProb <- matrix(c(h1ProbEmiss, h2ProbEmiss), 2)
#build model with the above inputs
hmm <- initHMM(States = states,
               Symbols = emissions,
               transProbs = transProb,
               emissionProbs = emissProb)
###### Function to run HMM on an input file
# Compute the inferred state using each sperm cell as the input
# (Input must be a vector)
runHMM <- function(sperm_dt, column_index) {
  original_obs <- sperm_dt[,column_index]
  inferred_state <- viterbi(hmm, na.omit(sperm_dt[, column_index]))
  original_obs[!is.na(original_obs)] <- inferred_state
  return(original_obs)
}

before_impute <- Sys.time()
imputed_sperm <- as_tibble(do.call(cbind, pbmclapply(1:ncol(dt),
                                                     function(x) runHMM(dt, x),
                                                     mc.cores = getOption("mc.cores", threads))))
after_impute <- Sys.time()
time_to_impute <- difftime(after_impute, before_impute, units = "mins")
message(paste0("Time used for imputing sperm haplotypes: ", time_to_impute))
colnames(imputed_sperm) <- colnames(dt)

# Works on our sperm! Need to make the function work on every sperm in test3
# and need to make fill up and down at the end
fill_NAs <- function(merged_sperm, col_index) {
  sperm_sample <- merged_sperm[,col_index] %>%
    rename(sperm = colnames(.)[1]) %>%
    mutate(sperm_up = sperm) %>%
    mutate(sperm_down = sperm) %>%
    fill(sperm_up, .direction = "up") %>%
    fill(sperm_down, .direction = "down") %>%
    mutate(is_match = (sperm_up == sperm_down)) %>%
    replace_na(list(is_match = FALSE))
  sperm_sample$sperm_imputed <- as.character(NA)
  sperm_sample[sperm_sample$is_match == TRUE,]$sperm_imputed <- sperm_sample[sperm_sample$is_match == TRUE,]$sperm_up
  #fill beginning of chromosome NA's
  first <- which(!is.na(sperm_sample$sperm_imputed))[1]
  sperm_sample$sperm_imputed[1:(first-1)] <- sperm_sample$sperm_imputed[first]
  #fill end of chromosome NA's
  sperm_sample$sperm_imputed <- rev(sperm_sample$sperm_imputed)
  first <- which(!is.na(sperm_sample$sperm_imputed))[1]
  sperm_sample$sperm_imputed[1:(first-1)] <- sperm_sample$sperm_imputed[first]
  #reverse chromosome imputation back so it faces the right way
  sperm_sample$sperm_imputed <- rev(sperm_sample$sperm_imputed)
  return(sperm_sample$sperm_imputed)
}

filled_sperm <- as_tibble(do.call(cbind, 
                                  pblapply(1:ncol(imputed_sperm),
                                           function(x) fill_NAs(imputed_sperm, x))))
colnames(filled_sperm) <- colnames(dt)
filename_fs <- paste0(outDir, sampleName, "_", chrom, "_filled_sperm_smoothed.csv")
write_csv(filled_sperm, filename_fs)


unsmooth <- function(original_gamete_df, filled_gamete_data){
  original_gamete_df[original_gamete_df == "h1"] <- "haplotype1"
  original_gamete_df[original_gamete_df == "h2"] <- "haplotype2"
  original_dt <- as.data.frame(original_gamete_df)
  filled_gamete_data <- as.data.frame(filled_gamete_data)
  filled_gamete_data[!is.na(original_dt)] <- original_dt[!is.na(original_dt)]
  filled_gamete_data <- as_tibble(filled_gamete_data)
  return (filled_gamete_data)
}

#find recombination spots
find_recomb_spots <- function(input_matrix, x, identities, genomic_positions){
  ident <- identities[x]
  input_tibble <- input_matrix[, x] %>%
    mutate(., index = row_number()) %>%
    mutate(., positions = genomic_positions)
  complete_cases_tibble <- input_tibble[complete.cases(input_tibble),]
  input_vec <- as.factor(complete_cases_tibble[[1]])
  switch_indices <- which(input_vec[-1] != input_vec[-length(input_vec)])
  switch_indices_input <- complete_cases_tibble[switch_indices,]$index
  crossover_start <- input_tibble[switch_indices_input,]$positions
  rev_input_tibble <- arrange(input_tibble, -index) %>%
    mutate(., index = row_number())
  complete_cases_rev_tibble <- rev_input_tibble[complete.cases(rev_input_tibble),]
  rev_input_vec <- as.factor(complete_cases_rev_tibble[[1]])
  rev_switch_indices <- which(rev_input_vec[-1] != rev_input_vec[-length(rev_input_vec)])
  rev_switch_indices_input <- complete_cases_rev_tibble[rev_switch_indices,]$index
  crossover_end <- rev(rev_input_tibble[rev_switch_indices_input,]$positions)
  recomb_spots <- tibble(Ident = ident, Genomic_start = crossover_start, Genomic_end = crossover_end)
  return(recomb_spots)
}

if (!smooth_crossovers){
  filled_sperm_recomb <- unsmooth(dt, filled_sperm)
  idents_for_csv <- paste0(paste0(sampleName, "_", chrom, "_"), colnames(filled_sperm_recomb))
  recomb_spots_all <- do.call(rbind, pbmclapply(1:ncol(filled_sperm_recomb),
                                                function(x) find_recomb_spots(filled_sperm_recomb, x, idents_for_csv, positions),
                                                mc.cores=getOption("mc.cores", threads))) %>%
    right_join(., tibble(Ident = idents_for_csv), by = "Ident")
} else {
  idents_for_csv <- paste0(paste0(sampleName, "_", chrom, "_"), colnames(filled_sperm))
  recomb_spots_all <- do.call(rbind, pbmclapply(1:ncol(filled_sperm),
                                                function(x) find_recomb_spots(filled_sperm, x, idents_for_csv, positions),
                                                mc.cores=getOption("mc.cores", threads))) %>%
    right_join(., tibble(Ident = idents_for_csv), by = "Ident")
}

filename_rs <- paste0(outDir, sampleName, "_", chrom, "_recombination_locs.csv")
write_csv(recomb_spots_all, filename_rs)

if (!smooth_imputed_genotypes & smooth_crossovers){
  filled_sperm <- unsmooth(dt, filled_sperm)
  filename_fs <- paste0(outDir, sampleName, "_", chrom, "_filled_sperm_unsmoothed.csv")
  write_csv(filled_sperm, filename_fs)
} else if(!smooth_imputed_genotypes & !smooth_crossovers){
  filled_sperm <- filled_sperm_recomb
}


td_test <- function(sperm_matrix, row_index) {
  test_row <- sperm_matrix[row_index,]
  gt_vector <- unlist(test_row)[-1]
  one_count <- sum(gt_vector == "haplotype1", na.rm = TRUE)
  two_count <- sum(gt_vector == "haplotype2", na.rm = TRUE)
  p_value <- binom.test(c(one_count, two_count))$p.value
  return(c(p_value, one_count, two_count))
}

df_counts_pvals <- do.call(rbind, pbmclapply(1:nrow(filled_sperm),
                                             function(x) td_test(filled_sperm, x),
                                             mc.cores=getOption("mc.cores", threads))) %>%
  as_tibble() %>%
  add_column(positions) #bind the positions vector to df_counts_pvals
colnames(df_counts_pvals) <- c("pval", "h1_count", "h2_count", "genomic_position")
filename_df <- paste0(outDir, sampleName, "_", chrom, "_pval.csv")
write_csv(df_counts_pvals, filename_df)

