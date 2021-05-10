library(data.table)
library(tidyverse)
#library(bedr)
library(stringr)
library(pbapply)
library(pbmcapply)
library(HMM)
library(ggplot2)

# args <- commandArgs(trailingOnly = TRUE)
# sampleName <- "sim3"
# chrom <- "chrT"
# outdir <- args[1]
# threads <- 8
# seqError <- 0.05
# hapProb <- 1 - seqError
# hmm_avg_recomb <- 1
# num_gametes <- as.integer(args[2])
# num_snps <- as.integer(args[3])
# coverage <- as.numeric(args[4])
# message(paste0("The coverage of this simulation is: ", coverage))
# num_genotypes <- num_gametes * num_snps
# missing_genotype_rate <- dpois(0, coverage)
# message(paste0("The missing genotype rate of this simulation is ", missing_genotype_rate))
# num_nas <- floor(num_genotypes * missing_genotype_rate)
# window_length <- 3000
# overlap_denom <- 2
# random_seed <- as.integer(args[5])
# set.seed(random_seed)
# recomb_lambda <- 1
# stopifnot(recomb_lambda > 0)
# num_recomb_sites <- rpois(num_gametes, recomb_lambda)
# message(paste0("Total number of recombination spots across gametes: ", sum(num_recomb_sites)))
# add_seq_error <- TRUE
# seqError_add <- 0.05
# add_de_novo_mut <- TRUE
# de_novo_lambda <- 5
# de_novo_alpha <- 7.5
# de_novo_beta <- 10
# stopifnot(de_novo_beta > 0)
# smooth_imputed_genotypes <- FALSE
# smooth_crossovers <- TRUE
# write_out_plot <- FALSE


###Argument inputs manual
sampleName <- "simTest"
chrom <- "chrT"
stopifnot(grepl("chr", chrom, fixed=TRUE))
outdir <- "~/tmp/"
threads <- 2L
seqError <- 0.005
hapProb <- 1 - seqError
hmm_avg_recomb <- 1
num_gametes <- 1000
num_snps <- 5000
coverage <- 0.1
message(paste0("The coverage of this simulation is: ", coverage))
num_genotypes <- num_gametes * num_snps
missing_genotype_rate <- dpois(0, coverage)
message(paste0("The missing genotype rate of this simulation is ", missing_genotype_rate))
num_nas <- floor(num_genotypes * missing_genotype_rate)
window_length <- 3000
overlap_denom <- 2
random_seed <- 42
set.seed(random_seed)
recomb_lambda <- 1
stopifnot(recomb_lambda > 0)
num_recomb_sites <- rpois(num_gametes, recomb_lambda)
message(paste0("Total number of recombination spots across gametes: ", sum(num_recomb_sites)))
add_seq_error <- TRUE
seqError_add <- 0.005
add_de_novo_mut <- FALSE
de_novo_lambda <- 5
de_novo_alpha <- 7.5
de_novo_beta <- 10
stopifnot(de_novo_beta > 0)
smooth_imputed_genotypes <- FALSE
smooth_crossovers <- TRUE
write_out_plot <- FALSE
###

# ###Argument inputs command line 
# args <- commandArgs(trailingOnly = TRUE)
# sampleName <- args[1]
# chrom <- args[2]
# stopifnot(grepl("chr", chrom, fixed=TRUE))
# outDir <- args[3]
# threads <- as.integer(args[4])
# seqError <- as.numeric(args[5])
# hapProb <- 1 - seqError
# hmm_avg_recomb <- as.numeric(args[6])
# num_gametes <- as.integer(args[7])
# num_snps <- as.integer(args[8])
# coverage <- as.numeric(args[9])
# message(paste0("The coverage of this simulation is: ", coverage))
# num_genotypes <- num_gametes * num_snps
# missing_genotype_rate <- dpois(0, coverage)
# message(paste0("The missing genotype rate of this simulation is ", missing_genotype_rate))
# num_nas <- floor(num_genotypes * missing_genotype_rate)
# window_length <- as.integer(args[10])
# overlap_denom <- as.integer(args[11])
# random_seed <- as.integer(args[12])
# set.seed(random_seed)
# recomb_lambda <- as.numeric(args[13])
# stopifnot(recomb_lambda > 0)
# num_recomb_sites <- rpois(num_gametes, recomb_lambda)
# message(paste0("Total number of recombination spots across gametes: ", sum(num_recomb_sites)))
# add_seq_error <- as.logical(args[14])
# seqError_add <- as.numeric(args[15])
# add_de_novo_mut <- as.logical(args[16])
# de_novo_lambda <- as.integer(args[17])
# de_novo_alpha <- as.numeric(args[18])
# de_novo_beta <- as.numeric(args[19])
# stopifnot(de_novo_beta > 0)
# smooth_imputed_genotypes <- as.logical(args[20])
# smooth_crossovers <- as.logical(args[21])
# write_out_plot <- as.logical(args[22])
# ###

####Alternative means to input coverage
##missing_genotype_rate <- as.numeric(args[8])
##if (missing_genotype_rate > 1){ #if entered as a percentage
##missing_genotype_rate <- missing_genotype_rate / 100
##}
##coverage <- -log(missing_genotype_rate)
####

###Genenrate 2 parental chromosomes of heterozygous sites
hap1 <- data.frame(V1 = sample(c(0, 1), size = num_snps, replace = TRUE)) #simulate first parental chromosome
hap2 <- abs(1-hap1) #switch bits to construct second parental
parental_haps <- data.frame(cbind(hap1, hap2))
colnames(parental_haps) <- c("Parental1", "Parental2")

###Generate gametes from these 2 parental chromosomes
generate_gam <- function(parental_haplotypes, n_crossovers){
  #this function generates a single gamete given parental haplotypes and the number of crossovers that should occurs for a given gamete
  # -- pick which parental haplotype corresponds to the beginning of the gamete's chromosome
  # -- assign the 0 & 1 genotypes of this first parental haplotype to the whole gamete
  # -- assign the 1 or 2 of reprenesting which haplotype to the whole gamete
  #     -- re-code this to 0 or 1 (just subtract 1)
  #     -- at the end re-recode this back to 1's and 2's (just add 1)
  # -- for each crossover index, i, between index i-1 and index i, switch the bits for the rest of the chromosome genotypes to the opposite parental haplotype  or genotype compared to the one just previous (abs(value-1))
  #input  -- parental_haplotypes: a dataframe with the two parental haplotypes
  #          n_crossovers: integer, the number of crossovers for this given gamete
  #output -- a list with 3 elements
  #             --crossover_indices, a vector of integers, such that each integer, wlog i, represents an index location of a crossover exchange point, such that the crossover occurred somewhere between SNP index locations i-1 and i
  #             --recombined_hap, a vector of 0's and 1's, representing the gamete's fully known generated genotype
  #             --recombined_hap_index, a vector of 1's and 2's, representing which parental haplotype that location's genotype originated from
  init_hap_index <- sample(1:2, 1) #which parental haplotype corresponds to the beginning of the gamete's chromosome
  init_hap <- parental_haplotypes[,init_hap_index]
  if (n_crossovers == 0) {
    return(list(NA, init_hap, rep(init_hap_index, length(init_hap))))
  } else {
    n_snps <- length(init_hap)
    crossover_indices <- sample(2:(n_snps-1), n_crossovers) #limiting it so that crossover locations can't be at the beginning or end to avoid the edge effect
    crossover_indices <- crossover_indices[order(crossover_indices)]
    recombined_hap <- init_hap
    recombined_hap_index <- rep(init_hap_index-1, n_snps) #recode 1's to 0's or 2's to 1's and store the initial haplotype at all SNP locations
    for (crossover in crossover_indices) { #switch haplotypes at some (unknown) location/exchange point between the indices crossover-1 and crossover
      recombined_hap <- c(recombined_hap[1:(crossover-1)], abs(1- recombined_hap[crossover:n_snps]))
      recombined_hap_index <- c(recombined_hap_index[1:(crossover-1)], abs(recombined_hap_index[crossover:n_snps]-1))
    }
    return(list(crossover_indices, recombined_hap, recombined_hap_index+1)) #re-recode recombined_hap_index so that 1's are 2's and 0' are 1's
  }
}

sim_gam <- lapply(1:num_gametes, function(x) generate_gam(parental_haps, num_recomb_sites[x]))
gam_mat <- sapply(sim_gam, "[[", 2)
gam_haps <- sapply(sim_gam, "[[", 3)
crossover_indices <- sapply(sim_gam, "[[", 1)

###Sparsify the knowledge of the gametes
add_to_na_flatten <- function(to_add_from, num_nas, num_gametes, num_snps){
  to_return <- rep(NA, (num_snps*num_gametes))
  coords_to_keep_genotype <- sample(1:(num_snps*num_gametes), size=((num_snps*num_gametes)-num_nas), replace = FALSE)
  to_return[coords_to_keep_genotype] <- as.vector(to_add_from[coords_to_keep_genotype])
  to_return <- matrix(to_return, ncol=num_gametes, nrow=num_snps)
  return (to_return)
}

gam_mat_with_na <- add_to_na_flatten(gam_mat, num_nas, num_gametes, num_snps)

###Include de novo mutations
if (add_de_novo_mut){
  new_rows <- c()
  num_dnm <- rpois(1, de_novo_lambda) + 1 #make sure this is greater than 0 by adding one
  message(paste0("Number of de novo mutations: ", num_dnm))
  num_gametes_affected_per_dnm <- ceiling(rgamma(num_dnm, de_novo_alpha, scale=de_novo_beta))
  parentals_with_dnm <- sample(1:2, num_dnm, replace=TRUE)
  for (i in 1:num_dnm){ #for every parental dnm
    parental_with_dnm <- parentals_with_dnm[i]
    row_loc <- sample(1:num_snps, 1) #pick random row which we'll add this new de novo mutation after
    gam_can_be_affected <- which(gam_haps[row_loc, ]==parental_with_dnm)
    num_gametes_affected <- min(num_gametes_affected_per_dnm[i], length(gam_can_be_affected))
    message(paste0("Number of gametes affected for de novo mutation ", i, ": ", num_gametes_affected))
    where_locs_greater <- which(new_rows >= (row_loc+1))
    new_rows[where_locs_greater] <- new_rows[where_locs_greater] + 1
    new_rows <- c(new_rows, (row_loc+1))
    message(paste0("new row location: ", row_loc+1))
    message(paste0("new rows vector: ", new_rows, collapse= " ; "))
    parental_haps <- rbind(parental_haps[1:row_loc,], rep(0, 2), parental_haps[(row_loc+1):num_snps, ])
    parental_haps[(row_loc+1), parental_with_dnm] <- 1
    affected_gam <- sort(sample(gam_can_be_affected, num_gametes_affected))
    #message(paste0("affected gametes: ", affected_gam, collapse= " ; "))
    new_row_haps <- rep(abs((parental_with_dnm - 1)-1)+1, num_gametes)
    new_row_haps[gam_can_be_affected] <- parental_with_dnm
    new_row_vals <- rep(0, num_gametes)
    new_row_vals[affected_gam] <- 1
    gam_mat <- rbind(gam_mat[1:row_loc, ], new_row_vals, gam_mat[(row_loc+1):num_snps, ])
    #add sparsity in
    num_nas_to_add <- floor(num_gametes * missing_genotype_rate)
    change_indices <- sample(1:num_gametes, num_nas_to_add)
    new_row_vals[change_indices] <- NA
    #add in new SNP line
    gam_mat_with_na <- rbind(gam_mat_with_na[1:row_loc, ], new_row_vals, gam_mat_with_na[(row_loc+1):num_snps, ])
    #add in new haps line
    gam_haps <- rbind(gam_haps[1:row_loc, ], new_row_haps, gam_haps[(row_loc+1):num_snps, ])
    #add in new indice and SNP
    #indices <- c(indices[1:row_loc], row_loc+1, indices[(row_loc+1):num_snps]+1)
    num_snps <- num_snps + 1 #increase num_snps by one
    #do we need to adjust recombination locations at all? Or are they fine being unchanged?
  }
}

td_test_truth <- function(gam_matrix, row_index) {
  test_row <- gam_matrix[row_index,]
  gt_vector <- unlist(test_row)[-1]
  one_count <- sum(gt_vector == 1, na.rm = TRUE)
  two_count <- sum(gt_vector == 2, na.rm = TRUE)
  p_value <- binom.test(c(one_count, two_count))$p.value
  return(c(p_value, one_count, two_count))
}

df_counts_pvals_truth <- do.call(rbind, pbmclapply(1:nrow(gam_haps),
                                                   function(x) td_test_truth(gam_haps, x),
                                                   mc.cores=getOption("mc.cores", threads))) %>%
  as_tibble() %>%
  add_column(1:nrow(gam_haps)) #bind the positions vector to df_counts_pvals
colnames(df_counts_pvals_truth) <- c("pval", "h1_count", "h2_count", "genomic_position")
if (write_out_plot){
  filename_df <- paste0(outDir, sampleName, "_", chrom, "_pval_sim_truth_full.csv")
  write_csv(df_counts_pvals_truth, filename_df)
}

###Include sequencing error  
if (add_seq_error){
  num_genotypes <- num_snps * num_gametes #updating since this may have changed while adding de novo mutations
  num_bits_to_flip <- as.integer(seqError_add * (num_genotypes - num_nas))
  switched_bit_mat <- c(abs(1-gam_mat_with_na)) #make a switched bit matrix compared to gam_mat_with_na
  where_locs <- sample(which(!is.na(switched_bit_mat)), size=num_bits_to_flip) #randomly pick num_bits_to_flip locations
  gam_mat_with_na <- c(gam_mat_with_na)
  gam_mat_with_na[where_locs] <- switched_bit_mat[where_locs] #take those locations from switched bit matrix and put them in place in the gam_mat_with_na matrix
  gam_mat_with_na <- matrix(gam_mat_with_na, nrow=num_snps, ncol=num_gametes)
}

##Transform into dataframes that can be passed to the rest of the pipeline -- genomic positions first, column names for each gamete
gam_na_df <- data.frame(pseudo_pos = 1:nrow(gam_mat_with_na), gam_mat_with_na)
colnames(gam_na_df) <- c("positions", paste0("gam", 1:num_gametes, "_"))
gam_full_df <- data.frame(pseudo_pos = 1:nrow(gam_mat), gam_mat)
colnames(gam_full_df) <- c("positions", paste0("gam", 1:num_gametes, "_"))

###Run the whole pipeline
##Ensure that each row has at least one 0 and one 1
keep_bool <- unname((rowSums(gam_na_df[, 2:ncol(gam_na_df)] == 0, na.rm = TRUE) > 0) & (rowSums(gam_na_df[, 2:ncol(gam_na_df)] == 1, na.rm = TRUE) > 0))
gam_na_df <- gam_na_df[keep_bool,]
gam_full_df <- gam_full_df[keep_bool,]
parental_haps <- parental_haps[keep_bool,]
num_snps <- sum(keep_bool)
message(paste0("new number of snps: ", num_snps))
if (add_de_novo_mut) {
  `%notin%` <- Negate(`%in%`)
  for (i in 1:length(new_rows)){
    message(paste0("dnm ", i, " is filtered out: ", new_rows[i] %notin% gam_na_df[,1]))
  }
}

if (write_out_plot){
  filename_df <- paste0(outDir, sampleName, "_", chrom, "_pval_sim_truth_filtered.csv")
  write_csv(df_counts_pvals_truth[keep_bool,], filename_df)
}

##Remove the first column (positions)
positions <- gam_na_df[, 1]
gam_na_df <- gam_na_df[,-1]

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
windows <- splitWithOverlap(rank(positions), window_length, overlap = window_length %/% overlap_denom)


# merge the last two windows to avoid edge effect
if (length(windows) > 1){
  combined <- unique(c(windows[[length(windows) - 1]], windows[[length(windows)]]))
  combined <- combined[order(combined)]
  total_combined <- windows[-c((length(windows) - 1), length(windows))]
  total_combined[[length(total_combined) + 1]] <- combined
  windows <- total_combined
}

message(paste0("Number of windows with overlap of ", window_length %/% overlap_denom , " and ", num_snps, " number of SNPs following filtering: ", length(windows)))


# function to reconstruct parental haplotypes
reconstruct_hap <- function(input_dt, input_positions, window_indices) {
  window_start <- min(window_indices)
  window_end <- max(window_indices)
  positions_for_window <- input_positions[window_start:window_end]
  # compute a distance matrix
  d <- dist(t(as.matrix(input_dt)[window_start:window_end,]), method = "binary")
  # put in 0.5 for any NA entries of the distance matrix
  d[is.na(d)] <- 0.5
  # cluster the distance matrix
  tree <- hclust(d, method = "ward.D2")
  # plot(tree, cex = 0.1) # uncomment to plot
  # cut the tree generated by clustering into two groups (haplotypes)
  haplotypes <- cutree(tree, k=2)
  # get the names of the gam cells falling into the two groups
  h1_gam <- names(haplotypes[haplotypes == 1])
  h2_gam <- names(haplotypes[haplotypes == 2])
  # reconstruct the original haplotypes by majority vote after inverting the opposite haplotype
  h1_inferred <- unname(apply(cbind(input_dt[window_start:window_end, h1_gam],
                                    invertBits(input_dt[window_start:window_end, h2_gam])),
                              1, function(x) getmode(x)))
  # h2_inferred <- unname(apply(cbind(input_dt[window_start:window_end, h2_gam],
  #                                   invertBits(input_dt[window_start:window_end, h1_gam])),
  #                             1, function(x) getmode(x)))
  return(tibble(index = window_indices, pos = positions_for_window, h1 = h1_inferred))
}

# infer the haplotypes within the overlapping windows
inferred_haplotypes <- pbmclapply(1:length(windows),
                                  function(x) reconstruct_hap(gam_na_df, positions, windows[[x]]),
                                  mc.cores = getOption("mc.cores", threads))

# stitch together the haplotypes
initial_haplotype <- inferred_haplotypes[[1]]
for (hap_window in 1:length(windows)) {
  olap_haps <- merge(initial_haplotype, inferred_haplotypes[[hap_window]], by = "index")
  olap_haps_complete <- merge(initial_haplotype, inferred_haplotypes[[hap_window]], by = "index", all = TRUE)
  mean_concordance <- mean(olap_haps$h1.x == olap_haps$h1.y, na.rm=TRUE)
  message(paste0("Window ", hap_window, "Mean Concordance: ", mean_concordance))
  if (mean_concordance < 0.1) {
    olap_haps_complete$h1.y <- invertBits(olap_haps_complete$h1.y)
  } else if (mean_concordance < 0.9) {
    message(paste0("Haplotypes within overlapping windows are too discordant to merge. Mean: ", mean_concordance, " Window: ", hap_window, "But continuing"))
    #stop(paste0("Haplotypes within overlapping windows are too discordant to merge. Mean: ", mean_concordance, " Window: ", hap_window))
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

hamming_distance_ignoreNA <- function(truth, predicted, num_snps){
  if (is.null(dim(truth))){
    num_mismatch <- sum((predicted - truth) != 0, na.rm=TRUE)
    to_return <- num_mismatch / num_snps * 100
  } else {
    num_mismatch_byCol <- colSums((predicted - truth) != 0, na.rm=TRUE)
    to_return <- num_mismatch_byCol / num_snps * 100
  }
  return (to_return)
}

completeness <- function(predicted, num_snps){
  if (is.null(dim(predicted))){
    num_nas <- sum(is.na(predicted))
    to_return <- 1 - (num_nas/num_snps)
  } else {
    num_nas_byCol <- colSums(is.na(predicted))
    to_return <- 1 - (num_nas_byCol/num_snps)
  }
  return (to_return)
}

lhs <- function(truth, predicted){ #need to figure out two dimensional case, call with apply?
  match01_encoding <- predicted - truth #match is a 0, mismatch is a 1 | -1
  rleres <- rle(match01_encoding) #any match location will be 0
  to_return <- max(rleres$lengths[which(rleres$values == 0)]) #find longest length of matches
  return (to_return)
}

ser <- function(truth, predicted, num_snps){ #also call with apply?
  match01_encoding <- predicted - truth #match is a 0, mismatch is a 1 | -1
  which_mismatch <- which(match01_encoding != 0)
  comp_with_before <- abs(diff(match01_encoding))[which_mismatch - 1]
  switch_errors <- sum(comp_with_before == 1, na.rm=TRUE) #if comp_with_before values are 0 or 2, then the value is a continuation following another error; if 1, it's a new switch error
  to_return <- switch_errors / num_snps
  return (to_return)
}

###Assessing the accuracy of parental haplotype reconstruction
accuracy_par_igna <- 100 - min(hamming_distance_ignoreNA(parental_haps$Parental1, complete_haplotypes$h1, num_snps),
                               hamming_distance_ignoreNA(parental_haps$Parental2, complete_haplotypes$h1, num_snps))

completeness_par <- completeness(complete_haplotypes$h1, num_snps)

lhs_par <- max(lhs(parental_haps$Parental1, complete_haplotypes$h1),
               lhs(parental_haps$Parental2, complete_haplotypes$h1))

ser_par <- min(ser(parental_haps$Parental1, complete_haplotypes$h1, num_snps),
               ser(parental_haps$Parental2, complete_haplotypes$h1, num_snps))


message(paste0("Phasing accuracy: ", accuracy_par_igna))
message(paste0("Phasing completeness: ", completeness_par))
message(paste0("Phasing LHS: ", lhs_par))
message(paste0("Phasing SER: ", ser_par))


# Going through each gamete, if an allele (0 or 1) in a gamete matches the allele (0 or 1)
# in h1 at that position, replace the allele with "h1". Do the same for h2.
for (i in 1:ncol(gam_na_df)) {
  gam_na_df[i][gam_na_df[i] == complete_haplotypes$h1] <- "h1"
  gam_na_df[i][gam_na_df[i] == complete_haplotypes$h2] <- "h2"
  gam_na_df[c(which(gam_na_df[,i] == 0 | gam_na_df[,i] == 1)),i] <- NA
}
####Might need to handle the 0's and 1's still in gam_na_df that don't match complete_haplotypes$h1 or complete_haplotypes$h2 because of NAs in complete_haplotypes

# Scan gamete by gamete to interpret state given emission
# First, we initialize our HMM
# set denominator for transition probability - one recombination event per chromosome
#num_snps <- nrow(complete_haplotypes)
# two states
states <- c("haplotype1", "haplotype2")
# probability of state at position x+1 given state at position x
hap1Prob <- c(1-(hmm_avg_recomb/num_snps), hmm_avg_recomb/num_snps)
hap2Prob <- c(hmm_avg_recomb/num_snps, 1-(hmm_avg_recomb/num_snps))
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
# Compute the inferred state using each gamete cell as the input
# (Input must be a vector)
runHMM <- function(gam_dt, column_index) {
  original_obs <- gam_dt[,column_index]
  inferred_state <- viterbi(hmm, na.omit(gam_dt[, column_index]))
  original_obs[!is.na(original_obs)] <- inferred_state
  return(original_obs)
}

imputed_gam <- as_tibble(do.call(cbind, pbmclapply(1:ncol(gam_na_df),
                                                   function(x) runHMM(gam_na_df, x),
                                                   mc.cores = getOption("mc.cores", threads))))
colnames(imputed_gam) <- colnames(gam_na_df)

# fill up and down at the end
fill_NAs <- function(merged_gam, col_index) {
  gam_sample <- merged_gam[,col_index] %>%
    rename(gam = colnames(.)[1]) %>%
    mutate(gam_up = gam) %>%
    mutate(gam_down = gam) %>%
    fill(gam_up, .direction = "up") %>%
    fill(gam_down, .direction = "down") %>%
    mutate(is_match = (gam_up == gam_down)) %>%
    replace_na(list(is_match = FALSE))
  gam_sample$gam_imputed <- as.character(NA)
  gam_sample[gam_sample$is_match == TRUE,]$gam_imputed <- gam_sample[gam_sample$is_match == TRUE,]$gam_up
  #fill beginning of chromosome NA's
  first <- which(!is.na(gam_sample$gam_imputed))[1]
  gam_sample$gam_imputed[1:(first-1)] <- gam_sample$gam_imputed[first]
  #fill end of chromosome NA's
  gam_sample$gam_imputed <- rev(gam_sample$gam_imputed)
  first <- which(!is.na(gam_sample$gam_imputed))[1]
  gam_sample$gam_imputed[1:(first-1)] <- gam_sample$gam_imputed[first]
  #reverse chromosome imputation back so it faces the right way
  gam_sample$gam_imputed <- rev(gam_sample$gam_imputed)
  return(gam_sample$gam_imputed)
}

filled_gam <- as_tibble(do.call(cbind,
                                pblapply(1:ncol(imputed_gam),
                                         function(x) fill_NAs(imputed_gam, x))))
colnames(filled_gam) <- colnames(gam_na_df)

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
  filled_gam_recomb <- unsmooth(gam_na_df, filled_gam)
  idents_for_csv <- paste0(paste0(sampleName, "_", chrom, "_"), colnames(filled_gam_recomb))
  recomb_spots_all <- do.call(rbind, pbmclapply(1:ncol(filled_gam_recomb),
                                                function(x) find_recomb_spots(filled_gam_recomb, x, idents_for_csv, positions),
                                                mc.cores=getOption("mc.cores", threads))) %>%
    right_join(., tibble(Ident = idents_for_csv), by = "Ident")
} else{
  idents_for_csv <- paste0(paste0(sampleName, "_", chrom, "_"), colnames(filled_gam))
  recomb_spots_all <- do.call(rbind, pbmclapply(1:ncol(filled_gam),
                                                function(x) find_recomb_spots(filled_gam, x, idents_for_csv, positions),
                                                mc.cores=getOption("mc.cores", threads))) %>%
    right_join(., tibble(Ident = idents_for_csv), by = "Ident")
}

###Assessing the accuracy of gamete haplotype reconstruction
re_recode_gametes <- function(dt, complete_haplotypes) {
  to_return <- data.frame(matrix(NA_real_, nrow=nrow(dt), ncol=ncol(dt)))
  for (i in 1:ncol(dt)) {
    locs_h1 <- dt[,i] == "haplotype1"
    locs_h1[which(is.na(locs_h1))] <- FALSE
    locs_h2 <- dt[,i] == "haplotype2"
    locs_h2[which(is.na(locs_h2))] <- FALSE
    to_return[locs_h1, i] <- complete_haplotypes$h1[locs_h1]
    to_return[locs_h2, i] <- complete_haplotypes$h2[locs_h2]
  }
  return(to_return)
}

if (!smooth_imputed_genotypes & smooth_crossovers){
  filled_gam <- unsmooth(gam_na_df, filled_gam)
  filled_gam_recode <- re_recode_gametes(filled_gam, complete_haplotypes)
} else if(!smooth_imputed_genotypes & !smooth_crossovers){
  filled_gam_recode <- re_recode_gametes(filled_gam_recomb, complete_haplotypes)
} else {
  filled_gam_recode <- re_recode_gametes(filled_gam, complete_haplotypes)
}

accuracy_gam <- 100-hamming_distance_ignoreNA(gam_full_df[,-1], filled_gam_recode, num_snps)
completeness_gam <- completeness(filled_gam_recode, num_snps)
lhs_gam <- do.call(rbind, lapply(1:num_gametes,function(x) lhs(gam_full_df[,-1][,x], filled_gam_recode[,x])))
ser_gam <- do.call(rbind, lapply(1:num_gametes,function(x) ser(gam_full_df[,-1][,x], filled_gam_recode[,x], num_snps)))

message(paste0("Mean gamete haplotype reconstruction accuracy: ", mean(accuracy_gam, na.rm=TRUE)))
message(paste0("Stdev gamete haplotype reconstruction accuracy: ", sd(accuracy_gam, na.rm=TRUE)))
message(paste0("Mean gamete completeness: ", mean(completeness_gam, na.rm=TRUE)))
message(paste0("Stdev gamete completeness: ", sd(completeness_gam, na.rm=TRUE)))
message(paste0("Mean gamete LHS: ", mean(lhs_gam, na.rm=TRUE)))
message(paste0("Stdev gamete LHS: ", sd(lhs_gam, na.rm=TRUE)))
message(paste0("Mean gamete SER: ", mean(ser_gam, na.rm=TRUE)))
message(paste0("Stdev gamete SER: ", sd(ser_gam, na.rm=TRUE)))


#### bedr --> foverlap

if (write_out_plot){
  filename = paste0(outDir, "simulated_gam_hap_reconstruction.pdf")
  pdf(file=filename)
  ggplot(accuracyDat, aes(x=name, y=val, fill=name)) + scale_x_discrete(limits=c("Raw", "Corrected")) + geom_violin(scale="width", adjust=1, width=0.5) + labs(y="accuracy", x="method", title="gam Haplotype Reconstruction") + scale_y_continuous(name="accuracy", breaks=c(98.75, 99, 99.25, 99.5, 99.75, 100), labels=c(98.75, 99, 99.25, 99.5, 99.75, 100), limits=c(98.75, 100))
  dev.off()
}

td_test <- function(gam_matrix, row_index) {
  test_row <- gam_matrix[row_index,]
  gt_vector <- unlist(test_row)[-1]
  one_count <- sum(gt_vector == "haplotype1", na.rm = TRUE)
  two_count <- sum(gt_vector == "haplotype2", na.rm = TRUE)
  p_value <- binom.test(c(one_count, two_count))$p.value
  return(c(p_value, one_count, two_count))
}

df_counts_pvals <- do.call(rbind, pbmclapply(1:nrow(filled_gam),
                                            function(x) td_test(filled_gam, x),
                                            mc.cores=getOption("mc.cores", threads))) %>%
 as_tibble() %>%
 add_column(positions) #bind the positions vector to df_counts_pvals
colnames(df_counts_pvals) <- c("pval", "h1_count", "h2_count", "genomic_position")
if (write_out_plot){
  filename_df <- paste0(outDir, sampleName, "_", chrom, "_pval_sim.csv")
  write_csv(df_counts_pvals, filename_df)
}

idents_for_csv <- paste0(paste0(sampleName, "_", chrom, "_"), colnames(filled_gam))
recomb_spots_all <- do.call(rbind, pbmclapply(1:ncol(filled_gam),
                                              function(x) find_recomb_spots(filled_gam, x, idents_for_csv, positions),
                                              mc.cores=getOption("mc.cores", threads))) %>%
  right_join(., tibble(Ident = idents_for_csv), by = "Ident")

metrics <- function(tp, fp, tn, fn){
    precision <- tp/(tp+fp)
    recall <- tp/(tp+fn)
    accuracy <- (tp + tn)/(tp + tn + fp + fn)
    f1 <- (2*precision*recall)/(precision + recall)
    specificity <- tn/(tn+fp)
    fdr <- fp/(tp+fp)
    fpr <- fp/(tn+fp)
    metric_list <- list(precision=precision,
                        recall=recall,
                        accuracy=accuracy,
                        specificity = specificity,
                        fdr = fdr,
                        fpr = fpr,
                        f1=f1)
    return (metric_list) }

find_tp <- function(truth_intersect_dt, cons=FALSE){
  tp <- nrow(truth_intersect_dt[!is.na(truth_intersect_dt$Predicted_Start),]) #want to count +1 for every truth that intersects a prediction, even if multiple truths intersect a single prediction; accomplish this by counting number of truths that intersect any prediction from the non NA intersection
  if (cons){
    tp <- tp - sum(duplicated(paste0(truth_intersect_dt[!is.na(truth_intersect_dt$Predicted_Start),]$gam, "_", truth_intersect_dt[!is.na(truth_intersect_dt$Predicted_Start),]$Predicted_Start, "_", truth_intersect_dt[!is.na(truth_intersect_dt$Predicted_Start),]$Predicted_End))) #want to count only one truth when multiple truths intersect a single prediction; accomplish this by subtracting the sum of the boolean reporting which items are duplicated in the prediction identities from the non NA intersection; note duplicated() only returns TRUE for the 2nd, 3rd, 4th, etc occurrences, not the 1st
  }
  return (tp)
}

find_tn <- function(truth_dt_na, pred_dt_na){
  tn <- nrow(merge(truth_dt_na, pred_dt_na, by="gam")) #both truth and predicted return NA for a given gamete meaning they both agree that gamete has no recombiantion events
  return (tn)
}

find_fp <- function(pred_intersect_dt){
  fp <- nrow(pred_intersect_dt[is.na(pred_intersect_dt$True_Start),]) #number predicted as recombination spots, but don't intersect with truth at all
  return (fp)
}

find_fn <- function(truth_intersect_dt, cons=FALSE){
  fn <- nrow(truth_intersect_dt[is.na(truth_intersect_dt$Predicted_Start),]) #number of true recombinantion spots that don't interesect any predictions
  if (cons){
    fn <- fn + sum(duplicated(paste0(truth_intersect_dt[!is.na(truth_intersect_dt$Predicted_Start),]$gam, "_", truth_intersect_dt[!is.na(truth_intersect_dt$Predicted_Start),]$Predicted_Start, "_", truth_intersect_dt[!is.na(truth_intersect_dt$Predicted_Start),]$Predicted_End))) #count the overflows when multiple truths intersect a single predicted
  }
  return (fn)
}

names(crossover_indices) <- paste0(rep("gam", num_gametes), 1:num_gametes, "_")
unlist_ci <- unlist(crossover_indices, use.names=TRUE)
tci_dt <- data.table(gam=str_split(names(unlist_ci), "_", simplify=TRUE)[,1], start =(unlist_ci-1), end=(unlist_ci))
tci_dt_nona <- tci_dt[!is.na(tci_dt$start),] %>% setkey()
tci_dt_na <- tci_dt[is.na(tci_dt$start),]

pci_dt <- data.table(gam=str_split(recomb_spots_all$Ident, "_", simplify=TRUE)[,3], start=recomb_spots_all$Genomic_start, end=recomb_spots_all$Genomic_end)
pci_dt_nona <- pci_dt[!is.na(pci_dt$start),] %>% setkey()
pci_dt_na <- pci_dt[is.na(pci_dt$start),]

if (nrow(tci_dt_nona) > 0 & nrow(pci_dt_nona) > 0){ #truths and predictions
  tci_intersect <- foverlaps(tci_dt_nona, pci_dt_nona) %>%  `colnames<-` (c("gam", "Predicted_Start", "Predicted_End", "True_Start", "True_End"))
  pci_intersect <- foverlaps(pci_dt_nona, tci_dt_nona) %>%  `colnames<-` (c("gam", "True_Start", "True_End", "Predicted_Start", "Predicted_End"))
  
  tp_lib <- find_tp(tci_intersect)
  tp_cons <- find_tp(tci_intersect, cons = TRUE)
  fn_lib <- find_fn(tci_intersect)
  fn_cons <- find_fn(tci_intersect, cons = TRUE)
  tn <- find_tn(tci_dt_na, pci_dt_na)
  fp <- find_fp(pci_intersect)
  
  metrics_cons <- metrics(tp_cons, fp, tn, fn_cons)
  message(paste0("Conservative recombination spot identification metrics\nPrecision: ", metrics_cons$precision,
                 "\nRecall: ", metrics_cons$recall,
                 "\nAccuracy: ", metrics_cons$accuracy,
                 "\nF1: ", metrics_cons$f1,
                 "\nSpecificity: ", metrics_cons$specificity,
                 "\nFDR: ", metrics_cons$fdr,
                 "\nFPR: ", metrics_cons$fpr))
  metrics_lib <- metrics(tp_lib, fp, tn, fn_lib)
  message(paste0("Liberal recombination spot identification metrics\nPrecision: ", metrics_lib$precision,
                 "\nRecall: ", metrics_lib$recall,
                 "\nAccuracy: ", metrics_lib$accuracy,
                 "\nF1: ", metrics_lib$f1,
                 "\nSpecificity: ", metrics_lib$specificity,
                 "\nFDR: ", metrics_lib$fdr,
                 "\nFPR: ", metrics_lib$fpr))
  
} else if (nrow(tci_dt_nona) == 0 & nrow(pci_dt_nona) > 0){ #no truths, but some predictions; i.e. no non-na truths, but there are non-na predictions, all are false positives
  fp <- nrow(pci_dt_nona)
  fn <- 0
  tp <- 0
  tn <- nrow(merge(tci_dt_na, pci_dt_na, by="gam")) #predicted and truth both return no crossovers for a given gamete
  message(paste0("recombination metrics match for liberal and conservative because no truths"))
  metrics_it <- metrics(tp, fp, tn, fn)
  message(paste0("Conservative recombination spot identification metrics\nPrecision: ", metrics_it$precision,
                 "\nRecall: ", NA, #metrics_it$recall,
                 "\nAccuracy: ", metrics_it$accuracy,
                 "\nF1: ", NA, #metrics_it$f1,
                 "\nSpecificity: ", metrics_it$specificity,
                 "\nFDR: ", metrics_it$fdr,
                 "\nFPR: ", metrics_it$fpr))
  message(paste0("Liberal recombination spot identification metrics\nPrecision: ", metrics_it$precision,
                 "\nRecall: ", NA, #metrics_it$recall,
                 "\nAccuracy: ", metrics_it$accuracy,
                 "\nF1: ", NA, #metrics_it$f1,
                 "\nSpecificity: ", metrics_it$specificity,
                 "\nFDR: ", metrics_it$fdr,
                 "\nFPR: ", metrics_it$fpr))
} else if (nrow(pci_dt_nona) == 0 & nrow(tci_dt_nona) >0){ #no predictions, but some truths, i.e. no non-na predictions but there are non-na truths, all are false negatives
  fn <- nrow(tci_dt_nona)
  fp <- 0
  tn <- nrow(merge(tci_dt_na, pci_dt_na, by="gam")) #predicted and truth both return no crossovers for a given gamete
  tp <- 0
  message(paste0("recombination metrics match for liberal and conservative because no predictions"))
  metrics_it <- metrics(tp, fp, tn, fn)
  message(paste0("Conservative recombination spot identification metrics\nPrecision: ", NA, #metrics_it$precision,
                 "\nRecall: ", metrics_it$recall,
                 "\nAccuracy: ", metrics_it$accuracy,
                 "\nF1: ", NA, #metrics_it$f1,
                 "\nSpecificity: ", metrics_it$specificity,
                 "\nFDR: ", NA, #metrics_it$fdr,
                 "\nFPR: ", metrics_it$fpr))
  message(paste0("Liberal recombination spot identification metrics\nPrecision: ", NA, #metrics_it$precision,
                 "\nRecall: ", metrics_it$recall,
                 "\nAccuracy: ", metrics_it$accuracy,
                 "\nF1: ", NA, #metrics_it$f1,
                 "\nSpecificity: ", metrics_it$specificity,
                 "\nFDR: ", NA, #metrics_it$fdr,
                  "\nFPR: ", metrics_it$fpr))
} else {
  message("no recombination metrics possible because there were no truths or predictions")
  message(paste0("Conservative recombination spot identification metrics\nPrecision: ", NA,
                                  "\nRecall: ", NA,
                                  "\nAccuracy: ", NA,
                                  "\nF1: ", NA,
                                  "\nSpecificity: ", NA,
                                  "\nFDR: ", NA,
                                  "\nFPR: ", NA))
                   message(paste0("Liberal recombination spot identification metrics\nPrecision: ", NA,
                                  "\nRecall: ", NA,
                                  "\nAccuracy: ", NA,
                                  "\nF1: ", NA,
                                  "\nSpecificity: ", NA,
                                  "\nFDR: ", NA,
                                  "\nFPR: ", NA))
}

# if (write_out_plot){
#   ###Plot the resolution of predicted recombination break point regions
#   filename = paste0(outDir, "simulated_recombination_resolution.pdf")
#   pdf(file=filename)
#   hist(recomb_spots_df$end - recomb_spots_df$start, xlab = "Differince in SNP index", breaks=100, main="Resolution of predicted recombination break point regions")
#   dev.off()
#   
#   ###Plot the distribution of sparsity by row/SNP in the simulated gamete data
#   real_reads <- rowSums(!is.na(gam_na_df))
#   filename = paste0(outDir, "simulated_notna_bysnp.pdf")
#   pdf(file=filename)
#   hist(real_reads, breaks=150, xlab="not NA by SNP", main="Distribution of simulated data not NA by SNP")
#   dev.off()
#   
#   ###Plot the distribution of sparsity by column/gamete in the simulated gamete data
#   real_reads2 <- colSums(!is.na(gam_na_df))
#   filename = paste0(outDir, "simulated_notna_bygam.pdf")
#   pdf(file=filename)
#   hist(real_reads2, breaks=100, xlab="not NA by gam", main="Distribution of simulated data not NA by gam")
#   dev.off()
#   
#   ###Plot where the de novo mutations are with respect to the mismatches in parental haplotypes and inference of parental haplotypes
#   if (add_de_novo_mut){
#     mismatches_parental <- which((complete_haplotypes$h1 -  parental_haps$Parental1)  != 0)
#     filename = paste0(outDir, "simulated_parental_mismatch_loc.pdf")
#     pdf(file=filename)
#     plot(c(mismatches_parental, new_rows), c(rep(1, length(mismatches_parental)), rep(1.05, length(new_rows))), cex=0.2, col=c(rep(rgb(red=0, green=0, blue=0, alpha=0.2), length(mismatches_parental)), rep(rgb(red=1, green=0, blue=0, alpha=0.8), length(new_rows))))
#     dev.off()
#   }
# }