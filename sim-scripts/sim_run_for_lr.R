library(data.table)
library(tidyverse)
library(pbapply)
library(pbmcapply)

args <- commandArgs(trailingOnly = TRUE)
sampleName <- "simLR"
chrom <- "chrT"
threads <- as.integer(args[1])

window_length <- as.integer(args[2])
overlap_denom <- as.numeric(args[3])
seqError <- 0.05
hapProb <- 1 - seqError

num_sperm <- as.integer(args[4])
num_snps <-  as.integer(args[5])
coverage <- as.numeric(args[6])

num_genotypes <- num_sperm * num_snps
missing_genotype_rate <- dpois(0, coverage)
message(paste0("The missing genotype rate of this simulation is ", missing_genotype_rate))
num_nas <- floor(num_genotypes * missing_genotype_rate)

random_seed <- as.integer(args[7])
set.seed(random_seed)

recomb_lambda <- as.numeric(args[8])
num_recomb_sites <- rpois(num_sperm, recomb_lambda)
message(paste0("Total number of recombination spots across gametes: ", sum(num_recomb_sites)))

add_seq_error <- TRUE
seqError_add <- 0.05
add_de_novo_mut <- TRUE
de_novo_lambda <- 5
de_novo_alpha <- 7.5
de_novo_beta <- 10
stopifnot(de_novo_beta > 0)

smooth_imputed_genotypes <- FALSE
smooth_crossovers <- TRUE

###Genenrate 2 parental chromosomes of heterozygous sites
hap1 <- data.frame(V1 = sample(c(0, 1), size = num_snps, replace = TRUE)) #simulate first parental chromosome
hap2 <- abs(1-hap1) #switch bits to construct second parental
parental_haps <- data.frame(cbind(hap1, hap2))
colnames(parental_haps) <- c("Parental1", "Parental2")

###Generate gametes from these 2 parental chromosomes
generate_sperm <- function(parental_haplotypes, n_crossovers){
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

sim_sperm <- lapply(1:num_sperm, function(x) generate_sperm(parental_haps, num_recomb_sites[x]))
sperm_mat <- sapply(sim_sperm, "[[", 2)
sperm_haps <- sapply(sim_sperm, "[[", 3)
crossover_indices <- sapply(sim_sperm, "[[", 1)

###Sparsify the knowledge of the gametes
add_to_na_flatten <- function(to_add_from, num_nas, num_sperm, num_snps){
  to_return <- rep(NA, (num_snps*num_sperm))
  coords_to_keep_genotype <- sample(1:(num_snps*num_sperm), size=((num_snps*num_sperm)-num_nas), replace = FALSE)
  to_return[coords_to_keep_genotype] <- as.vector(to_add_from[coords_to_keep_genotype])
  to_return <- matrix(to_return, ncol=num_sperm, nrow=num_snps)
  return (to_return)
}

sperm_mat_with_na <- add_to_na_flatten(sperm_mat, num_nas, num_sperm, num_snps)

###Include de novo mutations
if (add_de_novo_mut){
  new_rows <- c()
  num_dnm <- rpois(1, de_novo_lambda) + 1 #make sure this is greater than 0 by adding one
  message(paste0("Number of de novo mutations: ", num_dnm))
  num_sperm_affected_per_dnm <- ceiling(rgamma(num_dnm, de_novo_alpha, scale=de_novo_beta))
  parentals_with_dnm <- sample(1:2, num_dnm, replace=TRUE)
  for (i in 1:num_dnm){ #for every parental dnm
    parental_with_dnm <- parentals_with_dnm[i]
    row_loc <- sample(1:num_snps, 1) #pick random row which we'll add this new de novo mutation after
    sperm_can_be_affected <- which(sperm_haps[row_loc, ]==parental_with_dnm)
    num_sperm_affected <- min(num_sperm_affected_per_dnm[i], length(sperm_can_be_affected))
    message(paste0("Number of sperm affected for de novo mutation ", i, ": ", num_sperm_affected))
    where_locs_greater <- which(new_rows >= (row_loc+1))
    new_rows[where_locs_greater] <- new_rows[where_locs_greater] + 1
    new_rows <- c(new_rows, (row_loc+1))
    message(paste0("new row location: ", row_loc+1))
    message(paste0("new rows vector: ", new_rows, collapse= " ; "))
    parental_haps <- rbind(parental_haps[1:row_loc,], rep(0, 2), parental_haps[(row_loc+1):num_snps, ])
    parental_haps[(row_loc+1), parental_with_dnm] <- 1
    affected_sperm <- sort(sample(sperm_can_be_affected, num_sperm_affected))
    #message(paste0("affected sperm: ", affected_sperm, collapse= " ; "))
    new_row_haps <- rep(abs((parental_with_dnm - 1)-1)+1, num_sperm)
    new_row_haps[sperm_can_be_affected] <- parental_with_dnm
    new_row_vals <- rep(0, num_sperm)
    new_row_vals[affected_sperm] <- 1
    sperm_mat <- rbind(sperm_mat[1:row_loc, ], new_row_vals, sperm_mat[(row_loc+1):num_snps, ])
    #add sparsity in
    num_nas_to_add <- floor(num_sperm * missing_genotype_rate)
    change_indices <- sample(1:num_sperm, num_nas_to_add)
    new_row_vals[change_indices] <- NA
    #add in new SNP line
    sperm_mat_with_na <- rbind(sperm_mat_with_na[1:row_loc, ], new_row_vals, sperm_mat_with_na[(row_loc+1):num_snps, ])
    #add in new haps line
    sperm_haps <- rbind(sperm_haps[1:row_loc, ], new_row_haps, sperm_haps[(row_loc+1):num_snps, ])
    #add in new indice and SNP
    #indices <- c(indices[1:row_loc], row_loc+1, indices[(row_loc+1):num_snps]+1)
    num_snps <- num_snps + 1 #increase num_snps by one
    #do we need to adjust recombination locations at all? Or are they fine being unchanged?
  }
}

###Include sequencing error
if (add_seq_error){
  num_genotypes <- num_snps * num_sperm #updating since this may have changed while adding de novo mutations
  num_bits_to_flip <- as.integer(seqError_add * (num_genotypes - num_nas))
  switched_bit_mat <- c(abs(1-sperm_mat_with_na)) #make a switched bit matrix compared to sperm_mat_with_na
  where_locs <- sample(which(!is.na(switched_bit_mat)), size=num_bits_to_flip) #randomly pick num_bits_to_flip locations
  sperm_mat_with_na <- c(sperm_mat_with_na)
  sperm_mat_with_na[where_locs] <- switched_bit_mat[where_locs] #take those locations from switched bit matrix and put them in place in the sperm_mat_with_na matrix
  sperm_mat_with_na <- matrix(sperm_mat_with_na, nrow=num_snps, ncol=num_sperm)
}

##Transform into dataframes that can be passed to the rest of the pipeline -- genomic positions first, column names for each sperm
sperm_na_df <- data.frame(pseudo_pos = 1:nrow(sperm_mat_with_na), sperm_mat_with_na)
colnames(sperm_na_df) <- c("positions", paste0("sperm", 1:num_sperm, "_"))
sperm_full_df <- data.frame(pseudo_pos = 1:nrow(sperm_mat), sperm_mat)
colnames(sperm_full_df) <- c("positions", paste0("sperm", 1:num_sperm, "_"))

###Run the whole pipeline
##Ensure that each row has at least one 0 and one 1
keep_bool <- unname((rowSums(sperm_na_df[, 2:ncol(sperm_na_df)] == 0, na.rm = TRUE) > 0) & (rowSums(sperm_na_df[, 2:ncol(sperm_na_df)] == 1, na.rm = TRUE) > 0))
sperm_na_df <- sperm_na_df[keep_bool,]
sperm_full_df <- sperm_full_df[keep_bool,]
parental_haps <- parental_haps[keep_bool,]
num_snps <- sum(keep_bool)
#if (num_snps == 0) {quit with error message}
message(paste0("new number of snps: ", num_snps))
if (add_de_novo_mut) {
  `%notin%` <- Negate(`%in%`)
  for (i in 1:length(new_rows)){
    message(paste0("dnm ", i, " is filtered out: ", new_rows[i] %notin% sperm_na_df[,1]))
  }
}

##Remove the first column (positions)
positions <- sperm_na_df[, 1]
sperm_na_df <- sperm_na_df[,-1]

#new function but NAs likely in complete_haplotypes
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
  # get the names of the sperm cells falling into the two groups
  h1_sperm <- names(haplotypes[haplotypes == 1])
  h2_sperm <- names(haplotypes[haplotypes == 2])
  # reconstruct the original haplotypes by majority vote after inverting the opposite haplotype
  h1_inferred <- unname(apply(cbind(input_dt[window_start:window_end, h1_sperm],
                                    invertBits(input_dt[window_start:window_end, h2_sperm])),
                              1, function(x) getmode(x)))
  #h2_inferred <- unname(apply(cbind(input_dt[window_start:window_end, h2_sperm],
  #                                  invertBits(input_dt[window_start:window_end, h1_sperm])),
  #                            1, function(x) getmode(x)))
  return(tibble(index = window_indices, pos = positions_for_window, h1 = h1_inferred))
}

# infer the haplotypes within the overlapping windows
inferred_haplotypes <- pbmclapply(1:length(windows),
                                  function(x) reconstruct_hap(sperm_na_df, positions, windows[[x]]),
                                  mc.cores = getOption("mc.cores", threads))

# stitch together the haplotypes
initial_haplotype <- inferred_haplotypes[[1]]
for (hap_window in 1:length(windows)) {
  olap_haps <- merge(initial_haplotype, inferred_haplotypes[[hap_window]], by = "index")
  olap_haps_complete <- merge(initial_haplotype, inferred_haplotypes[[hap_window]], by = "index", all = TRUE)
  #mean_concordance <- mean(olap_haps$h1.x == olap_haps$h1.y)
  mean_concordance2 <- mean(olap_haps$h1.x == olap_haps$h1.y, na.rm=TRUE)
  message(paste0("Window ", hap_window, "Mean Concordance: ", mean_concordance2))
  if (mean_concordance2 < 0.1) {
    olap_haps_complete$h1.y <- invertBits(olap_haps_complete$h1.y)
  } else if (mean_concordance2 < 0.9) {
    message(paste0("Haplotypes within overlapping windows are too discordant to merge. Mean: ", mean_concordance2, " Window: ", hap_window, "But continuing"))
    #stop(paste0("Haplotypes within overlapping windows are too discordant to merge. Mean: ", mean_concordance2, " Window: ", hap_window))
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

hamming_distance_ignoreNA <- function(truth, predicted, num_snps){ #to use
  num_mismatch <- sum((predicted - truth) != 0, na.rm=TRUE)
  to_return <- num_mismatch / num_snps * 100
  #to_return <- (num_snps - num_mismatch) / num_snps * 100 #this is accuracy
  return (to_return)
}

# hamming_distance_penalizeNA <- function(truth, predicted, num_snps){
#   num_mismatch <- sum((predicted - truth) != 0, na.rm=TRUE)
#   num_nas_dif <- sum(is.na(predicted - truth))
#   to_return <- (num_mismatch + num_nas_dif) / num_snps * 100
#   #to_return <- (num_snps - (num_mismatch + num_nas_dif))/ num_snps * 100 #this is accuracy
#   return (to_return)
# }
#
# hamming_distance_forgiveNA <- function(truth, predicted, num_snps){
#   num_na_pred <- sum(is.na(predicted))
#   num_mismatch <- sum((predicted - truth) != 0, na.rm = TRUE)
#   to_return <- num_mismatch / (num_snps - num_na_pred) * 100
#   #to_return <- ((num_snps - num_na_pred) - num_mismatch) / (num_snps - num_na_pred) * 100 #this is accuracy
#   return (to_return)
# }

completeness <- function(predicted, num_snps){
  num_nas <- sum(is.na(predicted))
  to_return <- 1 - (num_nas/num_snps)
  return (to_return)
}



###Assessing the accuracy of parental haplotype reconstruction
accuracy_par_igna <- 100 - min(hamming_distance_ignoreNA(parental_haps$Parental1, complete_haplotypes$h1, num_snps),
                              hamming_distance_ignoreNA(parental_haps$Parental2, complete_haplotypes$h1, num_snps))

completeness_par <- completeness(complete_haplotypes$h1, num_snps)

message(paste0("Phasing accuracy: ", accuracy_par_igna))
message(paste0("Phasing completeness: ", completeness_par))
