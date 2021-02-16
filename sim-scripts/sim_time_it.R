library(data.table)
library(tidyverse)
library(bedr)
library(stringr)
library(pbapply)
library(pbmcapply)
library(HMM)
library(ggplot2)


start_script_time <- Sys.time()
sampleName <- "simulation"
chrom <- "chrT"
outDir <- "short_time_it/"
threads <- 4

seqError <- 0.005
hapProb <- 1 - seqError
window_length <- 3000

smooth <- TRUE

num_sperm <- 5000
num_snps <- 10000 
coverage <- 0.01

message(paste0("The coverage of this simulation is: ", coverage))

random_seed <- 42
set.seed(random_seed)


recomb_lambda <- 1
num_recomb_sites <- rpois(num_sperm, recomb_lambda)
message(paste0("Total number of recombination spots across gametes: ", sum(num_recomb_sites)))

add_seq_error <- TRUE
add_de_novo_mut <- TRUE
seqError_add <- 0.005
de_novo_lambda <- 15
de_novo_alpha <- 12.5
de_novo_beta <- 20
stopifnot(de_novo_beta > 0)

#add to the base num_not_nan_per_row_base to find out the actual number of not NAs for each row. We'll now have a vector rather than a single integer
num_genotypes <- num_sperm * num_snps
missing_genotype_rate <- dpois(0, coverage)
num_nas <- floor(num_genotypes * missing_genotype_rate)

#start with 2 parental chromosomes of heterozygous sites
#first parental chromosome
before_time <- Sys.time()
hap1 <- data.frame(V1 = sample(c(0, 1), size = num_snps, replace = TRUE))
hap2 <- abs(1-hap1)
parental_haps <- data.frame(cbind(hap1, hap2))
colnames(parental_haps) <- c("Parental1", "Parental2")
after_time <- Sys.time()
time_to_impute <- difftime(after_time, before_time, units="mins")
message(paste0("Time used for generating parental haplotypes: ", time_to_impute))


#indices <- 1:num_snps #we'll use these as both indices and a pseudo genomic position
#this data structure only appears to be used here and when we record the locs of the dnm's here later. Get rid of?
before_time <- Sys.time()
generate_sperm <- function(parental_haplotypes, n_crossovers){
  init_hap_index <- sample(1:2, 1)
  init_hap <- parental_haplotypes[,init_hap_index]
  if (n_crossovers == 0) {
    return(list(NA, init_hap, rep(init_hap_index, length(init_hap))))
  } else {
    n_snps <- length(init_hap)
    crossover_indices <- sample(1:n_snps, n_crossovers)
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
after_time <- Sys.time()
time_to_impute <- difftime(after_time, before_time, units="mins")
message(paste0("Time used for generating gametes and recording needed info: ", time_to_impute))

before_time <- Sys.time()
add_to_na_flatten <- function(to_add_from, num_nas, num_sperm, num_snps){
  to_return <- rep(NA, (num_snps*num_sperm))
  coords_to_keep_genotype <- sample(1:(num_snps*num_sperm), size=((num_snps*num_sperm)-num_nas), replace = FALSE)
  to_return[coords_to_keep_genotype] <- as.vector(to_add_from[coords_to_keep_genotype])
  to_return <- matrix(to_return, ncol=num_sperm, nrow=num_snps)
  return (to_return)
}

sperm_mat_with_na <- add_to_na_flatten(sperm_mat, num_nas, num_sperm, num_snps)
after_time <- Sys.time()
time_to_impute <- difftime(after_time, before_time, units="mins")
message(paste0("Time used for sparsifying gametes: ", time_to_impute))

before_time <- Sys.time()
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
    message(paste0("affected sperm: ", affected_sperm, collapse= " ; "))
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
after_time <- Sys.time()
time_to_impute <- difftime(after_time, before_time, units="mins")
message(paste0("Time used for adding dnm's: ", time_to_impute))

before_time <- Sys.time()
td_test_truth <- function(sperm_matrix, row_index) {
  test_row <- sperm_matrix[row_index,]
  gt_vector <- unlist(test_row)[-1]
  one_count <- sum(gt_vector == 1, na.rm = TRUE)
  two_count <- sum(gt_vector == 2, na.rm = TRUE)
  p_value <- binom.test(c(one_count, two_count))$p.value
  return(c(p_value, one_count, two_count))
}

df_counts_pvals_truth <- do.call(rbind, pbmclapply(1:nrow(sperm_haps),
                                             function(x) td_test_truth(sperm_haps, x),
                                             mc.cores=getOption("mc.cores", threads))) %>%
  as_tibble() %>%
  add_column(1:nrow(sperm_haps)) #bind the positions vector to df_counts_pvals
colnames(df_counts_pvals_truth) <- c("pval", "h1_count", "h2_count", "genomic_position")
filename_df <- paste0(outDir, sampleName, "_", chrom, "_pval_sim_truth_full.csv")
write_csv(df_counts_pvals_truth, filename_df)
after_time <- Sys.time()
time_to_impute <- difftime(after_time, before_time, units="mins")
message(paste0("Time used for finding and writing TD in simualted gametes: ", time_to_impute))

before_time <- Sys.time()
if (add_seq_error){
  num_genotypes <- num_snps * num_sperm #updating since this may have changed while adding de novo mutations
  num_bits_to_flip <- as.integer(seqError_add * (num_genotypes - num_nas))
  switched_bit_mat <- c(abs(1-sperm_mat_with_na)) #make a switched bit matrix compared to sperm_mat_with_na
  where_locs <- sample(which(!is.na(switched_bit_mat)), size=num_bits_to_flip) #randomly pick num_bits_to_flip locations
  sperm_mat_with_na <- c(sperm_mat_with_na)
  sperm_mat_with_na[where_locs] <- switched_bit_mat[where_locs] #take those locations from switched bit matrix and put them in place in the sperm_mat_with_na matrix
  sperm_mat_with_na <- matrix(sperm_mat_with_na, nrow=num_snps, ncol=num_sperm)
}
after_time <- Sys.time()
time_to_impute <- difftime(after_time, before_time, units="mins")
message(paste0("Time used for adding in sequencing error: ", time_to_impute))

before_time <- Sys.time()
#make it into a dataframe that I can give to the rest of the pipeline, so I need to have genomic positions first, column names for each sperm
sperm_na_df <- data.frame(pseudo_pos = 1:nrow(sperm_mat_with_na), sperm_mat_with_na)
colnames(sperm_na_df) <- c("positions", paste0("sperm", 1:num_sperm, "_"))
sperm_full_df <- data.frame(pseudo_pos = 1:nrow(sperm_mat), sperm_mat)
colnames(sperm_full_df) <- c("positions", paste0("sperm", 1:num_sperm, "_"))
after_time <- Sys.time()
time_to_impute <- difftime(after_time, before_time, units="mins")
message(paste0("Time used for converting to data frames: ", time_to_impute))

before_time <- Sys.time()
#run the whole pipeline
# ensure that each row has at least one 0 and one 1
keep_bool <- unname((rowSums(sperm_na_df[, 2:ncol(sperm_na_df)] == 0, na.rm = TRUE) > 0) & (rowSums(sperm_na_df[, 2:ncol(sperm_na_df)] == 1, na.rm = TRUE) > 0))
sperm_na_df <- sperm_na_df[keep_bool,]
sperm_full_df <- sperm_full_df[keep_bool,]
parental_haps <- parental_haps[keep_bool,]
num_snps <- sum(keep_bool)
message(paste0("new number of snps: ", num_snps))
`%notin%` <- Negate(`%in%`)
for (i in 1:length(new_rows)){
  message(paste0("dnm ", i, " is filtered out: ", new_rows[i] %notin% sperm_na_df[,1]))
}

filename_df <- paste0(outDir, sampleName, "_", chrom, "_pval_sim_truth_filtered.csv")
write_csv(df_counts_pvals_truth[keep_bool,], filename_df)
after_time <- Sys.time()
time_to_impute <- difftime(after_time, before_time, units="mins")
message(paste0("Time used for filtering gamete data: ", time_to_impute))

#stopifnot((floor(num_snps/window_length) >= 4))
#need to figure out how to adjust window_length based on number of snps after filtering

# remove the first column (positions)
positions <- sperm_na_df[, 1]
sperm_na_df <- sperm_na_df[,-1]

# this function gets the mode of a vector after removing the NAs
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv <- uniqv[!is.na(uniqv)]
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# this function replaces 0s with 1s and 1s with 0s in a data frame
invertBits <- function(df) {
  df[df == 0] <- -1
  df[df == 1] <- 0
  df[df == -1] <- 1
  return(df)
}
before_time <- Sys.time()
# overlapping window function from https://stackoverflow.com/questions/8872376/split-vector-with-overlapping-samples-in-r
splitWithOverlap <- function(vec, seg.length, overlap) {
  starts = seq(1, length(vec), by=seg.length-overlap)
  ends   = starts + seg.length - 1
  ends[ends > length(vec)] = length(vec)
  lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
}

# use overlaps of window length/2
windows <- splitWithOverlap(rank(positions), window_length, overlap = window_length / 2)

message(paste0("Number of windows with overlap of ", window_length / 2 , " and ", num_snps, " number of SNPs following filtering: ", length(windows)))

# merge the last two windows to avoid edge effect
combined <- unique(c(windows[[length(windows) - 1]], windows[[length(windows)]]))
combined <- combined[order(combined)]
total_combined <- windows[-c((length(windows) - 1), length(windows))]
total_combined[[length(total_combined) + 1]] <- combined
windows <- total_combined
after_time <- Sys.time()
time_to_impute <- difftime(after_time, before_time, units="mins")
message(paste0("Time used for generating windows: ", time_to_impute))

before_time <- Sys.time()
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
  h2_inferred <- unname(apply(cbind(input_dt[window_start:window_end, h2_sperm],
                                    invertBits(input_dt[window_start:window_end, h1_sperm])),
                              1, function(x) getmode(x)))
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
  if (mean_concordance2 < 0.1) {
    olap_haps_complete$h1.y <- invertBits(olap_haps_complete$h1.y)
  } else if (mean_concordance2 < 0.9) {
    stop(paste0("Haplotypes within overlapping windows are too discordant to merge. Mean: ", mean_concordance2, " Window: ", hap_window))
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
message(paste0("Time used for inferring simulated parental haplotypes: ", time_to_impute))

before_time <- Sys.time()
# Going through each sperm, if an allele (0 or 1) in a sperm matches the allele (0 or 1)
# in h1 at that position, replace the allele with "h1". Do the same for h2.
for (i in 1:ncol(sperm_na_df)) {
  sperm_na_df[i][sperm_na_df[i] == complete_haplotypes$h1] <- "h1"
  sperm_na_df[i][sperm_na_df[i] == complete_haplotypes$h2] <- "h2"
}

# Scan sperm by sperm to interpret state given emission
# First, we initialize our HMM
# set denominator for transition probability - one recombination event per chromosome
#num_snps <- nrow(complete_haplotypes)
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

imputed_sperm <- as_tibble(do.call(cbind, pbmclapply(1:ncol(sperm_na_df),
                                                     function(x) runHMM(sperm_na_df, x),
                                                     mc.cores = getOption("mc.cores", threads))))
after_time <- Sys.time()
time_to_impute <- difftime(after_time, before_time, units = "mins")
message(paste0("Time used for imputing simulated gamete haplotypes: ", time_to_impute))
colnames(imputed_sperm) <- colnames(sperm_na_df)

before_time <- Sys.time()
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
colnames(filled_sperm) <- colnames(sperm_na_df)
after_time <- Sys.time()
time_to_impute <- difftime(after_time, before_time, units="mins")
message(paste0("Time used for filling gametes: ", time_to_impute))

before_time <- Sys.time()
if (!smooth){
  original_dt <- dt %>% 
    mutate_all(funs(str_replace(., "h1", "haplotype1"))) %>%
    mutate_all(funs(str_replace(., "h2", "haplotype2")))
  original_dt <- as.data.frame(original_dt)
  filled_sperm <- as.data.frame(filled_sperm)
  filled_sperm[!is.na(original_dt)] <- original_dt[!is.na(original_dt)]
  filled_sperm <- as_tibble(filled_sperm)
}
after_time <- Sys.time()
time_to_impute <- difftime(after_time, before_time, units="mins")
message(paste0("Time used for un-smoothing: ", time_to_impute))

before_time <- Sys.time()
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
filename_df <- paste0(outDir, sampleName, "_", chrom, "_pval_sim.csv")
write_csv(df_counts_pvals, filename_df)
after_time <- Sys.time()
time_to_impute <- difftime(after_time, before_time, units="mins")
message(paste0("Time used for finding and writing TD in simualted gametes: ", time_to_impute))

before_time <- Sys.time()
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

idents_for_csv <- paste0(paste0(sampleName, "_", chrom, "_"), colnames(filled_sperm))
recomb_spots_all <- do.call(rbind, pbmclapply(1:ncol(filled_sperm),
                                              function(x) find_recomb_spots(filled_sperm, x, idents_for_csv, positions),
                                              mc.cores=getOption("mc.cores", threads))) %>%
  right_join(., tibble(Ident = idents_for_csv), by = "Ident")
after_time <- Sys.time()
time_to_impute <- difftime(after_time, before_time, units="mins")
message(paste0("Time used for finding recombination events: ", time_to_impute))

before_time <- Sys.time()
#Prepare a dataframe for the true crossover locations that can be passed to bedr
sperm_ident <- paste0(rep("sperm", num_sperm), 1:num_sperm, "_")
names(crossover_indices) <- sperm_ident
unlist_ci <- unlist(crossover_indices, use.names=TRUE) #just a list with a single row for each crossover, rowname corresponds to the sperm _1, _2, etc for sperm with more than one CO
#bedr has to have different starts and ends and cannot handle NAs
ci_df <- data.frame(chr=rep(as.character(chrom), length(unlist_ci)), start=(unlist_ci-1), end=(unlist_ci))
row.names(ci_df) <- names(unlist_ci)
truth_nona_df <- drop_na(ci_df) #this one has all valid regions for bedr
truth_onlyna_df <- ci_df[is.na(ci_df$start),] #this one is to check all the ones that don't have any recombination spots
truth_nona_df_sort <- bedr.sort.region(truth_nona_df) #these are the true recombination spots, with no NAs, lexicographically sorted


#make another bedr compatible dataframe for the predicted crossover locations
split_idents <- str_split(recomb_spots_all$Ident, "_", simplify=TRUE)
recomb_spots_df <- data.frame(chr=as.character(split_idents[,2]), start=recomb_spots_all$Genomic_start, end=recomb_spots_all$Genomic_end)
row.names(recomb_spots_df) <- make.names(paste0(split_idents[,3], "_"), unique=TRUE)
#bedr cannot handle NAs
pred_nona_df <- drop_na(recomb_spots_df) #this one has all valid regions for bedr
pred_onlyna_df <- recomb_spots_df[is.na(recomb_spots_df$start),] #this one is to check all the ones that don't have any recombination spots
pred_nona_df_sort <- bedr.sort.region(pred_nona_df) #these are the predicted recombination spots, with no NAs, lexicographically sorted
after_time <- Sys.time()
time_to_impute <- difftime(after_time, before_time, units="mins")
message(paste0("Time used for converting to bedr compatible data frames: ", time_to_impute))

before_time <- Sys.time()
##plot the resolution of predicted recombination break point regions
filename = paste0(outDir, "simulated_recombination_resolution.pdf")
pdf(file=filename)
hist(recomb_spots_df$end - recomb_spots_df$start, xlab = "Differince in SNP index", breaks=100, main="Resolution of predicted recombination break point regions")
dev.off()
after_time <- Sys.time()
plotting_time <- difftime(after_time, before_time, units="mins")


before_time <- Sys.time()
truth_pred_int <- bedr(input = list(a=truth_nona_df_sort, b=pred_nona_df_sort), 
                  method = "intersect",
                  params = "-loj -sorted") #those that aren't overlapping will be reported as NULL for b (. -1 -1). Those are fn
colnames(truth_pred_int) <- c("truth_chr", "truth_start", "truth_end", "pred_chr", "pred_start", "pred_end")
pred_truth_int <- bedr(input = list(a=pred_nona_df_sort, b=truth_nona_df_sort),
                       method = "intersect",
                       params = "-loj -sorted")
colnames(pred_truth_int) <- c("pred_chr", "pred_start", "pred_end", "truth_chr", "truth_start", "truth_end")
pred_v <- pred_truth_int[pred_truth_int$truth_chr == ".",] #nrow of this is fp

truth_pred_int_nona <- truth_pred_int[!truth_pred_int$pred_chr == ".",]
pred_idents <- paste0(truth_pred_int_nona$pred_chr, "_", truth_pred_int_nona$pred_start, "_", truth_pred_int_nona$pred_end)
after_time <- Sys.time()
time_to_impute <- difftime(after_time, before_time, units="mins")
message(paste0("Time used for running bedr: ", time_to_impute))

before_time <- Sys.time()
tp_cons <- length(unique(pred_idents))
tp_lib <- nrow(truth_pred_int)

fp <- nrow(pred_v)

fn_pt3 <- tp_lib - nrow(truth_pred_int_nona)
fn_pt1 <- nrow(truth_pred_int[truth_pred_int$pred_chr == ".",])
fn_pt2 <- length(setdiff(row.names(pred_onlyna_df), row.names(truth_onlyna_df))) #rownames in the prediction but not in the truth
fn_lib <- fn_pt1 + fn_pt2
fn_cons <- fn_pt1 + fn_pt2 + fn_pt3

tn <- nrow(merge(truth_onlyna_df, pred_onlyna_df, by=0))

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
after_time <- Sys.time()
time_to_impute <- difftime(after_time, before_time, units="mins")
message(paste0("Time used for finding and reporting recombination metrics: ", time_to_impute))

before_time <- Sys.time()
real_reads <- rowSums(!is.na(sperm_na_df))
filename = paste0(outDir, "simulated_notna_bysnp.pdf")
pdf(file=filename)
hist(real_reads, breaks=150, xlab="not NA by SNP", main="Distribution of simulated data not NA by SNP")
dev.off()

real_reads2 <- colSums(!is.na(sperm_na_df))
filename = paste0(outDir, "simulated_notna_bysperm.pdf")
pdf(file=filename)
hist(real_reads2, breaks=100, xlab="not NA by sperm", main="Distribution of simulated data not NA by sperm")
dev.off()
after_time <- Sys.time()
plotting_time <- plotting_time + difftime(after_time, before_time, units="mins")

before_time <- Sys.time()
## Assessing the accuracy of parental haplotype reconstruction
#complete_haplotypes compared to simulated hap1 and hap2
num_mismatch_parental1 <- min(sum((complete_haplotypes$h1 - parental_haps$Parental1) != 0, na.rm=TRUE), sum((complete_haplotypes$h1 - parental_haps$Parental2) != 0, na.rm=TRUE))
accuracy_parental1 <- (num_snps - num_mismatch_parental1) / num_snps * 100
message(paste0("Parental haplotype reconstruction accuracy: ", accuracy_parental1))
if (sum(is.na(complete_haplotypes$h1))> 0){
  num_na_h1 <- sum(is.na(complete_haplotypes$h1))
  accuracy_parental1_cor <- (((num_snps - num_na_h1) - num_mismatch_parental1)/(num_snps - num_na_hap1)) * 100
  message(past0("Parental haplotype reconstruction accuracy, corrected: ", accuracy_parental1_cor))
}
after_time <- Sys.time()
time_to_impute <- difftime(after_time, before_time, units="mins")
message(paste0("Time used for assessing parental haplotype reconstruction accuracy: ", time_to_impute))

before_time <- Sys.time()
mismatches_parental <- which((complete_haplotypes$h1 -  parental_haps$Parental1)  != 0)
filename = paste0(outDir, "simulated_parental_mismatch_loc.pdf")
pdf(file=filename)
plot(c(mismatches_parental, new_rows), c(rep(1, length(mismatches_parental)), rep(1.05, length(new_rows))), cex=0.2, col=c(rep(rgb(red=0, green=0, blue=0, alpha=0.2), length(mismatches_parental)), rep(rgb(red=1, green=0, blue=0, alpha=0.8), length(new_rows))))
dev.off()
after_time <- Sys.time()
plotting_time <- plotting_time + difftime(after_time, before_time, units="mins")

before_time <- Sys.time()
## Assessing the accuracy of gamete haplotype reconstruction
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

filled_sperm_recode <- re_recode_gametes(filled_sperm, complete_haplotypes)
num_mismatches_sperm_haplotype <- colSums((filled_sperm_recode - sperm_full_df[,-1]) != 0, na.rm = TRUE)
num_nas_byCol <- colSums(is.na(filled_sperm_recode))
rawAccuracy <- data.frame(val=((num_snps - num_mismatches_sperm_haplotype)/num_snps * 100), name="Raw") 
message(paste0("Mean sperm haplotype reconstruction accuracy (raw): ", mean(rawAccuracy$val)))
correctedAccuracy <- data.frame(val=(((num_snps - num_nas_byCol) - num_mismatches_sperm_haplotype)/(num_snps - num_nas_byCol) * 100), name="Corrected")
message(paste0("Mean sperm haplotype reconstruction accuracy (corrected): ", mean(correctedAccuracy$val, na.rm=TRUE)))
accuracyDat <- rbind(rawAccuracy, correctedAccuracy)
after_time <- Sys.time()
time_to_impute <- difftime(after_time, before_time, units="mins")
message(paste0("Time used for assessing gamete haplotype reconstruction accuracies: ", time_to_impute))

before_time <- Sys.time()
filename = paste0(outDir, "simulated_sperm_hap_reconstruction.pdf")
pdf(file=filename)
ggplot(accuracyDat, aes(x=name, y=val, fill=name)) + scale_x_discrete(limits=c("Raw", "Corrected")) + geom_violin(scale="width", adjust=1, width=0.5) + labs(y="accuracy", x="method", title="Sperm Haplotype Reconstruction") + scale_y_continuous(name="accuracy", breaks=c(98.75, 99, 99.25, 99.5, 99.75, 100), labels=c(98.75, 99, 99.25, 99.5, 99.75, 100), limits=c(98.75, 100))
dev.off()
after_time <- Sys.time()
plotting_time <- plotting_time + difftime(after_time, before_time, units="mins")
message(paste0("Time used for plotting: ", plotting_time))
end_script_time <- Sys.time()
script_time <- difftime(end_script_time, start_script_time, units = "mins")
message(paste0("Time used for the script: ", script_time))
