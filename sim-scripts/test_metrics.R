#date; time Rscript test_metrics.R 2 3000 1000 30000 0.01 42 1 &> run_sim3_20210218/wl3000_gam1000_snp30000_cov0.01/rsd42.txt
#date; time Rscript test_metrics.R 2 3000 1000 30000 0.01 386 1 &> run_sim3_20210218/wl3000_gam1000_snp30000_cov0.01/rsd386.txt
#date; time Rscript test_metrics.R 2 3000 1000 30000 0.01 1059 1 &> run_sim3_20210218/wl3000_gam1000_snp30000_cov0.01/rsd1059.txt
#date; time Rscript test_metrics.R 2 3000 1000 30000 0.01 27 1 &> run_sim3_20210218/wl3000_gam1000_snp30000_cov0.01/rsd27.txt
#date; time Rscript test_metrics.R 2 3000 1000 30000 0.01 651 1 &> run_sim3_202102018/wl3000_gam1000_snp30000_cov0.01/rsd651.txt
#date; time Rscript test_metrics.R 2 3000 1000 30000 0.01 2862 1 &> run_sim3_202102018/wl3000_gam1000_snp30000_cov0.01/rsd2862.txt
#date; time Rscript test_metrics.R 2 3000 1000 30000 0.01 2556 1 &> run_sim3_202102018/wl3000_gam1000_snp30000_cov0.01/rsd2556.txt
#date; time Rscript test_metrics.R 2 3000 1000 30000 0.01 2563 1 &> run_sim3_202102018/wl3000_gam1000_snp30000_cov0.01/rsd2563.txt
#date; time Rscript test_metrics.R 2 3000 1000 30000 0.01 3417 1 &> run_sim3_202102018/wl3000_gam1000_snp30000_cov0.01/rsd3417.txt
#date; time Rscript test_metrics.R 2 3000 1000 30000 0.01 4900 1 &> run_sim3_202102018/wl3000_gam1000_snp30000_cov0.01/rsd4900.txt



library(data.table)
library(tidyverse)
library(stringr)
library(pbapply)
library(pbmcapply)
library(HMM)
library(bedr)

args <- commandArgs(trailingOnly = TRUE)
sampleName <- "simTest"
chrom <- "chrT"
threads <- 2L #as.integer(args[1])

window_length <- 300 #as.integer(args[2])
seqError <- 0.05
hapProb <- 1 - seqError

num_sperm <- 3 #as.integer(args[3])
num_snps <- 100000 #as.integer(args[4])
coverage <- 0.1 #as.numeric(args[5])

num_genotypes <- num_sperm * num_snps
missing_genotype_rate <- dpois(0, coverage)
message(paste0("The missing genotype rate of this simulation is ", missing_genotype_rate))
num_nas <- floor(num_genotypes * missing_genotype_rate)

random_seed <- as.integer(args[1])
set.seed(random_seed)

recomb_lambda <- 1 #as.integer(args[7])
num_recomb_sites <- rpois(num_sperm, recomb_lambda)
message(paste0("Total number of recombination spots across gametes: ", sum(num_recomb_sites)))

add_seq_error <- TRUE
seqError_add <- 0.05

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

###Sparsify the knowledge of the gametes
add_to_na_flatten <- function(to_add_from, num_nas, num_sperm, num_snps){
  to_return <- rep(NA, (num_snps*num_sperm))
  coords_to_keep_genotype <- sample(1:(num_snps*num_sperm), size=((num_snps*num_sperm)-num_nas), replace = FALSE)
  to_return[coords_to_keep_genotype] <- as.vector(to_add_from[coords_to_keep_genotype])
  to_return <- matrix(to_return, ncol=num_sperm, nrow=num_snps)
  return (to_return)
}

sperm_mat_with_na <- add_to_na_flatten(sperm_mat, num_nas, num_sperm, num_snps)

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
message(paste0("new number of snps: ", num_snps))

##Remove the first column (positions)
positions <- sperm_na_df[, 1]
sperm_na_df <- sperm_na_df[,-1]

# this function gets the mode of a vector after removing the NAs
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv <- uniqv[!is.na(uniqv)]
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#getmode <- function(x) { #from https://stackoverflow.com/questions/56552709/r-no-mode-and-exclude-na?noredirect=1#comment99692066_56552709
#  ux <- unique(na.omit(x))
#  tx <- tabulate(match(x, ux))
#  if(length(ux) != 1 & sum(max(tx) == tx) > 1) {
#    if (is.character(ux)) return(NA_character_) else return(NA_real_)
#  }
#  max_tx <- tx == max(tx)
#  return(ux[max_tx])
#}



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
message(length(windows))


# merge the last two windows to avoid edge effect
if (length(windows) > 1){
  combined <- unique(c(windows[[length(windows) - 1]], windows[[length(windows)]]))
  combined <- combined[order(combined)]
  total_combined <- windows[-c((length(windows) - 1), length(windows))]
  total_combined[[length(total_combined) + 1]] <- combined
  windows <- total_combined
}

message(paste0("Number of windows with overlap of ", window_length / 2 , " and ", num_snps, " number of SNPs following filtering: ", length(windows)))

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
  message(paste0("Window ", hap_window, "Mean Concordance: ", mean_concordance2))
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

###Assessing the accuracy of parental haplotype reconstruction
#complete_haplotypes compared to simulated hap1 and hap2
num_mismatch_parental1 <- min(sum((complete_haplotypes$h1 - parental_haps$Parental1) != 0, na.rm=TRUE), sum((complete_haplotypes$h1 - parental_haps$Parental2) != 0, na.rm=TRUE))
num_nas_parental1 <- min(sum(is.na(complete_haplotypes$h1 - parental_haps$Parental1)), sum(is.na(complete_haplotypes$h1 - parental_haps$Parental2)))
accuracy_parental1 <- (num_snps - (num_mismatch_parental1 + num_nas_parental1))/ num_snps * 100
message(paste0("Parental haplotype reconstruction accuracy: ", accuracy_parental1))
if (sum(is.na(complete_haplotypes$h1))> 0){
  num_na_h1 <- sum(is.na(complete_haplotypes$h1))
  accuracy_parental1_cor <- (((num_snps - num_na_h1) - num_mismatch_parental1)/(num_snps - num_na_h1)) * 100
  message(paste0("Parental haplotype reconstruction accuracy, corrected: ", accuracy_parental1_cor))
}

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
colnames(imputed_sperm) <- colnames(sperm_na_df)

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

filled_sperm_recode <- re_recode_gametes(filled_sperm, complete_haplotypes)
num_mismatches_sperm_haplotype <- colSums((filled_sperm_recode - sperm_full_df[,-1]) != 0, na.rm = TRUE)
num_nas_byCol <- colSums(is.na(filled_sperm_recode))
rawAccuracy <- data.frame(val=((num_snps - (num_mismatches_sperm_haplotype + num_nas_byCol))/num_snps * 100), name="Raw") 
message(paste0("Mean sperm haplotype reconstruction accuracy (raw): ", mean(rawAccuracy$val, na.rm=TRUE)))
correctedAccuracy <- data.frame(val=(((num_snps - num_nas_byCol) - num_mismatches_sperm_haplotype)/(num_snps - num_nas_byCol) * 100), name="Corrected")
message(paste0("Mean sperm haplotype reconstruction accuracy (corrected): ", mean(correctedAccuracy$val, na.rm=TRUE)))
accuracyDat <- rbind(rawAccuracy, correctedAccuracy)

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

###Prepare a truth dataframe for the true crossover locations that is bedr compatible
##Note bedr has to have different start and end values and cannot handle NAs
canrun_bedr <- TRUE
no_truths <- FALSE
no_preds <- FALSE

names(crossover_indices) <- paste0(rep("sperm", num_sperm), 1:num_sperm, "_")
unlist_ci <- unlist(crossover_indices, use.names=TRUE) #just a list with a single row for each crossover, rowname corresponds to the sperm _1, _2, etc for sperm with more than one CO
ci_df <- data.frame(chr=rep(as.character(chrom), length(unlist_ci)), start=(unlist_ci-1), end=(unlist_ci))
row.names(ci_df) <- names(unlist_ci)
truth_onlyna_df <- ci_df[is.na(ci_df$start),] #this one is to check all the gametes that don't have any true recombination spots
if (nrow(drop_na(ci_df)) >=1 ){
  truth_nona_df_sort <- bedr.sort.region(drop_na(ci_df))
  truth_nona_df_sort <- add_column(truth_nona_df_sort, "gametes" = row.names(truth_nona_df_sort),  "gametes_simp"=gsub("_.*", "", row.names(truth_nona_df_sort)))  #these are the true recombination spots, valid regions with no NAs, lexicographically sorted
} else {
  canrun_bedr <- FALSE
  no_truths <- TRUE
}
  
###Prepare a prediction  dataframe for the predicted crossover locations that is bedr compatible
split_idents <- str_split(recomb_spots_all$Ident, "_", simplify=TRUE)
recomb_spots_df <- data.frame(chr=as.character(split_idents[,2]), start=recomb_spots_all$Genomic_start, end=recomb_spots_all$Genomic_end)
row.names(recomb_spots_df) <- make.names(paste0(split_idents[,3], "_"), unique=TRUE)
pred_onlyna_df <- recomb_spots_df[is.na(recomb_spots_df$start),] #this one is to check all the gametes that don't have any predicted recombination spots
#might need to add a check to make sure these are greater than 0 rows

if (nrow(drop_na(recomb_spots_df)) >= 1) {
  pred_nona_df_sort <- bedr.sort.region(drop_na(recomb_spots_df))
  pred_nona_df_sort <- add_column(pred_nona_df_sort, "gametes"= row.names(pred_nona_df_sort),  "gametes_simp"=gsub("_.*", "", row.names(pred_nona_df_sort)))#these are the predicted recombination spots, valid regions with no NAs, lexicographically sorted
} else {
  canrun_bedr <- FALSE
  no_preds <- TRUE
}

if (canrun_bedr) {
  truth_pred_int <- bedr(input = list(a=truth_nona_df_sort, b=pred_nona_df_sort), 
                         method = "intersect",
                         params = "-loj -sorted") #those that aren't overlapping will be reported as NULL for b (. -1 -1). Those are fn
  colnames(truth_pred_int) <- c("truth_chr", "truth_start", "truth_end", "truth_gametes", "truth_gametes_simp", "pred_chr", "pred_start", "pred_end", "pred_gametes",  "pred_gametes_simp")
  truth_pred_int_nona <- truth_pred_int[!truth_pred_int$pred_chr == ".",]
  rows_to_keep  <- which(truth_pred_int_nona$truth_gametes_simp == truth_pred_int_nona$pred_gametes_simp)
  truth_pred_int_nona  <- truth_pred_int_nona[rows_to_keep, ]

  fn_pt1 <- nrow(truth_pred_int[truth_pred_int$pred_chr == ".",]) #false negative because there is no overlapping predicted recombination for a true recombination spot, accomplished by finding the number of locations that are null in prediction column when truth intersected with prediction
  fn_pt2 <- length(setdiff(row.names(pred_onlyna_df), row.names(truth_onlyna_df))) #predicted returns only NA when there is at least one truth rownames in the prediction but not in the truth

  tp_cons <- length(unique(paste0(truth_pred_int_nona$pred_chr, "_", truth_pred_int_nona$pred_start, "_", truth_pred_int_nona$pred_end))) #want to count only one truth when multiple truths intersect a single prediction; accomplish this by the length of the unique prediction identities from the non NA intersection
  tp_lib <- nrow(truth_pred_int_nona) #want to count +1 for every truth that intersects a prediction, even if multiple truths intersect a single prediction. accomplish this by counting number of truths that intersect any prediction from the non NA intersection

  pred_truth_int <- bedr(input = list(a=pred_nona_df_sort, b=truth_nona_df_sort),
                       method = "intersect",
                       params = "-loj -sorted")
  colnames(pred_truth_int) <- c("pred_chr", "pred_start", "pred_end",  "pred_gametes",  "pred_gametes_simp", "truth_chr", "truth_start", "truth_end", "truth_gametes", "truth_gametes_simp")
  fp <- nrow(pred_truth_int[pred_truth_int$truth_chr == ".",]) #number predicted as recombination spots, but don't intersect with truth at all

  fn_lib <- fn_pt1 + fn_pt2
  fn_cons <- fn_pt1 + fn_pt2 + (tp_lib - nrow(truth_pred_int_nona)) #include overflow intersections, or when multiple truths intersect a single prediction 

  tn <- nrow(merge(truth_onlyna_df, pred_onlyna_df, by=0))
} else {
  if (no_truths & !no_preds){ #no non-na truths, but there are non-na predictions, all are false positives
    fp <- nrow(pred_nona_df_sort)
    fn <- length(setdiff(row.names(pred_onlyna_df), row.names(truth_onlyna_df)))
    tp <- 0
    tn <- nrow(merge(truth_onlyna_df, pred_onlyna_df, by=0))
  } 
  else if (no_preds & !no_truths){ #no non-na predictions but there are non-na truths, all are false negatives
    fn <- nrow(truth_nona_df_sort)
    fp <- 0
    tn <- nrow(merge(truth_onlyna_df, pred_onlyna_df, by=0))
    tp <- 0
  } 
  
}
  

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

if (no_truths & no_preds){
  message("no recombination metrics possible because there were no truths or predictions")
} else if (no_truths | no_preds){
  metrics_it <- metrics(tp, fp, tn, fn)
  message(paste0("Conservative recombination spot identification metrics\nPrecision: ", metrics_it$precision,
                 "\nRecall: ", metrics_it$recall,
                 "\nAccuracy: ", metrics_it$accuracy,
                 "\nF1: ", metrics_it$f1,
                 "\nSpecificity: ", metrics_it$specificity,
                 "\nFDR: ", metrics_it$fdr,
                 "\nFPR: ", metrics_it$fpr))
  message(paste0("Liberal recombination spot identification metrics\nPrecision: ", metrics_it$precision,
                 "\nRecall: ", metrics_it$recall,
                 "\nAccuracy: ", metrics_it$accuracy,
                 "\nF1: ", metrics_it$f1,
                 "\nSpecificity: ", metrics_it$specificity,
                 "\nFDR: ", metrics_it$fdr,
                 "\nFPR: ", metrics_it$fpr))
  
} else {
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
}