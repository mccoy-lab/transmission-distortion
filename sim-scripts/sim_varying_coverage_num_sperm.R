library(data.table)
library(tidyverse)
library(bedr)
library(stringr)
library(pbapply)
library(pbmcapply)
library(HMM)

#args <- commandArgs(trailingOnly = TRUE)
#sampleName <- args[1]
sampleName <- "test"

#chrom <- args[2]
chrom <- "chrT"

#outDir <- args[3]
outDir <- "tmp"

#seqError <- as.numeric(args[4]) #0.005
seqError <- 0.005
hapProb <- 1 - seqError

#threads <- as.integer(args[5])
threads <- 2L

#window_length <- as.integer(args[6]) #2500 default
window_length <- 2500

#num_sperm <- as.integer(args[7])
num_sperm <- 1000 #coverage decreases as this number increases

#num_snps <- as.integer(args[8])
num_snps <- 30000 #coverage doesn't seem to change with this number

#random_seed <- as.integer(args[9])
random_seed <- 42
set.seed(random_seed)

#lambda <- as.numeric(args[10])
lambda <- 1
num_recomb_sites <- rpois(num_sperm, lambda)

#num_not_nan_per_row <- as.integer(args[11])
num_not_nan_per_row <- 3 #coverage increases as this number increases
#set up a check that verifies that num_not_nan_per_row is at least 2, otherwise we're not working with heterozygous SNPs
stopifnot(num_not_nan_per_row >= 2)

#coverage <- as.numeric(args[11])
#coverage <- 0.001
#missing_genotype_rate <- dpois(0, coverage)

#the poisson distribution function, dpois, is P(x;mu) = e**(-mu) * mu**x / x! Therefore, when x=0 and P(x;mu) is known, mu is -ln(P(x;mu)). Here mu is coverage and we want to find it

num_nas_per_row <- num_sperm - num_not_nan_per_row
num_nas <- num_nas_per_row * num_snps
num_genotypes <- num_sperm * num_snps
missing_genotype_rate <- num_nas / num_genotypes

coverage <- -log(missing_genotype_rate)
message(paste0("The coverage of the this simulation is ", coverage))


#start with 2 parental chromosomes of heterozygous sites
#first parental chromosome
hap1 <- data.frame(V1 = sample(c(0, 1), size = num_snps, replace = TRUE))
hap2 <- abs(1-hap1)
parental_haps <- data.frame(cbind(hap1, hap2))
colnames(parental_haps) <- c("Parental1", "Parental2")

indices <- 1:num_snps #we'll use these as both indices and a pseudo genomic position

generate_sperm <- function(parental_haplotypes, n_crossovers){
  init_hap_index <- sample(1:2, 1)
  init_hap <- parental_haplotypes[,init_hap_index]
  if (n_crossovers == 0) {
    return(list(NA, init_hap, init_hap_index))
  } else {
    n_snps <- length(init_hap)
    crossover_indices <- sample(1:n_snps, n_crossovers)
    crossover_indices <- crossover_indices[order(crossover_indices)]
    recombined_hap <- init_hap
    for (crossover in crossover_indices) { #switch haplotypes at the crossover index
      recombined_hap <- c(recombined_hap[1:(crossover-1)], abs(1- recombined_hap[crossover:n_snps]))
    }
    return(list(crossover_indices, recombined_hap, init_hap_index))
  }
}

sim_sperm <- lapply(1:num_sperm, function(x) generate_sperm(parental_haps, num_recomb_sites[x]))
sperm_mat <- sapply(sim_sperm, "[[", 2)
first_haps <- sapply(sim_sperm, "[[", 3)
crossover_indices <- sapply(sim_sperm, "[[", 1)

add_to_na <- function(to_add_from, num_not_nas_perRow, num_sperm, x){
  to_return = rep(NA, num_sperm)
  where_0 = sample(which(to_add_from[x,]==0), size=1)
  where_1 = sample(which(to_add_from[x,]==1), size=1)
  to_return[c(where_0, where_1)] <- to_add_from[x,c(where_0, where_1)]
  if(num_not_nas_perRow -2 >0){
    to_keep = sample(c(1:num_sperm)[-c(where_0, where_1)], size=num_not_nas_perRow-2)
    to_return[to_keep] <- to_add_from[x,to_keep]
  }
  return (to_return)
}

sperm_mat_with_na <- do.call(rbind, pbmclapply(1:num_snps,
                                               function (x) add_to_na(sperm_mat, num_not_nan_per_row, num_sperm, x),
                                               mc.cores=getOption("mc.cores", threads)))

#make it into a dataframe that I can give to the rest of the pipeline, so I need to have genomic positions first, column names for each sperm
sperm_na_df <- data.frame(pseudo_pos = 1:nrow(sperm_mat_with_na), sperm_mat_with_na)
sperm_full_df <- data.frame(pseudo_pos = 1:nrow(sperm_mat), sperm_mat)
#I think sperm_na_df is the df that can be passed to assign_sperm_haplotypes_rm_kw.R directly

# #run the whole pipeline
# # remove the first column (positions)
# positions <- sperm_na_df[, 1]
# sperm_na_df <- sperm_na_df[,-1]
# 
# # this function gets the mode of a vector after removing the NAs
# getmode <- function(v) {
#   uniqv <- unique(v)
#   uniqv <- uniqv[!is.na(uniqv)]
#   uniqv[which.max(tabulate(match(v, uniqv)))]
# }
# 
# # this function replaces 0s with 1s and 1s with 0s in a data frame
# invertBits <- function(df) {
#   df[df == 0] <- -1
#   df[df == 1] <- 0
#   df[df == -1] <- 1
#   return(df)
# }
# 
# # overlapping window function from https://stackoverflow.com/questions/8872376/split-vector-with-overlapping-samples-in-r
# splitWithOverlap <- function(vec, seg.length, overlap) {
#   starts = seq(1, length(vec), by=seg.length-overlap)
#   ends   = starts + seg.length - 1
#   ends[ends > length(vec)] = length(vec)
#   lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
# }
# 
# # use overlaps of window length/2
# windows <- splitWithOverlap(rank(positions), window_length, overlap = window_length / 2)
# 
# # function to reconstruct parental haplotypes
# reconstruct_hap <- function(input_dt, input_positions, window_indices) {
#   window_start <- min(window_indices)
#   window_end <- max(window_indices)
#   positions_for_window <- input_positions[window_start:window_end]
#   # compute a distance matrix
#   d <- dist(t(as.matrix(input_dt)[window_start:window_end,]), method = "binary")
#   # plut in 0.5 for any NA entries of the distance matrix
#   d[is.na(d)] <- 0.5
#   # cluster the distance matrix
#   tree <- hclust(d, method = "ward.D2")
#   # plot(tree, cex = 0.1) # uncomment to plot
#   # cut the tree generated by clustering into two groups (haplotypes)
#   haplotypes <- cutree(tree, k=2)
#   # get the names of the sperm cells falling into the two groups
#   h1_sperm <- names(haplotypes[haplotypes == 1])
#   h2_sperm <- names(haplotypes[haplotypes == 2])
#   # reconstruct the original haplotypes by majority vote after inverting the opposite haplotype
#   h1_inferred <- unname(apply(cbind(input_dt[window_start:window_end, h1_sperm],
#                                     invertBits(input_dt[window_start:window_end, h2_sperm])),
#                               1, function(x) getmode(x)))
#   h2_inferred <- unname(apply(cbind(input_dt[window_start:window_end, h2_sperm],
#                                     invertBits(input_dt[window_start:window_end, h1_sperm])),
#                               1, function(x) getmode(x)))
#   return(tibble(index = window_indices, pos = positions_for_window, h1 = h1_inferred))
# }
# 
# # infer the haplotypes within the overlapping windows
# inferred_haplotypes <- pbmclapply(1:length(windows), 
#                                   function(x) reconstruct_hap(sperm_na_df, positions, windows[[x]]),
#                                   mc.cores = getOption("mc.cores", threads))
# 
# # stitch together the haplotypes
# initial_haplotype <- inferred_haplotypes[[1]]
# for (hap_window in 2:length(windows)) {
#   olap_haps <- merge(initial_haplotype, inferred_haplotypes[[hap_window]], by = "index")
#   olap_haps_complete <- merge(initial_haplotype, inferred_haplotypes[[hap_window]], by = "index", all = TRUE)
#   head(olap_haps$h1.x)
#   head(olap_haps$h1.y)
#   mean_concordance <- mean(olap_haps$h1.x == olap_haps$h1.y)
#   head(olap_haps$h1.x == olap_haps$h1.y)
#   message(mean_concordance)
# }
#   if (mean_concordance < 0.1) {
#     olap_haps_complete$h1.y <- invertBits(olap_haps_complete$h1.y)
#   } else if (mean_concordance < 0.9) {
#     error(paste0("Haplotypes within overlapping windows are too discordant to merge. Mean: ", mean_concordance))
#   }
#   initial_haplotype <- tibble(index = olap_haps_complete$index,
#                               pos = c(olap_haps_complete[is.na(olap_haps_complete$pos.y),]$pos.x,
#                                       olap_haps_complete[!is.na(olap_haps_complete$pos.x) &
#                                                            !is.na(olap_haps_complete$pos.y),]$pos.x,
#                                       olap_haps_complete[is.na(olap_haps_complete$pos.x),]$pos.y),
#                               h1 = c(olap_haps_complete[is.na(olap_haps_complete$pos.y),]$h1.x,
#                                      olap_haps_complete[!is.na(olap_haps_complete$pos.x) &
#                                                           !is.na(olap_haps_complete$pos.y),]$h1.x,
#                                      olap_haps_complete[is.na(olap_haps_complete$pos.x),]$h1.y))
# }
# 
# complete_haplotypes <- initial_haplotype %>%
#   mutate(h2 = invertBits(h1))
# # Going through each sperm, if an allele (0 or 1) in a sperm matches the allele (0 or 1)
# # in h1 at that position, replace the allele with "h1". Do the same for h2.
# for (i in 1:ncol(sperm_na_df)) {
#   sperm_na_df[i][sperm_na_df[i] == complete_haplotypes$h1] <- "h1"
#   sperm_na_df[i][sperm_na_df[i] == complete_haplotypes$h2] <- "h2"
# }
# 
# # Scan sperm by sperm to interpret state given emission
# # First, we initialize our HMM
# # set denominator for transition probability - one recombination event per chromosome
# num_snps <- nrow(complete_haplotypes)
# # two states
# states <- c("haplotype1", "haplotype2")
# # probability of state at position x+1 given state at position x
# hap1Prob <- c(1-(1/num_snps), 1/num_snps)
# hap2Prob <- c(1/num_snps, 1-(1/num_snps))
# transProb <- matrix(c(hap1Prob, hap2Prob), 2)
# # Two emissions (observations): an allele from h1 or an allele from h2
# emissions <- c("h1","h2")
# # Prob of emitting an h1 allele, prob of emitting an h2 allele in state `haplotype1`
# h1ProbEmiss <- c(hapProb, seqError)
# # Prob of emitting an h1 allele, prob of emitting an h2 allele in state `haplotype2`
# h2ProbEmiss <- c(seqError, hapProb)
# emissProb <- matrix(c(h1ProbEmiss, h2ProbEmiss), 2)
# #build model with the above inputs
# hmm <- initHMM(States = states,
#                Symbols = emissions,
#                transProbs = transProb,
#                emissionProbs = emissProb)
# ###### Function to run HMM on an input file
# # Compute the inferred state using each sperm cell as the input
# # (Input must be a vector)
# runHMM <- function(sperm_dt, column_index) {
#   original_obs <- sperm_dt[,column_index]
#   inferred_state <- viterbi(hmm, na.omit(sperm_dt[, column_index]))
#   original_obs[!is.na(original_obs)] <- inferred_state
#   return(original_obs)
# }
# 
# imputed_sperm <- as_tibble(do.call(cbind, pbmclapply(1:ncol(sperm_na_df),
#                                                      function(x) runHMM(sperm_na_df, x),
#                                                      mc.cores = getOption("mc.cores", threads))))
# colnames(imputed_sperm) <- colnames(sperm_na_df)
# 
# # Works on our sperm! Need to make the function work on every sperm in test3
# # and need to make fill up and down at the end
# fill_NAs <- function(merged_sperm, col_index) {
#   sperm_sample <- merged_sperm[,col_index] %>%
#     rename(sperm = colnames(.)[1]) %>%
#     mutate(sperm_up = sperm) %>%
#     mutate(sperm_down = sperm) %>%
#     fill(sperm_up, .direction = "up") %>%
#     fill(sperm_down, .direction = "down") %>%
#     mutate(is_match = (sperm_up == sperm_down)) %>%
#     replace_na(list(is_match = FALSE))
#   sperm_sample$sperm_imputed <- as.character(NA)
#   sperm_sample[sperm_sample$is_match == TRUE,]$sperm_imputed <- sperm_sample[sperm_sample$is_match == TRUE,]$sperm_up
#   #fill beginning of chromosome NA's
#   first <- which(!is.na(sperm_sample$sperm_imputed))[1]
#   sperm_sample$sperm_imputed[1:(first-1)] <- sperm_sample$sperm_imputed[first]
#   #fill end of chromosome NA's
#   sperm_sample$sperm_imputed <- rev(sperm_sample$sperm_imputed)
#   first <- which(!is.na(sperm_sample$sperm_imputed))[1]
#   sperm_sample$sperm_imputed[1:(first-1)] <- sperm_sample$sperm_imputed[first]
#   #reverse chromosome imputation back so it faces the right way
#   sperm_sample$sperm_imputed <- rev(sperm_sample$sperm_imputed)
#   return(sperm_sample$sperm_imputed)
# }
# 
# filled_sperm <- as_tibble(do.call(cbind, 
#                                   pblapply(1:ncol(imputed_sperm),
#                                            function(x) fill_NAs(imputed_sperm, x))))
# colnames(filled_sperm) <- colnames(sperm_na_df)
# 
# td_test <- function(sperm_matrix, row_index) {
#   test_row <- sperm_matrix[row_index,]
#   gt_vector <- unlist(test_row)[-1]
#   one_count <- sum(gt_vector == "haplotype1", na.rm = TRUE)
#   two_count <- sum(gt_vector == "haplotype2", na.rm = TRUE)
#   p_value <- binom.test(c(one_count, two_count))$p.value
#   return(c(p_value, one_count, two_count))
# }
# 
# #df_counts_pvals <- do.call(rbind, pbmclapply(1:nrow(filled_sperm), 
# #                                             function(x) td_test(filled_sperm, x),
# #                                             mc.cores=getOption("mc.cores", threads))) %>%
# #  as_tibble() %>% 
# #  add_column(positions) #bind the positions vector to df_counts_pvals
# #colnames(df_counts_pvals) <- c("pval", "h1_count", "h2_count", "genomic_position")
# #filename_df <- paste0(outDir, sampleName, "_", chrom, "_pval.csv")
# #write_csv(df_counts_pvals, filename_df)
# 
# #find recombination spots
# find_recomb_spots <- function(input_matrix, x, identities, genomic_positions){
#   ident <- identities[x]
#   input_tibble <- input_matrix[, x] %>%
#     mutate(., index = row_number()) %>%
#     mutate(., positions = genomic_positions)
#   complete_cases_tibble <- input_tibble[complete.cases(input_tibble),]
#   input_vec <- as.factor(complete_cases_tibble[[1]])
#   switch_indices <- which(input_vec[-1] != input_vec[-length(input_vec)])
#   switch_indices_input <- complete_cases_tibble[switch_indices,]$index
#   crossover_start <- input_tibble[switch_indices_input,]$positions
#   rev_input_tibble <- arrange(input_tibble, -index) %>%
#     mutate(., index = row_number())
#   complete_cases_rev_tibble <- rev_input_tibble[complete.cases(rev_input_tibble),]
#   rev_input_vec <- as.factor(complete_cases_rev_tibble[[1]])
#   rev_switch_indices <- which(rev_input_vec[-1] != rev_input_vec[-length(rev_input_vec)])
#   rev_switch_indices_input <- complete_cases_rev_tibble[rev_switch_indices,]$index
#   crossover_end <- rev(rev_input_tibble[rev_switch_indices_input,]$positions)
#   recomb_spots <- tibble(Ident = ident, Genomic_start = crossover_start, Genomic_end = crossover_end)
#   return(recomb_spots)
# }
# 
# idents_for_csv <- paste0(paste0(sampleName, "_", chrom, "_"), colnames(filled_sperm))
# recomb_spots_all <- do.call(rbind, pbmclapply(1:ncol(filled_sperm),
#                                               function(x) find_recomb_spots(filled_sperm, x, idents_for_csv, positions),
#                                               mc.cores=getOption("mc.cores", threads))) %>% 
#   right_join(., tibble(Ident = idents_for_csv), by = "Ident")

#Prepare a dataframe for the true crossover locations that can be passed to bedr
sperm_ident <- paste0(rep("sperm", num_sperm), 1:num_sperm, "_")
names(crossover_indices) <- sperm_ident
unlist_ci <- unlist(crossover_indices, use.names=TRUE) #just a list with a single row for each crossover, rowname corresponds to the sperm _1, _2, etc for sperm with more than one CO
#bedr has to have different starts and ends and cannot handle NAs
ci_df <- data.frame(chr=rep(as.character("chrT"), length(unlist_ci)), start=(unlist_ci-1), stop=(unlist_ci+1))
row.names(ci_df) <- names(unlist_ci)
truth_nona_df <- drop_na(ci_df) #this one has all valid regions for bedr
truth_onlyna_df <- ci_df[is.na(ci_df$start),] #this one is to check all the ones that don't have any recombination spots
truth_nona_df_sort <- bedr.sort.region(truth_nona_df) #these are the true recombination spots, with no NAs, lexicographically sorted


#make another bedr compatible dataframe for the predicted crossover locations
split_idents <- str_split(recomb_spots$Ident, "_", simplify=TRUE)
recomb_spots_df <- data.frame(chr=as.character(split_idents[,2]), start=recomb_spots$Genomic_start, end=recomb_spots$Genomic_end)
row.names(recomb_spots_df) <- make.names(paste0(split_idents[,3], "_"), unique=TRUE)
#bedr cannot handle NAs
pred_nona_df <- drop_na(recomb_spots_df) #this one has all valid regions for bedr
pred_onlyna_df <- recomb_spots_df[is.na(recomb_spots_df$start),] #this one is to check all the ones that don't have any recombination spots
pred_nona_df_sort <- bedr.sort.region(pred_nona_df) #these are the predicted recombination spots, with no NAs, lexicographically sorted

truth_pred_int <- bedr(input = list(a=truth_nona_df_sort, b=pred_nona_df_sort), 
                  method = "intersect",
                  params = "-loj -sorted") #those that aren't overlapping will be reported as NULL for b (. -1 -1). Those are fn
colnames(truth_pred_int) <- c("truth", "pred_chr", "pred_start", "pred_end")
pred_v <- bedr(input = list(a=pred_nona_df_sort, b=truth_nona_df_sort),
                method = "intersect",
                params = "-v -sorted") #length of this is fp

truth_pred_int_nona <- truth_pred_int[!truth_pred_int$pred_chr == ".",]
pred_idents <- paste0(truth_pred_int_nona$pred_chr, "_", truth_pred_int_nona$pred_start, "_", truth_pred_int_nona$pred_end)
tp_cons <- length(unique(pred_idents))
tp_lib <- dim(truth_pred_into_nona)[0]

fp <- length(pred_v)

fn_pt3 <- dim(truth_pred_int_nona)[0] - tp_lib
fn_pt1 <- dim(truth_pred_int[truth_pred_int$pred_chr == ".",])[0]
fn_pt2 <- length(setdiff(row.names(pred_onlyna_df), row.names(truth_onlyna_df))) #rownames in the prediction but not in the truth
fn_lib <- fn_pt1 + fn_pt2
fn_cons <- fn_pt1 + fn_pt2 + fn_pt3

tn <- dim(merge(truth_onlyna_df, pred_onlyna_df, by=0))[0]

metrics <- function(tp, fp, tn, fn){
  precision <- tp/(tp+fp)
  recall <- tp/(tp+fn)
  accuracy <- (tp + tn)/(tp + tn + fp + fn)
  f1 <- (2*precision*recall)/(precision + recall)
  specificty <- tn/(tn+fp)
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
metrics_lib <- metrics(tp_lib, fp, tn, fn_lib)
