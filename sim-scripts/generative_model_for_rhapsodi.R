library(data.table)
library(tidyverse)
library(stringr)
library(pbapply)
library(pbmcapply)

args <- commandArgs(trailingOnly = TRUE)
outDir <- args[1]
threads <- 8
# seqError <- 0.005
# hapProb <- 1 - seqError
# hmm_avg_recomb <- 1
num_gametes <- as.integer(args[2]) #3, 15, 50, 150, 500, 1000, 2500, 5000
num_snps <- as.integer(args[3]) #5000, 30000, 100000
coverage <- as.numeric(args[4]) #0.001, 0.01, 0.1, 0.223, 0.357, 0.511, 0.693, 1.204, 2.303
message(paste0("The coverage of this simulation is: ", coverage))
num_genotypes <- num_gametes * num_snps
missing_genotype_rate <- dpois(0, coverage)
message(paste0("The missing genotype rate of this simulation is ", missing_genotype_rate))
num_nas <- floor(num_genotypes * missing_genotype_rate)
window_length <- 3000
overlap_denom <- 2
random_seed <- as.integer(args[5]) #42, 347, 1848
set.seed(random_seed)
recomb_lambda <- as.numeric(args[6]) #0.6, 1, 3
stopifnot(recomb_lambda > 0)
num_recomb_sites <- rpois(num_gametes, recomb_lambda)
message(paste0("Total number of recombination spots across gametes: ", sum(num_recomb_sites)))
add_seq_error <- TRUE
seqError_add <- as.numeric(args[7]) #0.001, 0.005, 0.05
add_de_novo_mut <- FALSE
de_novo_lambda <- 5
de_novo_alpha <- 7.5
de_novo_beta <- 10
stopifnot(de_novo_beta > 0)
write_out_plot <- TRUE
sampleName <- paste("runGen", "gam", num_gametes, "snp",  num_snps, "cov", coverage, "seqerr", seqError_add, "avgr", recomb_lambda, "rs", random_seed, sep="_")

###Genenrate 2 donor chromosomes of heterozygous sites
hap1 <- data.frame(V1 = sample(c(0, 1), size = num_snps, replace = TRUE)) #simulate first donor chromosome
hap2 <- abs(1-hap1) #switch bits to construct second donor
donor_haps <- data.frame(cbind(hap1, hap2))
colnames(donor_haps) <- c("donor1", "donor2")

###Generate gametes from these 2 donor chromosomes
generate_gam <- function(donor_haplotypes, n_crossovers){
  #this function generates a single gamete given donor haplotypes and the number of crossovers that should occurs for a given gamete
  # -- pick which donor haplotype corresponds to the beginning of the gamete's chromosome
  # -- assign the 0 & 1 genotypes of this first donor haplotype to the whole gamete
  # -- assign the 1 or 2 of reprenesting which haplotype to the whole gamete
  #     -- re-code this to 0 or 1 (just subtract 1)
  #     -- at the end re-recode this back to 1's and 2's (just add 1)
  # -- for each crossover index, i, between index i-1 and index i, switch the bits for the rest of the chromosome genotypes to the opposite donor haplotype  or genotype compared to the one just previous (abs(value-1))
  #input  -- donor_haplotypes: a dataframe with the two donor haplotypes
  #          n_crossovers: integer, the number of crossovers for this given gamete
  #output -- a list with 3 elements
  #             --crossover_indices, a vector of integers, such that each integer, wlog i, represents an index location of a crossover exchange point, such that the crossover occurred somewhere between SNP index locations i-1 and i
  #             --recombined_hap, a vector of 0's and 1's, representing the gamete's fully known generated genotype
  #             --recombined_hap_index, a vector of 1's and 2's, representing which donor haplotype that location's genotype originated from
  init_hap_index <- sample(1:2, 1) #which donor haplotype corresponds to the beginning of the gamete's chromosome
  init_hap <- donor_haplotypes[,init_hap_index]
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

sim_gam <- lapply(1:num_gametes, function(x) generate_gam(donor_haps, num_recomb_sites[x]))
gam_mat <- sapply(sim_gam, "[[", 2)
gam_haps <- sapply(sim_gam, "[[", 3)

crossover_indices <- sapply(sim_gam, "[[", 1)
names(crossover_indices) <- paste0(rep("gam", num_gametes), 1:num_gametes, "_")
unlist_ci <- unlist(crossover_indices, use.names=TRUE)
tci_dt <- data.table(gam=str_split(names(unlist_ci), "_", simplify=TRUE)[,1], start =(unlist_ci-1), end=(unlist_ci))

gam_full_df <- data.frame(pseudo_pos = 1:nrow(gam_mat), gam_mat)
colnames(gam_full_df) <- c("positions", paste0("gam", 1:num_gametes, "_"))

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

###Sparsify the knowledge of the gametes
add_to_na_flatten <- function(to_add_from, num_nas, num_gametes, num_snps){
  to_return <- rep(NA, (num_snps*num_gametes))
  coords_to_keep_genotype <- sample(1:(num_snps*num_gametes), size=((num_snps*num_gametes)-num_nas), replace = FALSE)
  to_return[coords_to_keep_genotype] <- as.vector(to_add_from)[coords_to_keep_genotype]
  to_return <- matrix(to_return, ncol=num_gametes, nrow=num_snps)
  return (to_return)
}

add_na_flatten <- function(to_change, num_nas, num_gametes, num_snps){
  to_change <- as.vector(to_change)
  coords_to_change <- sample(1:(num_snps*num_gametes), size=num_nas, replace = FALSE)
  to_change[coords_to_change] <- NA
  to_return <- matrix(to_change, ncol=num_gametes, nrow=num_snps)
  return (to_return)
}

if (missing_genotype_rate > 0.5){
  gam_mat_with_na <- add_to_na_flatten(gam_mat, num_nas, num_gametes, num_snps)
} else if (missing_genotype_rate <= 0.5){
  gam_mat_with_na <- add_na_flatten(gam_mat, num_nas, num_gametes, num_snps)
}
gam_na_df <- data.frame(pseudo_pos = 1:nrow(gam_mat_with_na), gam_mat_with_na)
colnames(gam_na_df) <- c("positions", paste0("gam", 1:num_gametes, "_"))

if (write_out_plot){ #full truths, prior to filtering, DNMs, and sequencing error
  filename_dh <- paste0(outDir, sampleName, "_donorHaps_truth_full_ptfseqednm.csv")
  write_csv(donor_haps, filename_dh)

  filename_ci <- paste0(outDir, sampleName, "_crossoverIndices_truth_ptfseqednm.csv")
  write_csv(tci_dt, filename_ci)

  filename_gh <- paste0(outDir, sampleName, "_gameteDonorHap_byPosition_truth_ptfseqednm.csv")
  write_csv(as.data.frame(gam_haps), filename_gh)

  filename_gm <- paste0(outDir, sampleName, "_gametedf_full_truth_ptfseqednm.csv")
  write_csv(gam_full_df, filename_gm)

  filename_gn <- paste0(outDir, sampleName, "_gametedf_na_truth_ptfseqednm.csv")
  write_csv(gam_na_df, filename_gn)

  filename_td <- paste0(outDir, sampleName, "_pval_sim_truth_full_ptfseqednm.csv")
  write_csv(df_counts_pvals_truth, filename_td)

}

###Include de novo mutations
if (add_de_novo_mut){
  new_rows <- c()
  num_dnm <- rpois(1, de_novo_lambda) + 1 #make sure this is greater than 0 by adding one
  message(paste0("Number of de novo mutations: ", num_dnm))
  num_gametes_affected_per_dnm <- ceiling(rgamma(num_dnm, de_novo_alpha, scale=de_novo_beta))
  donors_with_dnm <- sample(1:2, num_dnm, replace=TRUE)
  for (i in 1:num_dnm){ #for every donor dnm
    donor_with_dnm <- donors_with_dnm[i]
    row_loc <- sample(1:num_snps, 1) #pick random row which we'll add this new de novo mutation after
    gam_can_be_affected <- which(gam_haps[row_loc, ]==donor_with_dnm)
    num_gametes_affected <- min(num_gametes_affected_per_dnm[i], length(gam_can_be_affected))
    message(paste0("Number of gametes affected for de novo mutation ", i, ": ", num_gametes_affected))
    where_locs_greater <- which(new_rows >= (row_loc+1))
    new_rows[where_locs_greater] <- new_rows[where_locs_greater] + 1
    new_rows <- c(new_rows, (row_loc+1))
    message(paste0("new row location: ", row_loc+1))
    message(paste0("new rows vector: ", new_rows, collapse= " ; "))
    donor_haps <- rbind(donor_haps[1:row_loc,], rep(0, 2), donor_haps[(row_loc+1):num_snps, ])
    donor_haps[(row_loc+1), donor_with_dnm] <- 1
    affected_gam <- sort(sample(gam_can_be_affected, num_gametes_affected))
    #message(paste0("affected gametes: ", affected_gam, collapse= " ; "))
    new_row_haps <- rep(abs((donor_with_dnm - 1)-1)+1, num_gametes)
    new_row_haps[gam_can_be_affected] <- donor_with_dnm
    new_row_vals <- rep(0, num_gametes)
    new_row_vals[affected_gam] <- 1
    #add sparsity in
    num_nas_to_add <- floor(num_gametes * missing_genotype_rate)
    if (missing_genotype_rate <= 0.5){
      change_indices <- sample(1:num_gametes, num_nas_to_add)
      new_row_vals[change_indices] <- NA
    } else if (missing_genotype_rate > 0.5){
      sparse_row <- rep(NA, num_gametes)
      keep_indices <- sample(1:num_gametes, num_gametes - num_nas_to_add)
      sparse_row[keep_indices] <- new_row_vals[keep_indices]
      new_row_vals <- sparse_row
    }
    #add in new SNP line to gamete data, automatically changing pseudo positions too for when we later make this a dataframe
    gam_mat <- rbind(gam_mat[1:row_loc, ], new_row_vals, gam_mat[(row_loc+1):num_snps, ])
    gam_mat_with_na <- rbind(gam_mat_with_na[1:row_loc, ], new_row_vals, gam_mat_with_na[(row_loc+1):num_snps, ])
    #add in new haps line
    gam_haps <- rbind(gam_haps[1:row_loc, ], new_row_haps, gam_haps[(row_loc+1):num_snps, ])
    #add in new SNP to number of SNPs
    num_snps <- num_snps + 1 #increase num_snps by one
    #Adjust recombination breakpoints
    bool_adj <- unlist_ci >= row_loc+1
    bool_adj[is.na(bool_adj)] <- FALSE
    unlist_ci[bool_adj] <- unlist_ci[bool_adj] + 1
  }
  tci_dt <- data.table(gam=str_split(names(unlist_ci), "_", simplify=TRUE)[,1], start =(unlist_ci-1), end=(unlist_ci))

  df_counts_pvals_truth <- do.call(rbind, pbmclapply(1:nrow(gam_haps),
                                                     function(x) td_test_truth(gam_haps, x),
                                                     mc.cores=getOption("mc.cores", threads))) %>%
    as_tibble() %>%
    add_column(1:nrow(gam_haps)) #bind the positions vector to df_counts_pvals
  colnames(df_counts_pvals_truth) <- c("pval", "h1_count", "h2_count", "genomic_position")

  gam_na_df <- data.frame(pseudo_pos = 1:nrow(gam_mat_with_na), gam_mat_with_na)
  colnames(gam_na_df) <- c("positions", paste0("gam", 1:num_gametes, "_"))

  gam_full_df <- data.frame(pseudo_pos = 1:nrow(gam_mat), gam_mat)
  colnames(gam_full_df) <- c("positions", paste0("gam", 1:num_gametes, "_"))

  if (write_out_plot){
    filename_df <- paste0(outDir, sampleName, "_pval_sim_truth_full_ptfseqe_adnm.csv")
    write_csv(df_counts_pvals_truth, filename_df)

    filename_dh <- paste0(outDir, sampleName, "_donorHaps_truth_full_ptfseqe_adnm.csv")
    write_csv(donor_haps, filename_dh)

    filename_ci <- paste0(outDir, sampleName, "_crossoverIndices_truth_ptfseqe_adnm.csv")
    write_csv(tci_dt, filename_ci)

    filename_gh <- paste0(outDir, sampleName, "_gameteDonorHap_byPosition_truth_ptfseqe_adnm.csv")
    write_csv(as.data.frame(gam_haps), filename_gh)

    filename_gm <- paste0(outDir, sampleName, "_gametedf_full_truth_ptfseqe_adnm.csv")
    write_csv(gam_full_df, filename_gm)

    filename_gn <- paste0(outDir, sampleName, "_gametedf_na_truth_ptfseqe_adnm.csv")
    write_csv(gam_na_df, filename_gn)
  }
}

#Should add something to increase number of NAs after DNMs

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

if (write_out_plot){
  filename_gn <- paste0(outDir, sampleName, "_gametedf_na_truth_ptf_aseqednm.csv")
  write_csv(gam_na_df, filename_gn)
}

###Run the filtering that the pipeline would do
##Ensure that each row has at least one 0 and one 1
keep_bool <- unname((rowSums(gam_na_df[, 2:ncol(gam_na_df)] == 0, na.rm = TRUE) > 0) & (rowSums(gam_na_df[, 2:ncol(gam_na_df)] == 1, na.rm = TRUE) > 0))
gam_na_df <- gam_na_df[keep_bool,]
gam_full_df <- gam_full_df[keep_bool,]
donor_haps <- donor_haps[keep_bool,]

# ci_post_filter <- ci_pre_filter[keep_bool,]
num_snps <- sum(keep_bool)
message(paste0("new number of snps: ", num_snps))
if (add_de_novo_mut) {
  `%notin%` <- Negate(`%in%`)
  for (i in 1:length(new_rows)){
    message(paste0("dnm ", i, " is filtered out: ", new_rows[i] %notin% gam_na_df[,1]))
  }
}

if (write_out_plot){
  filename_df <- paste0(outDir, sampleName, "_pval_sim_truth_afseqednm.csv")
  write_csv(df_counts_pvals_truth[keep_bool,], filename_df)

  filename_dh <- paste0(outDir, sampleName, "_donorHaps_truth_afseqednm.csv")
  write_csv(donor_haps, filename_dh)

  filename_gn <- paste0(outDir, sampleName, "_gametedf_na_truth_afseqednm.csv")
  write_csv(gam_na_df, filename_gn)

  filename_gm <- paste0(outDir, sampleName, "_gametedf_full_truth_afseqednm.csv")
  write_csv(gam_full_df, filename_gm)

  filename_gh <- paste0(outDir, sampleName, "_gameteDonorHap_byPosition_truth_ptfseqe_adnm.csv")
  write_csv(as.data.frame(gam_haps[keep_bool,]), filename_gh)

  ###Plot the distribution of sparsity by row/SNP in the simulated gamete data
  real_reads <- rowSums(!is.na(gam_na_df))
  filename = paste0(outDir, sampleName, "_simulated_notna_bysnp_filt.pdf")
  pdf(file=filename)
  hist(real_reads, breaks=150, xlab="not NA by SNP", main="Distribution of simulated data not NA by SNP")
  dev.off()

  ###Plot the distribution of sparsity by column/gamete in the simulated gamete data
  real_reads2 <- colSums(!is.na(gam_na_df))
  filename = paste0(outDir, sampleName, "_simulated_notna_bygam_filt.pdf")
  pdf(file=filename)
  hist(real_reads2, breaks=100, xlab="not NA by gam", main="Distribution of simulated data not NA by gam")
  dev.off()
}
