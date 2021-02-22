library(data.table)
library(tidyverse)
library(stringr)
library(pbapply)
library(pbmcapply)
library(HMM)
library(Hapi)
library(bedr)

args <- commandArgs(trailingOnly = TRUE)
sampleName <- "simTest"
chrom <- "chrT"
threads <- 2L #as.integer(args[1])

seqError <- 0.05
hapProb <- 1 - seqError

num_sperm <- 1000 #as.integer(args[2])
num_snps <- 30000 #as.integer(args[3])
coverage <- 0.01 #as.numeric(args[4])

num_genotypes <- num_sperm * num_snps
missing_genotype_rate <- dpois(0, coverage)
message(paste0("The missing genotype rate of this simulation is ", missing_genotype_rate))
num_nas <- floor(num_genotypes * missing_genotype_rate)

random_seed <- 27 #42 #as.integer(args[5])
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

chr <- rep(chrom, nrow(sperm_na_df))
positions <- sperm_na_df$positions
ref <- rep(0, nrow(sperm_na_df))
alt <- rep(1, nrow(sperm_na_df))
sperm_na_df_fHapi <- cbind(chr, positions, ref, alt, sperm_na_df[,-1])

###Should I filter out the non-het rows before passing to HAPI?

####INSERT HAPI STEPS
hapOutput <- hapiAutoPhase(gmt = sperm_na_df_fHapi)



