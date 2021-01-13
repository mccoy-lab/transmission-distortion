library(data.table)
library(tidyverse)
library(bedr)
library(stringr)

#args <- commandArgs(trailingOnly = TRUE)

#num_sperm <- as.integer(args[1])
num_sperm <- 1000

#num_snps <- as.integer(args[2])
num_snps <- 30000

#random_seed <- as.integer(args[3])
random_seed <- 42
set.seed(random_seed)

#coverage <- as.numeric(args[4])
coverage <- 0.001
missing_genotype_rate <- dpois(0, coverage)

#lambda <- as.numeric(args[5])
lambda <- 1
num_recomb_sites <- rpois(num_sperm, lambda)

num_genotypes <- num_sperm * num_snps
num_nas <- floor(num_genotypes * missing_genotype_rate)

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
#Got to make sure they'll be valid regions for bedr
sperm_ident <- paste0(rep("sperm", num_sperm), 1:num_sperm, "_")
names(crossover_indices) <- sperm_ident
unlist_ci <- unlist(crossover_indices, use.names=TRUE)
#have to have different starts and ends
ci_df <- data.frame(chr=rep(as.character("chrT"), length(unlist_ci)), start=(unlist_ci-1), stop=(unlist_ci+1))
row.names(ci_df) <- names(unlist_ci)
#bedr cannot handle NAs
ci_nona_df <- drop_na(ci_df) #this one has all valid regions for bedr
ci_onlyna_df <- ci_df[is.na(ci_df$start),] #this one is to check all the ones that don't have any recombination spots
#valid_ci_nona <- is.valid.region(ci_nona_df)
#sum(valid_ci_nona)
#ci_nona_df <- ci_nona_df[valid_ci_nona]

sperm_mat_with_na <- sperm_mat #copy matrix to a new matrix where we'll add in NAs so that we retain the full original knowledge
sperm_mat_with_na <- as.vector(sperm_mat_with_na)
coords_to_change <- sample(1:num_genotypes, size=num_nas, replace=FALSE)
sperm_mat_with_na[coords_to_change] <- NA
sperm_mat_with_na <- matrix(sperm_mat_with_na, ncol=num_sperm, nrow=num_snps)

#make it into a dataframe that I can give to the rest of the pipeline, so I need to have genomic positions first, column names for each sperm
sperm_na_df <- data.frame(pseudo_pos = 1:nrow(sperm_mat_with_na), sperm_mat_with_na)
sperm_full_df <- data.frame(pseudo_pos = 1:nrow(sperm_mat), sperm_mat)
#I think sperm_na_df is the df that can be passed to assign_sperm_haplotypes_rm_kw.R directly

#run the whole pipeline

#make another bedr compatible dataframe
split_idents <- str_split(recomb_spots$Ident, "_", simplify=TRUE)
recomb_spots_df <- data.frame(chr=as.character(split_idents[,2]), start=recomb_spots$Genomic_start, end=recomb_spots$Genomic_end)
row.names(recomb_spots_df) <- make.names(split_idents[,3], unique=TRUE)
#bedr cannot handle NAs
rs_nona_df <- drop_na(recomb_spots_df) #this one has all valid regions for bedr
rs_onlyna_df <- recomb_spots_df[is.na(recomb_spots_df$start),] #this one is to check all the ones that don't have any recombination spots




metrics <- function(tp, fp, tn, fn){
  precision <- tp/(tp+fp)
  recall <- tp/(tp+fn)
  accuracy <- (tp + tn)/(tp + tn + fp + fn)
  specificty <- tn/(tn+fp)
  f1 <- (2*precision*recall)/(precision + recall)
  metric_list <- list(precision=precision, recall=recall, accuracy=accuracy, specificty=specificty, f1=f1)
  return (metric_list) }
