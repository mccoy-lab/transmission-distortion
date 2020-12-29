library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
#num_sperm <- as.integer(args[1])
num_sperm <- 1000
#num_snps <- as.integer(args[2])
num_snps <- 30000
#coverage <- as.numeric(args[3])
coverage <- 0.001
missing_genotype_rate <- dpois(0, coverage)
num_genotypes <- num_sperm * num_snps
num_nas <- floor(num_genotypes * missing_genotype_rate)

#start with 2 parental chromosomes of heterozygous sites
#first parental chromosome
hap1 <- data.frame(V1 = sample(c(0, 1), size = num_snps, replace = TRUE))
hap2 <- abs(1-hap1)
parental_haps <- data.frame(cbind(hap1, hap2))
colnames(parental_haps) <- c("Parental1", "Parental2")

indices <- 1:num_snps
pseudo_positions <- (indices+40)*1000 #expect a SNP every 1000 bp or so 

#store full genotypes in a matrix
sperm_mat <- matrix(ncol=num_sperm, nrow=num_snps)

#need to store the recombination sites so that we can compare later results with ground truth
recomb_df <- matrix(nrow=num_sperm, ncol=3, data=NA) %>% data.frame()
colnames(recomb_df) <- c('Recombination_spot_pseudo_position_1', 'Recombination_spot_pseudo_position_2', 'hap_first')

num_recomb_sites <- sample(c(0,1,2), size=num_sperm, replace = TRUE, prob=c(0.05, 0.75, 0.20))

for (i in 1:num_sperm) {
  num_rs <- num_recomb_sites[i]
  rs_sites <- sort(sample(1:(num_snps-1), size=num_rs, replace=FALSE, prob=rep(1/(num_snps-1), num_snps-1)), decreasing=FALSE)
  
  if (num_rs > 0){
    hap_first <- sample(c(1,2), size=1, prob=rep(0.5,2))
    if (hap_first == 1){
      hap_second <- 2
    } else {hap_second <- 1}
    if (num_rs > 1) { #two recombination sites, so we need hap_first, hap_second, and hap_first again to build the sperm matrix column
      recomb_df[i,1:2] <- pseudo_positions[rs_sites]
      sperm_mat[1:rs_sites[1], i] = parental_haps[1: rs_sites[1], hap_first]
      sperm_mat[(rs_sites[1]+1) : (rs_sites[2]-1), i] = parental_haps[(rs_sites[1]+1) : (rs_sites[2]-1), hap_second]
      sperm_mat[rs_sites[2]:num_snps, i] = parental_haps[rs_sites[2]:num_snps, hap_first]
    } else { #only one recombination site, so we only need hap_first and hap_second to build the sperm matrix column
      recomb_df[i, 1] <- pseudo_positions[rs_sites]
      sperm_mat[1:rs_sites[1], i] = parental_haps[1: rs_sites[1], hap_first]
      sperm_mat[(rs_sites[1]+1): num_snps, i] = parental_haps[(rs_sites[1]+1): num_snps, hap_second]
    }
    recomb_df[i, 3] <- hap_first
  } else { #no recombination sites, so we pick a parental chromosome and completely copy to build the sperm matrix column
    hap_single <- sample(c(1,2), size=1, prob=rep(0.5,2))
    sperm_mat[,i] = parental_haps[, hap_single]
    recomb_df[i, 3] <- hap_single 
  }
}

sperm_mat_with_na <- sperm_mat #copy matrix to a new matrix where we'll add in NAs so that we retain the full original knowledge
sperm_mat_with_na <- as.vector(sperm_mat_with_na)
coords_to_change <- sample(1:num_genotypes, size=num_nas, replace=FALSE)
sperm_mat_with_na[coords_to_change] <- NA
sperm_mat_with_na <- matrix(sperm_mat_with_na, ncol=num_sperm, nrow=num_snps)

#make it into a dataframe that I can give to the rest of the pipeline, so I need to have genomic positions first, column names for each sperm
colnames_df <- c("Pseudo_genomic_positions", paste0("test_sperm_", 1:num_sperm))
sperm_na_df <- cbind(pseudo_positions, sperm_mat_with_na) %>% as.data.frame()
colnames(sperm_na_df) <- colnames_df
sperm_full_df <- cbind(pseudo_positions, sperm_mat) %>% as.data.frame()
colnames(sperm_full_df) <- colnames_df

#I think sperm_na_df is the df that can be passed to assign_sperm_haplotypes_rm_kw.R directly