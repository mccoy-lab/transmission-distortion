library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
#num_sperm <- as.integer(args[1])
num_sperm <- 15
#num_snps <- as.integer(args[2])
num_snps <- 30000
#coverage <- as.numeric(args[3])
coverage <- 0.001

#start with 2 parental chromosomes of heterozygous sites
#first parental chromosome
hap1 <- data.frame(V1 = sample(c(0, 1), size = num_snps, replace = TRUE))
hap2 <- abs(1-hap1)
parental_haps <- data.frame(cbind(hap1, hap2))
colnames(parental_haps) <- c("Parental1", "Parental2")

indices <- 1:num_snps
pseudo_positions <- (indices+40)*1000 #expect a SNP every 1000 bp or so 


sperm_mat <- matrix(ncol=num_sperm, nrow=num_snps)
num_recomb_sites <- sample(c(0,1,2), size=num_sperm, replace = TRUE, prob=c(0.05, 0.75, 0.20))
#need to store the recombination sites so that we can compare later results with ground truth
for (i in 1:num_sperm) {
  num_rs <- num_recomb_sites[i]
  rs_sites <- sample(1:(num_snps-1), size=num_rs, replace=FALSE, prob=rep(1/(num_snps-1), num_snps-1))
  if (num_rs > 0){
    hap_first <- sample(c(1,2), size=1, prob=rep(0.5,2))
    if (hap_first == 1){
      hap_second <- 2
    } else {hap_second <- 1}
    if (num_rs > 1) { #two recombination sites, so we need hap_first, hap_second, and hap_first again to build the sperm matrix column
      sperm_mat[1:rs_sites[1], i] = parental_haps[1: rs_sites[1], hap_first]
      sperm_mat[(rs_sites[1]+1) : (rs_sites[2]-1), i] = parental_haps[(rs_sites[1]+1) : (rs_sites[2]-1), hap_second]
      sperm_mat[rs_sites[2]:num_snps, i] = parental_haps[rs_sites[2]:num_snps, hap_first]
    } else { #only one recombination site, so we only need hap_first and hap_second to build the sperm matrix column
      sperm_mat[1:rs_sites[1], i] = parental_haps[1: rs_sites[1], hap_first]
      sperm_mat[(rs_sites[1]+1): num_snps, i] = parental_haps[(rs_sites[1]+1): num_snps, hap_second]
    }
    
  } else { #no recombination sites, so we pick a parental chromosome and completely copy to build the sperm matrix column
    hap_single <- sample(c(1,2), size=1, prob=rep(0.5,2))
    sperm_mat[,i] = parental_haps[, hap_single] 
  }
}

sperm_mat_with_na <- sperm_mat #copy matrix to a new matrix where we'll add in NAs so that we retain the full original knowledge