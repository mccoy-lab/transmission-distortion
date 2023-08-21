# Using an individual from the 1000 genomes project, individual # NA0606985
# process vcf file (link in shell-scripts) using /Users/saracarioscia/mccoy-lab/transmission-distortion/shell-scripts/vcf-sperm-sim.sh

# read in hap file into a dataframe: two columns (one for each parental chromosome) x many rows (# of SNPs) 
# 0 means haplotype had the same allele as the reference, 1 means the haplotype had an alternative allele
# the two parental haplotypes may be homozygous for 1 at a position (i.e., both parents had an alternative allele)
# there are no sites where the parental haplotypes are homozygous for 0 (they would not be SNPs for this individual and are thus not incl. in the .hap file)
hap_NA06985 = read.table("/Users/saracarioscia/mccoy-lab/transmission-distortion/vcf_data/NA06985_donor.hap", header = FALSE, sep = " ") # loads a dataframe

# remove the homozygous sites (where both parental haplotypes have 1) - since we are interested in differential transmission of an allele to the sperm
hap_NA06985 <- hap_NA06985[(hap_NA06985$V1 + hap_NA06985$V2) == 1,]
# For individual NA06985, at heterozygous sites (i.e., having removed sites where both are 1) hap1 has 15,182 alleles with 1 and hap2 has 16,552 alleles with 1


# we simulate a matrix with n_sperm columns and n_SNPs rows
rows <- nrow(hap_NA06985)
M_NA06985 <- matrix(ncol = 10000, nrow = rows)

# using x, assign recombination breakpoints randomly in the parental haplotypes to generate each sperm
recomb_spot = sample(1:(rows-1), 10000, replace=T)

# for each of the sperm, choosing which parental haplotype to use for the top of each sperm (hap_first) and then use the other haplotype as the "bottom"
for(i in 1:10000) {
  haps <- c(1, 2)
  hap_first = sample(haps, size = 1)
  if (hap_first < 2) {
    hap_second <- 2 
  } else {
    hap_second <- 1
  }
  # use hap_first for the top of each sperm, from row 1 until row recomb_spot[i] 
  M_NA06985[1:recomb_spot[i], i] = hap_NA06985[1:recomb_spot[i], hap_first]
  M_NA06985[(recomb_spot[i] + 1):rows, i] = hap_NA06985[(recomb_spot[i] + 1):rows, hap_second]
}

#deprecated: now, remove homozygous 1 sites at the beginning
#het_index = which(hap_NA06985$V1 + hap_NA06985$V2 == 1) # deprecated: this was how to subset the matrix to only incl het
#p_vals <- apply(M_NA06985[het_index,], MARGIN = 1, FUN = function(x) binom.test(table(x))$p.value)


p_vals <- apply(M_NA06985, MARGIN = 1, FUN = function(x) binom.test(table(x))$p.value)
hist(p_vals, breaks = 100)
plot(-log10(p_vals), ylim = c(0, 8))

# choose a SNP to experience biased transmission
snp_of_interest <- 24377
# grab rows where the sperm have a 0 or 1 at the snp of interest
soi_zero <- which(M_NA06985[snp_of_interest,] == 0)
soi_one <- which(M_NA06985[snp_of_interest,] == 1)

# make biased matrix based on a set ratio of zeros and ones for that SNP (here, 300 and 200)
M_NA06985_biased <- M_NA06985[, c(sample(soi_zero, 300, replace = FALSE), sample(soi_one, 200, replace = FALSE))]
#deprecated: now, remove homozygous sites immediately:  p_vals_biased <- apply(M_NA06985_biased[het_index,], MARGIN = 1, FUN = function(x) binom.test(table(x))$p.value)
p_vals_biased <- apply(M_NA06985_biased, MARGIN = 1, FUN = function(x) binom.test(table(x))$p.value)
hist(p_vals_biased, breaks = 100)
plot(-log10(p_vals_biased), ylim = c(0, 8)) 
# output named sperm-sim-log-pvalues-biased.png
