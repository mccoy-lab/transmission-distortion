# Trying to format as a .ped file to use in PLINK 
# Each row is an individual (sperm), each column is a SNP 
# Each SNP needs to be biallelic, even if it's haploid 

# Usually a .ped file has six initial columns: 
# Family ID, Individual ID, Paternal ID, Maternal ID, Sex (1=male; 2=female; other=unknown), Phenotype

# Use flags to avoid the "obligatory" columns, which I don't have: 
# --no-fid for no Family ID; 
# --no-parents for no  paternal and maternal ID codes; 
# --no-sex (can be used with --allow-no-sex or they will all be marked as missing); 
# --no-pheno for no phenotype 
# Using all these flags together would specify the most simple kind of file: a single, unique ID code followed by all genotype data
# Command: plink --file /Users/saracarioscia/mccoy-lab/transmission-distortion/raw_data_for_tests/nc10oldoil_13_filled_sperm_converted --no-fid --no-parents --no-sex 
# --no-pheno --indep 50 5 2 --out nc10oldoil_13

library(dplyr)


args <- commandArgs(trailingOnly = TRUE)
sperm_file <- args[1]
haps_file <- args[2]
filename_ped <- args[3]
filename_map <- args[4]
chrom <- args[5]

filled_sperm <- read.csv(sperm_file)
parental_haps <- read.csv(haps_file)
colnames(parental_haps) <- c("index","pos","haplotype1","haplotype2")

# replace NA with haplotype1 to allow conversion 
filled_sperm[is.na(filled_sperm)] <- "haplotype1"

# function to replace value in sperm matrix with numeric in parental haplotypes 
re_recode_gametes <- function(dt, complete_haplotypes) {
  to_return <- data.frame(matrix(NA_real_, nrow=nrow(dt), ncol=ncol(dt)))
  for (i in 1:ncol(dt)) {
    locs_h1 <- dt[,i] == "haplotype1"
    locs_h1[which(is.na(locs_h1))] <- FALSE
    locs_h2 <- dt[,i] == "haplotype2"
    locs_h2[which(is.na(locs_h2))] <- FALSE
    to_return[locs_h1, i] <- complete_haplotypes$haplotype1[locs_h1]
    to_return[locs_h2, i] <- complete_haplotypes$haplotype2[locs_h2]
    colnames(to_return) <- colnames(dt)
  }
  return (to_return)
}
# run function 
output_recode_gametes <- re_recode_gametes(filled_sperm, parental_haps)
# add row and column names to dataframe to allow operations based on unique rows/columns 
colnames(output_recode_gametes) <- colnames(filled_sperm)
rownames(output_recode_gametes) <- 1:nrow(filled_sperm)

# fill in 0 for NA - we may want to change this in the future to more accurately reflect the recombination break   
output_recode_gametes[is.na(output_recode_gametes)] <- 0

# replace ref/alt with SNP letters - we may want to change this to reflect the true major and minor allele at these positions 
output_recode_gametes_nucleotides <- output_recode_gametes %>% 
  mutate_if(is.numeric,funs(ifelse(. == 1, 'AA', 'GG')))


# remove names and transpose into dataframe. Plink wants each row to be an individual and each column to be a genomic position 
names(output_recode_gametes_nucleotides) <- NULL
output_recode_gametes_nucleotides_plink = t(output_recode_gametes_nucleotides) %>% as.data.frame()

# add column for the individual ID - each individual must have a unique ID 
output_recode_gametes_nucleotides_plink <- cbind(1:nrow(output_recode_gametes_nucleotides_plink),output_recode_gametes_nucleotides_plink)
# remove row and column names again for writing to file 
names(output_recode_gametes_nucleotides_plink) <- NULL


# Write converted SNP file to ped file 
write.table(output_recode_gametes_nucleotides_plink, filename_ped, col.names = FALSE, row.names = FALSE, quote = FALSE)


# Make map file - chr, SNP ID, centimorgans, basepair
pos <- parental_haps$pos
cm_pos <- cbind(a=0, b=pos)
ID_cm_pos <- cbind(a=1:nrow(filled_sperm), b=cm_pos)
map_full <- cbind(a=chrom, b=ID_cm_pos)
colnames(map_full) <- NULL
# write to map file 
write.table(map_full, filename_map, col.names = FALSE, row.names = FALSE)

