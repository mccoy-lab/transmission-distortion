# Using genotype data from Avery Bell (in MARCC at `/work/rmccoy22/mccarroll-sperm-seq-data`) 

- First, transform SNP data into CSV file with genotype position in column 1 and each sperm cell in columns 2-n. Each cell will include 0 (reference allele), 1 (alternate allele) or NA (no read for that SNP in that sperm cell), using `gt_data_processing.R`. 
Tables output from each (22 tables from each of 20 donors) are found in `/work/scarios1/transmission-distortion/genotype_tables`

- Next, impute parental haplotypes haplotypes 
