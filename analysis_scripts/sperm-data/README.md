# Raw genotype data from Avery Bell 
- Data for heterozygous snps of each donor (on each chromosome) are in MARCC at `/work/rmccoy22/mccarroll-sperm-seq-data`

# Subsetting data 
- To subset the raw data, we used two metrics. First, we selected only those sperm cells that were included in Bell's `autosomalploidy.txt` file (i.e,. those that were not filtered out before she investigated ploidy). Next, we selected only those that were marked euploid for the given donor and chromosome of interest. 
- The subsetted data are in MARCC at `work/scarios1/bell_filtered_data` with a directory for each donor. The naming convention is `donor_euploid_chr.txt`. 
- To generate, use the bash script `subset_bash.sh` to access each file and execute the R script `subset_data.R`. 

# Using the data to impute parental haplotypes 
- The imputation script either takes the raw data with one `dplyr` formatting step, or can take the subsetted data (only requires the notice that it's tab-delimited)