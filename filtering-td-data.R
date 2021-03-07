## Trying to filter the data to remove SNPs in the genome in a bottle consortium 
library(tidyverse)
library(data.table)

# Format blacklists, format giab, format raw data 
# Look for matches between blacklist and raw data; remove snps in the matched list 

# This file is the union of non-problematic regions of the genome 
#giab_file <- "/Users/saracarioscia/mccoy-lab/GRCh38_notinalldifficultregions.bed"
giab_file <- "/Users/saracarioscia/mccoy-lab/transmission-distortion/filtering_bell_data/giab_union_regions.bed"
giab_union <- read.table(giab_file, sep = "\t", header = FALSE) %>% as.data.frame()
#giab_union <- as.data.frame(giab_union)
colnames(giab_union) <- c("chr", "start", "end")

# to check what percent of the genome is included 
# sum(giab_union$end - giab_union$start)
# Separate giab_union into chr of interest based on input file 
#giab_chr3 <- giab_union[giab_union$chr == "chr3",] 
#giab_chr3 <- as.data.table(giab_chr3)
giab_chr3 <- giab_union[giab_union$chr == "chr3",] %>% as.data.table()


# Read in and format our raw data - hg38
nc26_chr3_file <- "/Users/saracarioscia/mccoy-lab/transmission-distortion/raw_data_for_tests/nc26abcd_euploid_3.txt"
input_dt_nc26chr3 <- read_delim(nc26_chr3_file, delim = "\t") %>% as.data.frame()

# Use positions from input file to make a bam-style file that can be used in foverlap 
positions_nc26chr3 <- fread(cmd = "cat /Users/saracarioscia/mccoy-lab/transmission-distortion/raw_data_for_tests/nc26abcd_euploid_3.txt | cut -f 1")
positions_nc26chr3[, chr := "chr3"]
positions_nc26chr3[, start := positions]
positions_nc26chr3[, end := positions]


#blacklist_file <- "/Users/saracarioscia/Downloads/ENCFF356LFX.bed.gz"
blacklist_file <- "/Users/saracarioscia/mccoy-lab/transmission-distortion/filtering_bell_data/encode_blacklist.bed.gz"
blacklist <- read_delim(blacklist_file, col_names = FALSE, delim = "\t") %>% 
  as.data.table()
#blacklist <- as.data.table(blacklist)
colnames(blacklist) <- c("chr", "start", "end")
# we want to exclude these 

# Separate blacklist into chr of interest based on input file
blacklist_chr3 <- blacklist[blacklist$chr == "chr3"]

# if overlap in blacklist, want to remove
# if overlap in giab, want to keep 

# first, take raw and find matches with blacklist
setkey(blacklist_chr3, start, end)
blacklisted_nc26chr3 <- data.table::foverlaps(positions_nc26chr3, blacklist_chr3, type="any", nomatch=NULL)
# Then, remove those rows (i.e., snps) from the positions matrix we're using to compare 
'%!in%' <- function(x,y)!('%in%'(x,y))
pos_not_blacklist_nc26chr3 <- positions_nc26chr3[positions_nc26chr3$positions %!in% blacklisted_nc26chr3$positions,]

# these are the positions that are okay to use. now `pos_not_blacklist_nc26chr3` has the info but WITHOUT the snps that were blacklisted
# now compare `pos_not_blacklist_nc26chr3` to the giab 
setkey(giab_chr3, start, end)
in_union_nc26chr3 <- data.table::foverlaps(pos_not_blacklist_nc26chr3, giab_chr3, type="any", nomatch=NULL)

snps_to_use_nc26chr3 <- input_dt_nc26chr3[input_dt_nc26chr3$positions %in% in_union_nc26chr3$positions,]

# put this into a function that we can put in the package
# so that users can have their data filtered too! if they want 


##### holding this here - it returns a formatting error in mclapply 
# Read in and format our raw data - hg38
nc26_chr3_file <- "/Users/saracarioscia/mccoy-lab/transmission-distortion/raw_data_for_tests/nc26abcd_euploid_3.txt"
dt <- read_delim(nc26_chr3_file, delim = "\t") %>% as.data.table()

dt <- dt[((rowSums(dt[, 2:ncol(dt)] == 0, na.rm = TRUE) > 0) & (rowSums(dt[, 2:ncol(dt)] == 1, na.rm = TRUE) > 0)),]

# get the positions from the inpute data
positions <- dt[, 1]


# Read in materials necessary for filtering 
# First, the encode blacklist; we want to exclude the positions here  
#blacklist_file <- "/Users/saracarioscia/Downloads/ENCFF356LFX.bed.gz"
blacklist_file <- "/Users/saracarioscia/mccoy-lab/transmission-distortion/filtering_bell_data/encode_blacklist.bed.gz"
blacklist <- read_delim(blacklist_file, col_names = FALSE, delim = "\t") %>% 
  as.data.table()
colnames(blacklist) <- c("chr", "start", "end")

# Next, the genome in a bottle (giab); we want to include only the positions here 
giab_file <- "/Users/saracarioscia/mccoy-lab/transmission-distortion/filtering_bell_data/giab_union_regions.bed"
giab_union <- read.table(giab_file, sep = "\t", header = FALSE) %>% as.data.table()
colnames(giab_union) <- c("chr", "start", "end")


filter_problem_genome <- function(dt, positions, blacklist, giab_union) {
  col_chr_name <- paste0("chr", chrom)
  # Format the raw data info in bam file format for comparison 
  new_positions <- dt[,1]
  new_positions[, chr := col_chr_name]
  new_positions[, start := positions]
  new_positions[, end := positions]
  
  # Find regions in the raw data that are in the blacklist 
  blacklist_chr <- blacklist[blacklist$chr == col_chr_name] # select chr of interest 
  setkey(blacklist_chr, start, end)
  blacklisted <- data.table::foverlaps(new_positions, blacklist_chr, type="any", nomatch=NULL)
  
  # Then, remove those rows (i.e., snps) from the list of positions we're planning to include 
  '%!in%' <- function(x,y)!('%in%'(x,y))
  pos_not_blacklist <- new_positions[new_positions$positions %!in% blacklisted$positions,]
  
  # Select entries in giab from chr of interest 
  giab_chr <- giab_union[giab_union$chr == col_chr_name,] %>% as.data.table()
  
  # Find overlaps between giab and our data table (excluding the blacklisted sites) 
  setkey(giab_chr, start, end)
  in_union <- data.table::foverlaps(pos_not_blacklist, giab_chr, type="any", nomatch=NULL)
  
  # Export only the rows from the input dataframe that are from positions in the approved positions list
  #dt <- dt[dt$positions %in% in_union$positions,]
  return(in_union)
}
in_union_list <- filter_problem_genome(dt, positions, blacklist, giab_union)
new_dt <- dt[dt$positions %in% in_union_list$positions,]












