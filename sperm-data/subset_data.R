## This script is called by subset_bash.sh. It acts on each file from the McCarroll lab to generate subsetted files containing 
## only the cells that a) were in Bell's autosomal ploidy analysis `autosomalploidy.txt` and b) are euploid. 
## It outputs this subsetted data to a new file with naming convention `donor_euploid_chr.txt`. 
## The subsetted data are on MARCC at `/work/scarios1/bell_filtered_data` with one folder per donor.  

#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

library(tidyverse)
library(dplyr)

donor <- args[3]
donor <- readr::parse_number(donor)
chr <- args[4]
chr <- readr::parse_number(chr)


aneuploid_info <- read.table("autosomalploidy.txt", header = TRUE, sep = "\t")
colnames(aneuploid_info) <- c("donor", "cell", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
"18", "19", "20", "21", "22")
colnames(aneuploid_info) <- as.character(colnames(aneuploid_info))
aneuploid_info$donor <- gsub("^.{0,2}", "", aneuploid_info$donor)

donor_of_interest <- aneuploid_info[which(aneuploid_info$donor == donor),]

input_file <- args[1]

full_data <- read_delim(input_file, delim = "\t") %>%
  pivot_wider(., names_from = "cell", values_from = "gt") %>%
  arrange(., pos) %>%
  as.data.frame()


# Remove positions from data because they isn't relevant for this subsetting
positions <- full_data[,1]
full_data <- full_data[,-1]

# Subset the full data to include only the cells that made it to Bell's aneuploidy analysis
cells_in_full_data <- colnames(full_data)
subsetted_data <- full_data[,which(cells_in_full_data %in% donor_of_interest$cell)]
cells_in_subset <- colnames(subsetted_data)

# Subset the subsetted data to include only the cells that were euploid in Bell's analysis
euploid_cells <- which(donor_of_interest[, which(colnames(donor_of_interest) == chr)] == "_")
euploid_cell_IDs <- donor_of_interest[euploid_cells,2]

euploid_data <- subsetted_data[,which(cells_in_subset %in% euploid_cell_IDs)]
euploid_data <- cbind(positions, euploid_data)


# Write to file
write.table(euploid_data, file=args[2], sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)