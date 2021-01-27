# Code to plot our haplotypes with our recombination breaks and Bell's recombination breaks 
# Code to compare our recombination breaks with those of Bell (location on chromosome, resolution)
# Code to compare accuracy, precision, etc. for recombination between our method and Bell's method 
library(ggplot2)
library(cowplot)

##### 
# Format our recombination breaks, sample n sperm, graph our haplotypes with our recombination and Bell's recombination break points 

# Where `outcome` is the full matrix of imputed sperm 
outcome <- fill_gametes(dt, complete_haplotypes, sequencing_error=0.005, threads=1)

# Format the sperm matrix and randomly choose 20 sperm to consider 
outcome_new <- as.data.frame(outcome)
outcome_new <- outcome_new[,sample(1:ncol(outcome_new), 20, replace=FALSE)]

# Where positions is the SNP positions from the original data, i.e., original_dataframe[,-1]
rownames(outcome_new) <- positions

# Format sperm data 
outcome_new4 <- outcome_new %>% 
  tibble::rownames_to_column('location') %>%
  transform(location = as.numeric(location)) %>%
  pivot_longer(-'location', names_to = 'chromosome', values_to = 'haplotype')

# Grab the cell IDs of the chromosomes randomly chosen above 
inferred_cell_IDs <- unique(outcome_new4$chromosome)


# Read in Bell's recombination data and assign names to each column 
# Bell's recombination data are on Sara's local in `full_donors` with a folder for each donor and a .txt file for each fhr
# e.g., 2115 rows = 2115 recombination spots noted among the 1176 sperm cells 
bell_data <- read.table("/Users/saracarioscia/mccoy-lab/transmission-distortion/full_donors/NC10/NC10_chr1.txt", header=FALSE)

# Rename Bell's recombination data 
cell_ID <- bell_data$cell_ID
colnames(bell_data) <- c("donor", "cell_ID", "chr", "genomic_start", "genomic_end")
bell_data$Ident <- paste0(sampleName, "_", chrom, "_", cell_ID)

#Select rows (i.e., cells) from Bell data that match those in our random sample from above 
bell_cells_ofinterest2 <- subset(bell_data, cell_ID %in% inferred_cell_IDs)

# This plots the haplotype blocks, along with Bell's recombination points. It's a little hard to see because there are so many points on the x axis   
outcome_new4 %>% ggplot(aes(x=location, y = chromosome)) +
  geom_point(aes(color=haplotype)) +
  geom_point(data = bell_cells_ofinterest2, aes(x = genomic_start, y = cell_ID), shape = 18, color = "#00AFBB") + 
  geom_point(data = bell_cells_ofinterest2, aes(x = genomic_end, y = cell_ID), shape = 17, color = "#C4961A") +
  ggtitle("donor 10 chromosome 1") + 
  theme_minimal()


# Read in our exact recombination spots (generated via `find_recomb_spots` below)
# E.g., we found 2363 recombination spots (relative to Bell's 2115) and we found no recombinant sperm (i.e., all of the 1176 cells are here)
our_recomb_spots <- read.csv("/Users/saracarioscia/mccoy-lab/transmission-distortion/nc10oldoil_chr1_recombination_locs.csv", header=TRUE, sep = ",")
# Rename our identifier to match the format of the other files here 
our_recomb_spots$Ident <- str_remove_all(our_recomb_spots$Ident, "[nc10oldoil_chr1_]")
# Select cells that are in our random sample from above 
our_recomb_spots_ofinterest <- subset(our_recomb_spots, Ident %in% inferred_cell_IDs)


# These both work! First one graphs the start points and second graphs the end points of each recombination breakpoint! 
ggplot() + 
  geom_point(data = bell_cells_ofinterest2, aes(x = genomic_start, y = cell_ID), shape = 18, color = "#0072B2", alpha = 1) + 
  geom_point(data = our_recomb_spots_ofinterest, aes(x = Genomic_start, y = Ident), color = "#E69F00", alpha = 0.4) + 
  ggtitle("compare recomb spots - donor 10 chromosome 1") + 
  theme_minimal()

ggplot() + 
  geom_point(data = bell_cells_ofinterest2, aes(x = genomic_end, y = cell_ID), shape = 18, color = "#0072B2", alpha = 1) + 
  geom_point(data = our_recomb_spots_ofinterest, aes(x = Genomic_end, y = Ident), color = "#E69F00", alpha = 0.4) + 
  ggtitle("compare recomb spots - donor 10 chromosome 1") + 
  theme_minimal()



#################
# Code to calculate our recombination breakpoint exact positions (requires `outcomes` aka the full sperm matrix)

# Find our recombination breaks
find_recomb_spots <- function(input_matrix, x, identities, genomic_positions){
  ident <- identities[x]
  input_tibble <- input_matrix[, x] %>%
    mutate(., index = row_number()) %>%
    mutate(., positions = genomic_positions)
  complete_cases_tibble <- input_tibble[complete.cases(input_tibble),]
  input_vec <- as.factor(complete_cases_tibble[[1]])
  switch_indices <- which(input_vec[-1] != input_vec[-length(input_vec)])
  switch_indices_input <- complete_cases_tibble[switch_indices,]$index
  crossover_start <- input_tibble[switch_indices_input,]$positions
  rev_input_tibble <- arrange(input_tibble, -index) %>%
    mutate(., index = row_number())
  complete_cases_rev_tibble <- rev_input_tibble[complete.cases(rev_input_tibble),]
  rev_input_vec <- as.factor(complete_cases_rev_tibble[[1]])
  rev_switch_indices <- which(rev_input_vec[-1] != rev_input_vec[-length(rev_input_vec)])
  rev_switch_indices_input <- complete_cases_rev_tibble[rev_switch_indices,]$index
  crossover_end <- rev(rev_input_tibble[rev_switch_indices_input,]$positions)
  recomb_spots <- tibble(Ident = ident, Genomic_start = crossover_start, Genomic_end = crossover_end)
  return(recomb_spots)
}

threads <- 1
positions=input_datatable[, 1]
idents_for_csv <- paste0(paste0(sampleName, "_", chrom, "_"), colnames(outcome))
recomb_spots_all <- do.call(rbind, pbmclapply(1:ncol(outcome),
                                              function(x) find_recomb_spots(outcome, x, idents_for_csv, positions),
                                              mc.cores=getOption("mc.cores", threads))) %>% 
  right_join(., tibble(Ident = idents_for_csv), by = "Ident")

## Information to compare the recombination breakpoint accuracy 

# Format our data for comparison to Bell's
# This is for use in the bedr files. `recomb_spots_all` our the predicted crossover locations
# Uses Kate's script to make a bedr-compatible dataframe for the predicted crossover locations (bedr cannot handle NAs)
split_idents <- str_split(recomb_spots_all$Ident, "_", simplify=TRUE)
recomb_spots_df <- data.frame(chr=as.character(split_idents[,2]), start=recomb_spots_all$Genomic_start, end=recomb_spots_all$Genomic_end)
recomb_spots_df$chr <- paste0("chr", chr=as.character(split_idents[,2]))
row.names(recomb_spots_df) <- make.names(paste0(split_idents[,3], "_"), unique=TRUE)

pred_nona_df <- drop_na(recomb_spots_df) #this one has all valid regions for bedr
pred_onlyna_df <- recomb_spots_df[is.na(recomb_spots_df$start),] #this one is to check all the ones that don't have any recombination spots
pred_nona_df_sort <- bedr.sort.region(pred_nona_df) #these are the predicted recombination spots, with no NAs, lexicographically sorted

# Use `bell_data` from above with modifications in naming. Note `bell_data` already excludes non-recombinant sperm 
num_sperm = ncol(dt) # where dt was the original file, aka this is the number of sperm our study considered 
bell_split_idents <- str_split(bell_data$Ident, "_", simplify=TRUE)
# Format to data frame with columns `chr`, `genomic_start`, `genomic_end`
bell_recomb_spots <- data.frame(chr=as.character(bell_split_idents[,2]), start=bell_data$genomic_start, end=bell_data$genomic_end)
bell_recomb_spots$chr <- paste0("chr", chr=as.character(bell_split_idents[,2]))
row.names(bell_recomb_spots) <- make.names(paste0(bell_split_idents[,3], "_"), unique=TRUE)

bell_nona_df <- drop_na(bell_recomb_spots) #this one has all valid regions for bedr
bell_onlyna_df <- bell_recomb_spots[is.na(bell_recomb_spots$start),] #this one is to check all the ones that don't have any recombination spots
bell_nona_df_sort <- bedr.sort.region(bell_nona_df) #these are the predicted recombination spots, with no NAs, lexicographically sorted


# Plot location of recombination called by our method and Bell's 
# This will show the two plots on an x axis of the chromosome 
our_method_hist <- ggplot(data = recomb_spots_all, aes(x = Genomic_start)) + geom_histogram(bins = 150) + ggtitle("our method donor 10 chr 1")
sperm_seq_hist <- ggplot(data = bell_recomb_spots, aes(x = start)) + geom_histogram(bins = 150) + ggtitle("sperm seq")
plot_grid(our_method_hist, sperm_seq_hist, nrow = 2)


# Compare resolution of recombination offered by our method and Bell's 
our_method_resolution <- ggplot(data = recomb_spots_all, aes(x = Genomic_end - Genomic_start)) +
  geom_histogram() +
  scale_x_log10() +
  ggtitle("our method donor 10 chr 1")
sperm_seq_resolution <- ggplot(data =  bell_recomb_spots, aes(x = end - start)) +
  geom_histogram() +
  scale_x_log10() + 
  ggtitle("sperm seq")
plot_grid(our_method_resolution, sperm_seq_resolution, nrow = 2)



#########
# bedr code (run on command line) to compare the overlap in our recombination breakpoints with those from Bell 

# The `bedr` package requires access to the $PATH; this section of the script must be called from the command line. 
#Set our predictions versus Bell's 
bell_pred_intersections <- bedr(engine = "bedtools", input = list(a=bell_nona_df_sort, b=pred_nona_df_sort), 
                                method = "intersect",
                                params = "-loj -sorted") #those that aren't overlapping will be reported as NULL for b (. -1 -1). Those are fn

colnames(bell_pred_intersections) <- c("bell_chr", "bell_start", "bell_end", "pred_chr", "pred_start", "pred_end")

#Set Bell's predictions versus ours 
pred_bell_intersections <- bedr(input = list(a=pred_nona_df_sort, b=bell_nona_df_sort),
                                method = "intersect",
                                params = "-loj -sorted") #those that aren't overlapping will be reported as NULL for b (. -1 -1). Those are fn
colnames(pred_bell_intersections) <- c("pred_chr", "pred_start", "pred_end", "bell_chr", "bell_start", "bell_end")

# Calculate metrics (definitions and explanations in my 01-2021 lab notebook) 
pred_v <- pred_bell_intersections[pred_bell_intersections$bell_chr == ".",] #nrow of this is fp
fp <- nrow(pred_v)

bell_pred_int_nona <- bell_pred_intersections[!bell_pred_intersections$pred_chr == ".",]
pred_idents <- paste0(bell_pred_int_nona$pred_chr, "_", bell_pred_int_nona$pred_start, "_", bell_pred_int_nona$pred_end)
tp_cons <- length(unique(pred_idents))
tp_lib <- nrow(bell_pred_intersections)

fn_pt3 <- tp_lib - nrow(bell_pred_int_nona)
fn_pt1 <- nrow(bell_pred_intersections[bell_pred_intersections$pred_chr == ".",])
fn_pt2 <- length(setdiff(row.names(pred_onlyna_df), rownames(bell_onlyna_df)))

fn_lib <- fn_pt1 + fn_pt2
fn_cons <- fn_pt1 + fn_pt2 + fn_pt3

tn <- nrow(merge(bell_onlyna_df, pred_onlyna_df, by=0))

metrics <- function(tp, fp, tn, fn){
  precision <- tp/(tp+fp)
  recall <- tp/(tp+fn)
  accuracy <- (tp + tn)/(tp + tn + fp + fn)
  f1 <- (2*precision*recall)/(precision + recall)
  specificity <- tn/(tn+fp)
  fdr <- fp/(tp+fp)
  fpr <- fp/(tn+fp)
  metric_list <- list(precision=precision, 
                      recall=recall, 
                      accuracy=accuracy, 
                      specificity = specificity,
                      fdr = fdr,
                      fpr = fpr,
                      f1=f1)
  return (metric_list) }

metrics_cons <- metrics(tp_cons, fp, tn, fn_cons)
metrics_lib <- metrics(tp_lib, fp, tn, fn_lib)
print(metrics_cons)
print(metrics_lib)


