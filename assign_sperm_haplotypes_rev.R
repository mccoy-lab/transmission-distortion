library(data.table)
library(tidyverse)
library(pbapply)
library(pbmcapply)
library(rhapsodi)

#args <- commandArgs(trailingOnly = TRUE)
input_file <- "~/mccoylab_withOthers/transmission-distortion/test_data_to_save/nc18ab_2_full_filtered_dt.csv" #args[1]
input_file_type <- tail(unlist(strsplit(input_file, "[.]")), n=1)
sampleName <- "nc18ab" #args[2]
chrom <- "2" #args[3]
outDir <-  "~/mccoyLab_withOthers/transmission-distortion/test_data_to_save/" #args[4]
threads <- 2L #s.integer(args[5])
seqError <- 0.005
avgRecomb <- 1
window_length <- 3000
smooth_imputed_genotypes <- FALSE
smooth_crossovers <- TRUE

if (input_file_type == "txt"){
  dt <- read_delim(input_file, delim = "\t", col_types = cols(.default = "d")) %>%
    as.data.frame()
}

if (input_file_type == "csv"){
  dt <- read_csv(input_file) %>% as.data.frame()
}

#Run rhapsodi
standard_input_out <- rhapsodi::standard_input(NULL, use_dt = TRUE, input_dt = dt)

# Reverse rows of input dt so that end of chromosome is first and start of chromosome is last
# since the rev function acts column-wise, we transpose the dataframe first, moving the rows to the columns, take the rev, and then re transpose moving the dt back to original
standard_input_out$dt <- t(rev(data.frame(t(standard_input_out$dt)))) %>%  as.data.frame()

#phase donor haplotypes and impute gamete genotypes 
complete_haplotypes <- rhapsodi::impute_donor_haplotypes(standard_input_out$dt, standard_input_out$positions, threads = threads, window_length = window_length)
filled_gametes <- rhapsodi::fill_gametes(standard_input_out$dt, complete_haplotypes, threads = threads, sequencing_error = seqError, avg_recomb = avgRecomb)
## Re-reverse two outputs and one input back to normal orientation
complete_haplotypes$h1 <- t(rev(data.frame(t(complete_haplotypes$h1))))
complete_haplotypes$h2 <- t(rev(data.frame(t(complete_haplotypes$h2)))) 
complete_haplotypes <- as.tibble(complete_haplotypes)
colnames(complete_haplotypes) <- c("index", "pos", "h1", "h2")
filled_gametes <- t(rev(data.frame(t(filled_gametes)))) %>% as.tibble()
standard_input_out$dt <- t(rev(data.frame(t(standard_input_out$dt)))) %>%  as.data.frame()

#finish with discovery of recombination 
rhapsodi_out <- rhapsodi::report_gametes(smooth_crossovers, smooth_imputed_genotypes, complete_haplotypes, standard_input_out$dt, filled_gametes, standard_input_out$positions, sampleName, chrom, threads=threads)

#Report Donor Haps
filename_hap <- paste0(outDir, sampleName, "_", chrom, "_parental_hap.csv")
write_csv(rhapsodi_out$donor_haps, filename_hap)

#Report filled sperm
filename_fs <- paste0(outDir, sampleName, "_", chrom, "_filled_sperm_smoothed.csv")
write_csv(filled_gametes, filename_fs)
filename_fs <- paste0(outDir, sampleName, "_", chrom, "_filled_sperm_unsmoothed.csv")
write_csv(rhapsodi_out$gamete_haps, filename_fs)

# Report recombination
filename_rs <- paste0(outDir, sampleName, "_", chrom, "_recombination_locs.csv")
write_csv(rhapsodi_out$recomb_breaks, filename_rs)

td_test <- function(sperm_matrix, row_index) {
  test_row <- sperm_matrix[row_index,]
  gt_vector <- unlist(test_row)[-1]
  one_count <- sum(gt_vector == "haplotype1", na.rm = TRUE)
  two_count <- sum(gt_vector == "haplotype2", na.rm = TRUE)
  p_value <- binom.test(c(one_count, two_count))$p.value
  return(c(p_value, one_count, two_count))
}

df_counts_pvals <- do.call(rbind, pbmclapply(1:nrow(rhapsodi_out$gamete_haps),
                                             function(x) td_test(rhapsodi_out$gamete_haps, x),
                                             mc.cores=getOption("mc.cores", threads))) %>%
  as_tibble() %>%
  add_column(standard_input_out$positions) #bind the positions vector to df_counts_pvals
colnames(df_counts_pvals) <- c("pval", "h1_count", "h2_count", "genomic_position")
filename_df <- paste0(outDir, sampleName, "_", chrom, "_pval.csv")
write_csv(df_counts_pvals, filename_df)

read_pvals <- function(sample_id, chrom, fname) {
  message(paste(sample_id, chrom))
  tryCatch(
    {
      data <- fread(paste0(fname), header = TRUE) %>%
        .[, sample_id := sample_id] %>%
        .[, chrom := chrom]
    },
    error = function(e){
      data <- data.table(pval = NA, h1_count = NA, h2_count = NA, 
                         genomic_position = NA, sample_id = NA, chrom = NA)
    })
}

pvals <- read_pvals("nc18ab", 2, filename_df)

# Take all p values below 0.01 and then sample from those above 0.01
dt_pvals_downsample <- rbind(pvals[pval < 0.01],
                             pvals[pval >= 0.01] %>%
                               .[sample(1:nrow(.), 100000),])

ggplot(data = dt_pvals_downsample, 
       aes(x = genomic_position, y = -log10(pval), color = (pmax(h1_count,h2_count) / (h2_count + h1_count))*100)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c() +
  theme_bw() +
  facet_grid(sample_id ~ factor(chrom), scales = "free_x", space = "free") + 
  theme(panel.spacing.x = unit(0, "lines"),
        axis.text.x = element_blank()) +
  xlab("Genomic Coordinate (bp)") +
  ylab(expression(-log[10](p)))
