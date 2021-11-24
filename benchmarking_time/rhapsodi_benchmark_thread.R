library(rhapsodi)

# Prevent scientific notation (for file names that have higher snp counts)
options(scipen = 999)

# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
# number of gametes
n_gam <- as.numeric(args[1])
# number of snps
n_snp <- as.numeric(args[2])
# coverage
cov <- as.numeric(args[3])
# sequencing error rate
seqerr <- as.numeric(args[4])
# average recombination rate
avgr <- as.numeric(args[5])
# random seed
rs <- as.numeric(args[6])
# Task ID 
task_ID <- as.numeric(args[7])


# Specify directory and file names and generate input file name
gen_model_dir <- "/scratch/groups/rmccoy22/kweave23/sc_transmission_distortion/generative_model_noDNM/gen_model_results_noDNM"
file_dir <- paste0("g", n_gam, "_s", n_snp, "_c", cov, "_se", seqerr, "_r", avgr)
fname_prefix <- paste0("runGen_gam_", n_gam, "_snp_", n_snp, "_cov_", cov, "_seqerr_", seqerr, "_avgr_", avgr, "_rs_", rs)
file_name <- paste0(fname_prefix, "_gametedf_na_truth_afseqednm.csv")
input_name <- paste0(gen_model_dir, '/', file_dir, '/', file_name)
input_dt <- read.csv(input_name)

rhapsodi_out <- rhapsodi::rhapsodi_autorun(NULL, input_dt = input_dt, use_dt = TRUE, threads = 8, 
                                           seqError_model = seqerr ,
                                           avg_recomb_model = avgr,
                                           mcstop = FALSE)