library(HMM)
library(Hapi)
library(knitr)
library(rhapsodi)
library(tidyverse)

# Prevent scientific notation (for file names that have higher snp counts)
options(scipen = 999)

# Create a log that will store messages about script progress across all files
log <- data.frame(n_gam=c(), n_snp=c(), n_het_snp=c(), cov=c(), seqerr=c(), avgr=c(), rs=c(),
                  gmt_import=c(), filter_error=c(), frame_select=c(),impute_and_phase=c(),final_draft=c(),consenus=c(), consensus_end=c(),CO=c())

# Read in command line arguments
args = commandArgs(trailingOnly=TRUE)
# number of gametes
n_gam = as.numeric(args[1])
# number of snps
n_snp = as.numeric(args[2])
# coverage
cov = as.numeric(args[3])
# sequencing error rate
seqerr = as.numeric(args[4])
# average recombination rate
avgr = as.numeric(args[5])
# random seed
rs = as.numeric(args[6])
# Task ID 
task_ID = as.numeric(args[7])


# Specify directory and file names and generate input file name
gen_model_dir <- "~/work/kweave23/sc_transmission_distortion/generative_model_noDNM/gen_model_results_noDNM"
file_dir <- paste0("g", n_gam, "_s", n_snp, "_c", cov, "_se", seqerr, "_r", avgr)
fname_prefix <- paste0("runGen_gam_", n_gam, "_snp_", n_snp, "_cov_", cov, "_seqerr_", seqerr, "_avgr_", avgr, "_rs_", rs)
file_name <- paste0(fname_prefix, "_gametedf_na_truth_afseqednm.csv")
input_name <- paste0(gen_model_dir, '/', file_dir, '/', file_name)
print(fname_prefix)

# Specify output names 
log_name <- paste0('~/work/abortvi2/rhapsodi/command_line_results/logs/log_', task_ID, '.csv')
SAGI_name <- paste0('~/work/abortvi2/rhapsodi/command_line_results/out/SAGI/SAGI_test', task_ID, '.csv')
SAP_name <- paste0('~/work/abortvi2/rhapsodi/command_line_results/out/SAP/SAP_test_', task_ID, '.csv')
SAR_name <- paste0('~/work/abortvi2/rhapsodi/command_line_results/out/SAR/SAR_test_', task_ID, '.csv')

# Read in data. Add chr ref and alt columns. Reorder columns to fit hapi's intended format.
gmt <- read.csv(input_name)
 
# Break if no observations
if (nrow(gmt) == 0) {
  log <- rbind(log, list(n_gam=n_gam, n_snp=n_snp, n_het_snp=0, cov=cov, seqerr=seqerr, avgr=avgr, rs=rs,
                         gmt_import=0,
                         filter_error=NA,
                         frame_select=NA,
                         impute_and_phase=NA,
                         final_draft=NA,
                         consenus=NA,
                         consensus_end=NA,
                         CO=NA), stringsAsFactors=FALSE)
  write.table(log, log_name , append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = TRUE)
  print(paste0(file_name, " has no hetSNPs."))
  q()
}

rownames(gmt) <- gmt$position

gmt['ref']=0
gmt['alt']=1
gmt['chr']=1

gmt <- gmt[c(1, ncol(gmt), ncol(gmt)-2, ncol(gmt)-1, 2:(ncol(gmt)-3))]

# Split gmt into hetDA and gmtDA
hetDa <- gmt[,1:4]
ref <- hetDa$ref
alt <- hetDa$alt

gmtDa <- gmt[,-(1:4)]
# Define HMM
co_rate = avgr/(nrow(gmtDa))
hmm = initHMM(States=c("S","D"), Symbols=c("s","d"),
              transProbs=matrix(c(1-co_rate,co_rate,co_rate,1-co_rate),2),
              emissionProbs=matrix(c(1-seqerr,seqerr,seqerr,1-seqerr),2),
              startProbs = c(0.5,0.5))

# Implement HMM; Filter Error
gmtDA <- tryCatch(hapiFilterError(gmt = gmtDa, hmm = hmm), error= function(e) {return(2)}  )
if (length(gmtDA) == 1){
  if (gmtDA==2) {
    log <- rbind(log, list(n_gam=n_gam, n_snp=n_snp, n_het_snp=nrow(gmtDa), cov=cov, seqerr=seqerr, avgr=avgr, rs=rs,
                           gmt_import=1,
                           filter_error=0,
                           frame_select=NA,
                           impute_and_phase=NA,
                           final_draft=NA,
                           consenus=NA,
                           consensus_end=NA,
                           CO=NA),
                 stringsAsFactors=FALSE)
    write.table(log, log_name , append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = TRUE)
    print(paste0(file_name, " could not filter error"))
    q()
  }
}

# Frame Selection
gmtFrame <- hapiFrameSelection(gmt = gmtDA, n = 3)

if (nrow(gmtFrame)==0){
  log <- rbind(log, list(n_gam=n_gam, n_snp=n_snp, n_het_snp=nrow(gmtDa), cov=cov, seqerr=seqerr, avgr=avgr, rs=rs,
                         gmt_import=1,
                         filter_error=1,
                         frame_select=0,
                         impute_and_phase=NA,
                         final_draft=NA,
                         consenus=NA,
                         consensus_end=NA,
                         CO=NA),
               stringsAsFactors=FALSE)
  write.table(log, log_name , append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = TRUE)
  print(paste0(file_name, " has no hetSNPs after frame selection"))
  q()
}

imputedFrame <- hapiImupte(gmt = gmtFrame, nSPT = 2, allowNA = n_gam*2)
draftHap <- tryCatch(hapiPhase(gmt = imputedFrame), error= function(e) {return(2)}  )
if (length(draftHap) == 1){
  if (draftHap==2) {
    log <- rbind(log, list(n_gam=n_gam, n_snp=n_snp, n_het_snp=nrow(gmtDa), cov=cov, seqerr=seqerr, avgr=avgr, rs=rs,
                           gmt_import=1,
                           filter_error=1,
                           frame_select=1,
                           impute_and_phase=0,
                           final_draft=NA,
                           consenus=NA,
                           consensus_end=NA,
                           CO=NA),
                 stringsAsFactors=FALSE)
    write.table(log, log_name , append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = TRUE)
    print(paste0(file_name, " could not phase haplotypes"))
    q()
  }
}

draftHap[draftHap$cvlink>=1,]
cvCluster <- hapiCVCluster(draftHap = draftHap, cvlink = 2)

# Filtering
filter <- c()
for (i in 1:nrow(cvCluster)) {
  filter <- c(filter, which (rownames(draftHap) >= cvCluster$left[i] &
                               rownames(draftHap) <= cvCluster$right[i]))
}

if (length(filter) > 0) {
  imputedFrame <- imputedFrame[-filter, ]
  draftHap <- hapiPhase(imputedFrame)
}

### Assess gam imputation
# Read in true data
true_gf <- read.delim(paste0(gen_model_dir, '/', file_dir, '/', fname_prefix, '_gametedf_full_truth_afseqednm.csv'),
                      sep = ",", na.strings = c("NA"))
true_dh <- read.delim(paste0(gen_model_dir, '/', file_dir, '/', fname_prefix, '_donorHaps_truth_afseqednm.csv'),
                      sep=",", na.strings = c("NA"))

### Assess gam imputation

# Find which rows were removed during imputation
missing_rows <- setdiff(true_gf$positions, rownames(imputedFrame))

# Make a copy of imputed Frame and add in rows of NAs to the end of the dataframe with rownames corresponding to missing rows
hapi_gf <- imputedFrame
for (i in missing_rows){
  hapi_gf <- rbind(hapi_gf, new_row=NA)
  rownames(hapi_gf)[rownames(hapi_gf) == "new_row"] <- i
}

# Reorder so that all positions are in the correct row
hapi_gf <- hapi_gf[ order(as.numeric(row.names(hapi_gf))), ]

# Assess impuation
SAGI_hapi <- sim_assess_gam_imputation(true_gf, hapi_gf, nrow(true_dh), (ncol(true_gf) -1 ))
#Write to file
SAGI_hapi <-  append(SAGI_hapi, list(n_gametes = n_gam, n_snp = n_snp, cov = cov, seqerr=seqerr, avgr=avgr, rs=rs ), 0)
write.table(SAGI_hapi, SAGI_name , append = TRUE, sep = ",", dec = ".", row.names = TRUE, col.names = TRUE)


# Assemble final draft
finalDraft <- tryCatch(hapiBlockMPR(draftHap = draftHap, gmtFrame = gmtFrame, cvlink = 1), error= function(e) {return(2)}  )
if (length(finalDraft) == 1){
  if (finalDraft==2) {
    log <- rbind(log, list(n_gam=n_gam, n_snp=n_snp, n_het_snp=nrow(gmtDa), cov=cov, seqerr=seqerr, avgr=avgr, rs=rs,
                           gmt_import=1,
                           filter_error=1,
                           frame_select=1,
                           impute_and_phase=1,
                           final_draft=0,
                           consenus=NA,
                           consensus_end=NA,
                           CO=NA),
                 stringsAsFactors=FALSE)
    write.table(log, log_name , append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = TRUE)
    print(paste0(file_name, " could not form final draft haplotype"))
    q()
  }
}

consensusHap <- tryCatch(hapiAssemble(draftHap = finalDraft, gmt = gmtDA), error= function(e) {return(2)}  )
if (length(consensusHap) == 1){
  if (consensusHap==2) {
    log <- rbind(log, list(n_gam=n_gam, n_snp=n_snp, n_het_snp=nrow(gmtDa), cov=cov, seqerr=seqerr, avgr=avgr, rs=rs,
                           gmt_import=1,
                           filter_error=1,
                           frame_select=1,
                           impute_and_phase=1,
                           final_draft=1,
                           consenus=0,
                           consensus_end=NA,
                           CO=NA),
                 stringsAsFactors=FALSE)
    write.table(log, log_name , append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = TRUE)
    print(paste0(file_name, " could not form consensus haplotype"))
    q()
  }
}

consensusHap <- tryCatch(hapiAssembleEnd(gmt = gmtDA, draftHap = finalDraft,
                                         consensusHap = consensusHap, k = 300), error= function(e) {return(2)}  )
if (length(consensusHap) == 1){
  if (consensusHap==2) {
    log <- rbind(log, list(n_gam=n_gam, n_snp=n_snp, n_het_snp=nrow(gmtDa), cov=cov, seqerr=seqerr, avgr=avgr, rs=rs,
                           gmt_import=1,
                           filter_error=1,
                           frame_select=1,
                           impute_and_phase=1,
                           final_draft=1,
                           consenus=1,
                           consensus_end=0,
                           CO=NA),
                 stringsAsFactors=FALSE)
    write.table(log, log_name , append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = TRUE)
    print(paste0(file_name, " could not assemble end of final consensus haplotype"))
    q()
  }
}

# Remove the row of NAs from consensushap
consensusHap2 <- head(consensusHap, length(consensusHap[,1]) - 1)

# manipulate data to fit rhapsodi assess function input format
int_df <- cbind(pos=rownames(consensusHap2), consensusHap2[,c(1,2)])
fin_df <- cbind(index=rownames(int_df), int_df)
colnames(fin_df) <- c('index', 'pos', 'h1', 'h2')

# Getting CO events
snp <- which(rownames(hetDa) %in% rownames(consensusHap2))
hapOutput <- data.frame(gmt[snp,], consensusHap2)

# Separate data into haplotypes and gametes data
hap <- hapOutput[,(ncol(hapOutput)-4):(ncol(hapOutput)-3)]
gmt <- hapOutput[,5:(ncol(hapOutput)-5)]


# Rhapsodi assess

# Data import
# Read in true files

### Assess phasing

# Find which rows were removed during haplotype construction
df_missing_rows <- setdiff(rownames(gmtDa), rownames(fin_df))

# Add in new rows
hapi_dh <- fin_df %>%
  mutate(h1 = na_if(h1, 7)) %>%
  mutate(h2 = na_if(h2, 7))
hapi_dh$index <- as.character(hapi_dh$index)
hapi_dh$pos <- as.character(hapi_dh$pos)
for (i in df_missing_rows){
  hapi_dh <- rbind(hapi_dh, new_row=list(i, i,NA,NA))
  rownames(hapi_dh)[rownames(hapi_dh) == "new_row"] <- i
}

# Reorder so that all positions are in the correct row
hapi_dh <- hapi_dh[ order(as.numeric(row.names(hapi_dh))), ]

hapi_dh$index <- factor(hapi_dh$index)
hapi_dh$pos <- factor(hapi_dh$pos)

SAP_hapi <- sim_assess_phasing(true_dh, hapi_dh, nrow(true_dh))

# Modify rhapsodi analysis results for printing to file
SAP_hapi <-  append(SAP_hapi, list(n_gametes = n_gam, n_snp = n_snp, cov = cov, seqerr=seqerr, avgr=avgr, rs=rs ), 0)
# Write data to file
write.table(SAP_hapi, SAP_name, append = TRUE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)

# Assess recombination

#Identify COs
cvOutput <- tryCatch(hapiIdentifyCV(hap = hap, gmt = gmt), error= function(e) {return(2)}  )
if (length(cvOutput) == 1){
  if (cvOutput==2) {
    log <- rbind(log, list(n_gam=n_gam, n_snp=n_snp, n_het_snp=nrow(gmtDa), cov=cov, seqerr=seqerr, avgr=avgr, rs=rs,
                           gmt_import=1,
                           filter_error=1,
                           frame_select=1,
                           impute_and_phase=1,
                           final_draft=1,
                           consenus=1,
                           consensus_end=1,
                           CO=0),
                 stringsAsFactors=FALSE)
    write.table(log, log_name , append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = TRUE)
    print(paste0(file_name, " could not calculate COs"))
    ### NO NEXT STATEMENT HERE!!!!
  }
 } else {
  # Write to log
  log <- rbind(log, list(n_gam=n_gam, n_snp=n_snp, n_het_snp=nrow(gmtDa), cov=cov, seqerr=seqerr, avgr=avgr, rs=rs,
                         gmt_import=1,
                         filter_error=1,
                         frame_select=1,
                         impute_and_phase=1,
                         final_draft=1,
                         consenus=1,
                         consensus_end=1,
                         CO=1),
               stringsAsFactors=FALSE)
  write.table(log, log_name , append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = TRUE)

  true_ci <- read.delim(paste0(gen_model_dir, '/', file_dir, '/', fname_prefix, '_crossoverIndices_truth_ptfseqednm.csv'),
                        sep = ",", na.strings = c("NA"))
  hapi_ci <-  cvOutput[, c(1:3)]
  colnames(hapi_ci) <- c('Ident', 'Genomic_start', 'Genomic_end')
  hapi_ci <- mutate(hapi_ci, Ident = paste0("sampleT_chrT_", Ident))

  # Get full list of what the rownames to Ident should be
  true_Ident_list = c()
  for (i in 1:(ncol(true_gf) -1 )){
    true_Ident_list <- c(true_Ident_list, paste0("sampleT_chrT_gam", i, "_"))
  }

  # Get list of gametes that do not have crossovers
  gametes_no_COs <- setdiff(true_Ident_list, hapi_ci$Ident)

  # Add rows corresponding to crossover-less gametes
  for (i in gametes_no_COs){
    hapi_ci <- rbind(hapi_ci, list(i, NA,NA))
  }

  # Assess reocmbination prediction accuracy
  SAR_hapi <- sim_assess_recomb(true_ci, hapi_ci, cons = FALSE)

  SAR_hapi <- append(SAR_hapi, list(n_gametes = n_gam, n_snp = n_snp, cov = cov, seqerr=seqerr, avgr=avgr, rs=rs ), 0)
  write.table(SAR_hapi, SAR_name , append = TRUE, sep = ",", dec = ".", row.names = TRUE, col.names = TRUE)
}