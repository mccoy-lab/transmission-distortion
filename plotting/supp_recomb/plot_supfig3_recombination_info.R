library(tidyverse)
#library(ggplot2)
library(data.table)
#library(gridExtra)
library(itertools)
library(rhapsodi)

dir_base <- "/home/kweave23/gamete_data/gen_model_results_noDNM/"
#dir_base <- "/home-3/kweave23@jhu.edu/work-rmccoy22/kweave23/sc_transmission_distortion/generative_model_noDNM/rerun_with_rhapsodi_out_from_gen_model_results_noDNM/"
#dir_base2 <- "/home-3/kweave23@jhu.edu/work-rmccoy22/kweave23/sc_transmission_distortion/generative_model_noDNM/"
true_snp_file <- "true_nsnp_gen_model_noDNM.Rdata"

#load(paste0(dir_base2, true_snp_file))
load(paste0(dir_base, true_snp_file))

nsnps <- c(5000, 30000, 100000)
nsnps_iter <- as.list(itertools::enumerate(nsnps))
ngams <- c(3, 15, 50, 150, 500, 1000, 2500, 5000)
ngams_iter <- as.list(itertools::enumerate(ngams))
covs <- c(0.001, 0.01, 0.1, 0.223, 0.357, 0.511, 0.693, 1.204, 2.303)
covs_iter <- as.list(itertools::enumerate(covs))
rsds <- c(42, 357, 1848)
gen_seqerror <- 0.005
#mod_seqerror <- 0.005
gen_avgrecomb <- 1
#mod_avgrecomb <- 1

find_fp <- function(pred_intersect_dt, pred_dt_nona, no_truths=FALSE, no_preds=FALSE){
  if (!no_truths & !no_preds){
    fp <- pred_intersect_dt[is.na(pred_intersect_dt$True_Start),] #number predicted as recombination spots, but don't intersect with truth at all
  } else if(no_truths & !no_preds){
    fp <- pred_dt_nona #no truths, but some predictions; i.e. no non-na truths, but there are non-na predictions, all are false positives
  } else if(!no_truths & no_preds){
    fp <- NA
  } else { fp <- NA}
  return (fp)
}

find_fn <- function(truth_intersect_dt, truth_dt_nona, no_truths=FALSE, no_preds=FALSE){
  if (!no_truths & !no_preds){
    fn <- truth_intersect_dt[is.na(truth_intersect_dt$Predicted_Start),] #number of true recombinantion spots that don't interesect any predictions
  } else if (no_truths & !no_preds){
    fn <- NA
  } else if (!no_truths & no_preds){
    fn <- truth_dt_nona ##no predictions, but some truths, i.e. no non-na predictions but there are non-na truths, all are false negatives
  } else {fn <- NA}
  return (fn)
}

assess_recomb <- function(true_recomb, pred_recomb){
  ##Recombination Discovery assessment
  true_recomb <- data.table::data.table(true_recomb)
  true_recomb_nona <- true_recomb[!is.na(true_recomb$start),] %>% data.table::setkey()
  true_recomb_na <- true_recomb[is.na(true_recomb$start),]
  
  pred_recomb_dt <- data.table(gam=sapply(strsplit(pred_recomb$Ident, "_"), `[`, 3), start=pred_recomb$Genomic_start, end=pred_recomb$Genomic_end)
  pred_recomb_nona <- pred_recomb_dt[!is.na(pred_recomb_dt$start),] %>% data.table::setkey()
  pred_recomb_na <- pred_recomb_dt[is.na(pred_recomb_dt$start),]
  
  if (nrow(true_recomb_nona) > 0){
    no_truths = FALSE
  } else {no_truths = TRUE}
  if (nrow(pred_recomb_nona) > 0){
    no_preds = FALSE
  } else {no_preds = TRUE}
  
  if (!no_truths & !no_preds){
    truth_intersect <- data.table::foverlaps(true_recomb_nona, pred_recomb_nona) %>% `colnames<-`(c("gam", "Predicted_Start", "Predicted_End", "True_Start", "True_End"))
    pred_intersect <- data.table::foverlaps(pred_recomb_nona, true_recomb_nona) %>% `colnames<-`(c("gam", "True_Start", "True_End", "Predicted_Start", "Predicted_End"))
  } else{
    truth_intersect <- NULL
    pred_intersect <- NULL
  }
  
  fn <- find_fn(truth_intersect, true_recomb_nona, no_truths = no_truths, no_preds = no_preds)
  fp <- find_fp(pred_intersect, pred_recomb_nona, no_truths = no_truths, no_preds = no_preds)
  
  return(list(fn_instances = fn, fp_instances = fp))
}

resolution_vecs <- c()
resolution_vecs_norm <- c()
fns_locs <- c(starts = c(), ends = c())
fps_locs <- c(starts = c(), ends = c())
name_gams <- c()
name_snps <- c()
name_covs <- c()
for (nsnp in nsnps_iter){
  for (ngam in ngams_iter){
    for (cov in covs_iter){
      dir_data <- paste0("g", ngam$value, "_s", as.character(as.integer(nsnp$value)), "_c", cov$value, "_se", gen_seqerror, "_r", gen_avgrecomb, "/")
      resolution_vec <- c()
      resolution_vec_norm <- c()
      locs_of_fps <- c(starts = c(), ends = c())
      locs_of_fns <- c(starts = c(), ends = c())
      for (rsd in rsds){
        sim_base <- paste0("/home/kweave23/gamete_data/gen_model_results_noDNM/g",ngam$value,"_s",as.character(as.integer(nsnp$value)),"_c",cov$value,"_se",gen_seqerror,"_r",gen_avgrecomb,"/runGen_gam_",ngam$value,"_snp_", as.character(as.integer(nsnp$value)), "_cov_", cov$value, "_seqerr_", gen_seqerror, "_avgr_", gen_avgrecomb, "_rs_", rsd)
        #sim_base <- paste0("/home-3/kweave23@jhu.edu/work-rmccoy22/kweave23/sc_transmission_distortion/generative_model_noDNM/rerun_with_rhapsodi_out_from_gen_model_results_noDNM/g",ngam$value,"_s",nsnp$value,"_c",cov$value,"_se",gen_seqerror,"_r",gen_avgrecomb,"/runGen_gam_",ngam$value,"_snp_", nsnp$value, "_cov_", cov$value, "_seqerr_", gen_seqerror, "_avgr_", gen_avgrecomb, "_rs_", rsd)
        input_file <- paste0(sim_base,"_gametedf_na_truth_afseqednm.csv")
        dt <- read.delim(input_file, sep=",", na.strings = c("NA"))
        standard_out <- rhapsodi::standard_input(NULL, use_dt = TRUE, input_dt = dt)
        
        
        data_file_pred <- paste0("rhapsodi_out_rs_", rsd ,".Rdata")
        data_file_true <- paste0("runGen_gam_", ngam$value, "_snp_", as.character(as.integer(nsnp$value)), "_cov_", cov$value, "_seqerr_", gen_seqerror, "_avgr_", gen_avgrecomb, "_rs_", rsd, "_crossoverIndices_truth_ptfseqednm.csv")
        fileOIpred <- paste0(dir_base, dir_data, data_file_pred) #predictions
        fileOItrue <- paste0(dir_base, dir_data, data_file_true) #truths
        if (file_test("-f", fileOIpred) & file_test("-f", fileOItrue)){
          true_snp_number <- true_nsnp_dict[[format(nsnp$value, scientific = FALSE, trim = TRUE)]][[as.character(ngam$value)]][[as.character(cov$value)]][[as.character(gen_seqerror)]][[as.character(gen_avgrecomb)]][[as.character(rsd)]]
          #load fileOIs
          load(fileOIpred)
          pred_dt <- rhapsodi_out$recomb_breaks
          new_starts <- unlist(lapply(1:length(pred_dt$Genomic_start), function(x) match(pred_dt$Genomic_start[x], standard_out$positions)))
          new_ends <- unlist(lapply(1:length(pred_dt$Genomic_end), function(x) match(pred_dt$Genomic_end[x], standard_out$positions)))
          #compute resolution/length of predictions
          resolution_vec <- c(resolution_vec, c(new_ends[!is.na(new_ends)] - new_starts[!is.na(new_starts)] + 1))
          resolution_vec_norm <- c(resolution_vec_norm, c((new_ends[!is.na(new_ends)] - new_starts[!is.na(new_starts)] + 1)/true_snp_number))
          #find locs for which ones are missed
          true_dt <- read.delim(fileOItrue, sep = ",", na.strings = c("NA"))
          #find locs for which ones are added (fps)
          ar_out <- assess_recomb(true_dt, pred_dt)
          fp_obs <- ar_out$fp_instances
          if (!is.null(dim(fp_obs))){
            new_fp_starts <- unlist(lapply(1:length(fp_obs$Predicted_Start), function(x) match(fp_obs$Predicted_Start[x], standard_out$positions)))
            new_fp_ends <- unlist(lapply(1:length(fp_obs$Predicted_End), function(x) match(fp_obs$Predicted_End[x], standard_out$positions)))
            locs_of_fps[["starts"]] <- c(locs_of_fps[["starts"]], new_fp_starts / true_snp_number)
            locs_of_fps[["ends"]] <- c(locs_of_fps[["ends"]], new_fp_ends / true_snp_number)
          }
          fn_obs <- ar_out$fn_instances
          if (!is.null(dim(fn_obs))){
            new_fn_starts <- unlist(lapply(1:length(fn_obs$True_Start), function(x) match(fn_obs$True_Start[x], standard_out$positions)))
            new_fn_ends <- unlist(lapply(1:length(fn_obs$True_End), function(x) match(fn_obs$True_End[x], standard_out$positions)))
            locs_of_fns[["starts"]] <- c(locs_of_fns[["starts"]], new_fn_starts / true_snp_number)
            locs_of_fns[["ends"]] <- c(locs_of_fns[["ends"]], new_fn_ends / true_snp_number)
          }
        }
      }
      name_gams <- c(name_gams, ngam$value)
      name_snps <- c(name_snps, nsnp$value)
      name_covs <- c(name_covs, cov$value)
      resolution_vecs <- c(resolution_vecs, list(resolution_vec))
      resolution_vecs_norm <- c(resolution_vecs_norm, list(resolution_vec_norm))
      fns_locs[["starts"]] <- c(fns_locs[["starts"]], list(locs_of_fns[["starts"]]))
      fns_locs[["ends"]] <- c(fns_locs[["ends"]], list(locs_of_fns[["ends"]]))
      fps_locs[["starts"]] <- c(fps_locs[["starts"]], list(locs_of_fps[["starts"]]))
      fps_locs[["ends"]] <- c(fps_locs[["ends"]], list(locs_of_fps[["ends"]]))
    }
  }
}

roi <- lapply(1:length(resolution_vecs), function(x) identical(resolution_vecs[[x]], numeric(0)))
resolution_vecs[unlist(roi)] <- NA
resolution_vecs_norm[unlist(roi)] <- NA

save(name_gams, file='supfig3_name_gams.Rdata')
save(name_snps, file='supfig3_name_snps.Rdata')
save(name_covs, file='supfig3_name_covs.Rdata')
save(resolution_vecs, file='supfig3_res_vecs.Rdata')
save(resolution_vecs_norm, file='supfig3_res_vecs_norm.Rdata')
save(fps_locs, file='supfig3_fp_locs.Rdata')
save(fns_locs, file='supfig3_fn_locs.Rdata')

# df_res <- data.frame(name_gams, name_snps, name_covs, list(resolution_vecs))
# df_res_norm <- data.frame(name_gams, name_snps, name_covs, list(resolution_vecs_norm))
# df_fp <- data.frame(name_gams, name_snps, name_covs, fps_locs[["starts"]], fps_locs[["ends"]])
# df_fn <- data.frame(name_gams, name_snps, name_covs, fns_locs[["starts"]], fns_locs[["ends"]])
# 
# save(df_res, file = "supfig3_res.Rdata")
# save(df_res_norm, file = "supfig3_res_norm.Rdata")
# save(df_fp, file = "supfig3_fp_locs.Rdata")
# save(df_fn, file = "supfig3_fn_locs.Rdata")
