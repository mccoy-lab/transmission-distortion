library(tidyverse)
library(data.table)
library(rhapsodi)
library(ggplot2)

## are these fp and fn locations co-occurring in the same gametes?
covs <- c(0.001, 0.01, 0.1, 0.223, 0.357, 0.511, 0.693, 1.204, 2.303)
rsds <- c(42, 357, 1848)
nsnp <- 30000
ngam <- 1000
gen_seqerror <- 0.005
#mod_seqerror <- 0.005
gen_avgrecomb <- 1
#mod_avgrecomb <- 1

#dir_base <- "/home/kweave23/gamete_data/gen_model_results_noDNM/"
dir_base <- "/Users/kateweaver/mccoyLab_withOthers/transmission-distortion/test_data_rhapsodi_gen/"
true_snp_file <- "true_nsnp_gen_model_noDNM.Rdata"
load(paste0(dir_base, true_snp_file))

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

storage_list <- list()
#dir_data <- paste0("g", ngam, "_s", as.character(as.integer(nsnp)), "_c", cov, "_se", gen_seqerror, "_r", gen_avgrecomb, "/")
for (cov in covs){
  storage_list[[paste0("c", cov)]] <- list()
  for (rsd in rsds){
    #sim_base <- paste0("/home/kweave23/gamete_data/gen_model_results_noDNM/g",ngam,"_s",as.character(as.integer(nsnp)),"_c",cov,"_se",gen_seqerror,"_r",gen_avgrecomb,"/runGen_gam_",ngam,"_snp_", as.character(as.integer(nsnp)), "_cov_", cov, "_seqerr_", gen_seqerror, "_avgr_", gen_avgrecomb, "_rs_", rsd)
    sim_base <- paste0("bell_sim_data/")
    file_base <- paste0("runGen_gam_",ngam,"_snp_", as.character(as.integer(nsnp)), "_cov_", cov, "_seqerr_", gen_seqerror, "_avgr_", gen_avgrecomb, "_rs_", rsd)
    input_file <- paste0(dir_base, sim_base, file_base, "_gametedf_na_truth_afseqednm.csv")
    dt <- read.delim(input_file, sep=",", na.strings = c("NA"))
    standard_out <- rhapsodi::standard_input(NULL, use_dt = TRUE, input_dt = dt)
    data_file_pred <- paste0("rhapsodi_out_rs_", rsd ,".Rdata")
    data_file_true <- paste0("_crossoverIndices_truth_ptfseqednm.csv")
    #fileOIpred <- paste0(dir_base, dir_data, data_file_pred) #predictions
    fileOIpred <- paste0(dir_base, sim_base, "c", cov, "/", data_file_pred)
    #fileOItrue <- paste0(dir_base, dir_data, data_file_true) #truths
    fileOItrue <- paste0(dir_base, sim_base, file_base, data_file_true)
    if (file_test("-f", fileOIpred) & file_test("-f", fileOItrue)){
      true_snp_number <- true_nsnp_dict[[format(nsnp, scientific = FALSE, trim = TRUE)]][[as.character(ngam)]][[as.character(cov)]][[as.character(gen_seqerror)]][[as.character(gen_avgrecomb)]][[as.character(rsd)]]
      #load fileOIs
      load(fileOIpred)
      pred_dt <- rhapsodi_out$recomb_breaks
      true_dt <- read.delim(fileOItrue, sep = ",", na.strings = c("NA"))
      #find locs for which ones are added (fps)
      ar_out <- assess_recomb(true_dt, pred_dt)
      fp_obs <- ar_out$fp_instances
      if (!is.null(dim(fp_obs))){
        new_fp_starts <- lapply(1:length(fp_obs$Predicted_Start), function(x) match(fp_obs$Predicted_Start[x], standard_out$positions) / true_snp_number )
        fp_df <- do.call(rbind.data.frame, new_fp_starts) %>% `rownames<-`(make.unique(fp_obs$gam)) %>% `colnames<-`("fp")
        fp_df$gam <- sapply(strsplit(rownames(fp_df), "[.]"), `[`, 1)
      } else { fp_df <- data.frame()}
      fn_obs <- ar_out$fn_instances
      if (!is.null(dim(fn_obs))){
        new_fn_starts <- lapply(1:length(fn_obs$True_Start), function(x) which(abs(fn_obs$True_Start[x] - standard_out$positions) == min(abs(fn_obs$True_Start[x] - standard_out$positions))) / true_snp_number )
        names(new_fn_starts) <- sapply(strsplit(make.unique(fn_obs$gam), "[.]"), `[`, 1)
        fn_df <- data.frame(gam=names(unlist(new_fn_starts)), fn=unlist(new_fn_starts))
      } else { fn_df <- data.frame() }
      storage_name <- paste0("df", rsd)
      storage_list[[paste0("c", cov)]][[storage_name]] <- list("fp"=fp_df, "fn"=fn_df)
    } else {message("no files")}
  }
}

bell_named_list_fn <- list()
for (i in 1:length(covs)){
  bell_long_df <- data.frame()
  cov <- covs[i]
  for (rsd in rsds){
    bell_long_df <- rbind(bell_long_df, data.frame(fn=storage_list[[paste0("c",cov)]][[paste0("df",rsd)]]$fn, cov=cov))
  }
  bell_named_list_fn[[paste0("c", cov)]] <- bell_long_df
}

bell_long_fn <- rbindlist(bell_named_list_fn)

p <- ggplot(bell_long_fn, aes(as.factor(cov), fn.fn)) + theme_bw() + theme(panel.grid = element_blank()) 
p + geom_point() + labs(x = 'Coverage (x)', y='relative location of false negatives')
ggsave('bell_sim_breakpoint_fn_full.png')


bell_named_list_fp <- list()
for (i in 1:length(covs)){
  bell_long_df <- data.frame()
  cov <- covs[i]
  for (rsd in rsds){
    bell_long_df <- rbind(bell_long_df, data.frame(fp=storage_list[[paste0("c",cov)]][[paste0("df",rsd)]]$fp, cov=cov))
  }
  bell_named_list_fp[[paste0("c", cov)]] <- bell_long_df
}

bell_long_fp <- rbindlist(bell_named_list_fp)

p <- ggplot(bell_long_fp, aes(as.factor(cov), fp.fp)) + theme_bw() + theme(panel.grid = element_blank()) 
p + geom_point() + labs(x = 'Coverage (x)', y='relative location of false positives')
ggsave('bell_sim_breakpoint_fp_full.png')

plot_fn_fp <- function(covs, rsds, storage_list){
  tall_df <- data.frame()
  for (cov in covs){
    for (rsd in rsds){
      melted_df <- melt(merge(storage_list[[paste0("c",cov)]][[paste0("df",rsd)]]$fn, storage_list[[paste0("c",cov)]][[paste0("df",rsd)]]$fp, all = TRUE, by = "gam"))
      melted_df$coverage <- paste0("Coverage: ", cov)
      tall_df <- rbind(tall_df, melted_df)  
    }
  }
  p <- ggplot(tall_df, aes(as.factor(gam), value, shape=factor(variable), size=factor(variable))) + scale_size_manual(values=c("fn"=1.25, "fp"=1.05)) + scale_shape_manual(values=c("fn"=19, "fp"=3)) + theme_bw() + theme(panel.grid = element_blank(), axis.text.x=element_blank()) + facet_wrap(~coverage, ncol=length(covs), scales="free_x")
  new_p <- p + geom_point(aes(colour = factor(variable)), size=1.25, stroke=1.5) + scale_colour_manual(values=c("fn"="tomato", "fp"="black")) + labs(y='relative location', x="gamete")
  
  ggsave("fn_fp_by_gamete.png")
  return(new_p)
}

covs2 <- c(0.01, 0.1, 0.693, 1.204)
new_p <- plot_fn_fp(covs2, rsds, storage_list)
new_p


