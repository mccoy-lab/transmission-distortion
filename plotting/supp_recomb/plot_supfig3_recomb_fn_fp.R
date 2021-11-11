library(tidyverse)
library(data.table)
library(rhapsodi)
library(ggplot2)

load_data <- TRUE

## are these fp and fn locations co-occurring in the same gametes?
covs <- c(0.001, 0.01, 0.1, 0.223, 0.357, 0.511, 0.693, 1.204, 2.303)
rsds <- c(42, 357, 1848)
nsnp <- 30000
ngam <- 1000
gen_seqerror <- 0.005
#mod_seqerror <- 0.005
gen_avgrecomb <- 1
#mod_avgrecomb <- 1

if(!load_data){
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
      standard_out <- rhapsodi::read_data(NULL, use_dt = TRUE, input_dt = dt)
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
  save(storage_list, file="fn_fp_storage_list.Rdata")
} else {
  load("fn_fp_storage_list.Rdata")
} 
  
plot_fn_fp <- function(covs, rsds, storage_list){
  tall_df <- data.frame()
  for (cov in covs){
    num_uniq_gametes <- 0
    for (rsd in rsds){
      fp_info <- storage_list[[paste0("c",cov)]][[paste0("df",rsd)]]$fp
      fp_info$gam <- unlist(lapply(1:length(storage_list[[paste0("c",cov)]][[paste0("df",rsd)]]$fp$gam), function(x) paste0(storage_list[[paste0("c",cov)]][[paste0("df",rsd)]]$fp$gam[x], "_", rsd)))
      fn_info <- storage_list[[paste0("c",cov)]][[paste0("df",rsd)]]$fn
      fn_info$gam <- unlist(lapply(1:length(storage_list[[paste0("c",cov)]][[paste0("df",rsd)]]$fn$gam), function(x) paste0(storage_list[[paste0("c",cov)]][[paste0("df",rsd)]]$fn$gam[x], "_", rsd)))
      melted_df <- melt(merge(as.data.table(fn_info), as.data.table(fp_info), all = TRUE, by = "gam"))
      melted_df$coverage <- paste0("Coverage: ", cov)
      tall_df <- rbind(tall_df, melted_df)
      num_uniq_gametes <- num_uniq_gametes + length(unique(melted_df$gam))
    }
    message(paste0("Coverage: ", cov, " num_gametes: ", num_uniq_gametes))
  }
  
  tall_df$variable <- toupper(tall_df$variable)
  #save(tall_df, file="fn_fp_tall_df.Rdata")
  
  p <- ggplot(tall_df, aes(value, as.factor(gam), shape=factor(variable), size=factor(variable))) + scale_size_manual(values=c("FN"=1.5, "FP"=0.75)) + scale_shape_manual(values=c("FN"=19, "FP"=20)) + theme_bw() + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line( size=0.05, color="gray"), panel.grid = element_blank(), axis.text.y=element_blank()) + facet_wrap(~coverage, ncol=2, scales="free_y")
  new_p <- p + geom_point(aes(colour = factor(variable))) + scale_colour_manual(values=c("FN"="tomato", "FP"="black")) + labs(x='Relative chromosomal location', y="Gamete", color="prediction class", shape="prediction class", size="prediction class")
  return(new_p)
}



covs2 <- c(0.01, 0.1, 0.693, 1.204)
new_p <- plot_fn_fp(covs2, rsds, storage_list)
new_p
ggsave( "flipped_facet_wrap_2col_fn_fp.pdf" , plot=new_p, width = 10, height = 7, units = "in")

covs3 <- c(0.01, 0.1)
p3 <- plot_fn_fp(covs3, rsds, storage_list)
p3
ggsave('bioinfo_fn_fp.pdf', plot=p3)
