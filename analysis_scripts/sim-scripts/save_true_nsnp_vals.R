true_nsnp_dict = list()
for (nsnp_start in c("5000","30000","100000")){
  true_nsnp_dict[[nsnp_start]] = list()
  for (ngam in c("3", "15", "50", "150", "500", "1000", "2500", "5000")){
    true_nsnp_dict[[nsnp_start]][[ngam]] = list()
    for (cov in c("0.001", "0.01", "0.1", "0.223", "0.357", "0.511", "0.693", "1.204", "2.303")){
      true_nsnp_dict[[nsnp_start]][[ngam]][[cov]] = list()
      for (seqe in c("0.001", "0.005", "0.05")){
        true_nsnp_dict[[nsnp_start]][[ngam]][[cov]][[seqe]] = list()
        for (avgr in c("0.6", "1", "3")){
          true_nsnp_dict[[nsnp_start]][[ngam]][[cov]][[seqe]][[avgr]] = list()
          for (rsd in c("42", "357", "1848")){
            input_file <- paste0("/home/kweave23/gamete_data/gen_model_results_noDNM/g", ngam, "_s", nsnp_start, "_c", cov, "_se", seqe, "_r", avgr, "/runGen_gam_", ngam, "_snp_", nsnp_start, "_cov_", cov, "_seqerr_", seqe, "_avgr_", avgr, "_rs_", rsd, "_gametedf_na_truth_afseqednm.csv")
            dt <- read.delim(input_file, sep=",", na.strings = c("NA"))
            nsnp_true <- nrow(dt)
            true_nsnp_dict[[nsnp_start]][[ngam]][[cov]][[seqe]][[avgr]][[rsd]] = nsnp_true
            rm(dt)
          }
        }
      }
    }
  }
}

save(true_nsnp_dict, file="true_nsnp_gen_model_noDNM.Rdata")
