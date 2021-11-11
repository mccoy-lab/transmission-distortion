library(tidyverse)
library(ggplot2)
library(pbmcapply)
library(data.table)
library(plyr)

load_data <- TRUE

if (!load_data){
  load("/Users/kateweaver/mccoyLab_withOthers/transmission-distortion/plotting/supfig3_rdata/supfig3_name_covs.Rdata")
  load("/Users/kateweaver/mccoyLab_withOthers/transmission-distortion/plotting/supfig3_rdata/supfig3_name_gams.Rdata")
  load("/Users/kateweaver/mccoyLab_withOthers/transmission-distortion/plotting/supfig3_rdata/supfig3_name_snps.Rdata")
  load("/Users/kateweaver/mccoyLab_withOthers/transmission-distortion/plotting/supfig3_rdata/supfig3_res_vecs.Rdata")

  meta_df <- data.frame(name_gams, name_snps, name_covs)
  roi <- meta_df$name_gams == 1000 & meta_df$name_snps == 30000
  bell_df <- meta_df[roi,]
  bell_res_snp <- resolution_vecs[roi]
  bell_named_list_res_snp <- list()
  for (i in 1:length(bell_res_snp)){
    bell_long_snp <- data.frame()
    ngam <- bell_df$name_gams[i]
    nsnp <- bell_df$name_snps[i]
    cov <- bell_df$name_covs[i]
    res_snp <- bell_res_snp[i][[1]]
    bell_named_list_res_snp[[i]] <- rbindlist(pbmclapply(1:length(res_snp),
                                                        function(x) if (!is.na(res_snp[x])) rbind(bell_long_snp, list(name_ngam=ngam, name_nsnp=nsnp, name_cov=cov, snp_res=res_snp[x])),
                                                        mc.cores = getOption("mc.cores", 4)))
  }

  bell_long_snp_final <- rbindlist(bell_named_list_res_snp)
  bell_long_snp_final$bp <- (bell_long_snp_final$snp_res - 1) * 1000

  save(bell_long_snp_final, file="bell_long_snp_final.Rdata")
} else {
  load("bell_long_snp_final.Rdata")
}
medians <- ddply(bell_long_snp_final, .(as.factor(name_cov)), summarise, snp_med = median(snp_res), bp_med = median(bp)) %>% `colnames<-`(c("cov", "snp_med", "bp_med"))
      
p <- ggplot(bell_long_snp_final, aes(as.factor(name_cov), bp)) + theme_bw() + theme(panel.grid = element_blank())
p2 <- p + geom_boxplot() + geom_hline(yintercept = (2-1)*1000, linetype = 'dashed') + labs(x = 'Coverage (x)', y="Crossover breakpoint resolution (bp)") + annotate("text", x = "0.001", y = (2-1)*1000, label = "2 SNP resolution", vjust = -0.5, hjust = "center", size=2.5)
p3 <- p2 + scale_y_log10() + annotation_logticks(scaled=TRUE, sides='l')
show(p3)
ggsave('bell_sim_breakpoint_res_bp_10_5.pdf', plot = p3, width=10, height=5, units="in")
