library(tidyverse)
library(ggplot2)
library(pbmcapply)
library(data.table)
library(plyr)

load("/Users/kateweaver/mccoyLab_withOthers/transmission-distortion/plotting/supfig3_rdata/supfig3_name_covs.Rdata")
load("/Users/kateweaver/mccoyLab_withOthers/transmission-distortion/plotting/supfig3_rdata/supfig3_name_gams.Rdata")
load("/Users/kateweaver/mccoyLab_withOthers/transmission-distortion/plotting/supfig3_rdata/supfig3_name_snps.Rdata")
load("/Users/kateweaver/mccoyLab_withOthers/transmission-distortion/plotting/supfig3_rdata/supfig3_res_vecs_norm.Rdata")
load("/Users/kateweaver/mccoyLab_withOthers/transmission-distortion/plotting/supfig3_rdata/supfig3_res_vecs.Rdata")


meta_df <- data.frame(name_gams, name_snps, name_covs)
roi <- meta_df$name_gams == 1000 & meta_df$name_snps == 30000
bell_df <- meta_df[roi,]
bell_res <- resolution_vecs_norm[roi]
bell_res_snp <- resolution_vecs[roi]

bell_named_list_res <- list()
bell_named_list_res_snp <- list()
for (i in 1:length(bell_res)){
  bell_long_df <- data.frame()
  bell_long_snp <- data.frame()
  ngam <- bell_df$name_gams[i]
  nsnp <- bell_df$name_snps[i]
  cov <- bell_df$name_covs[i]
  res_vec <-bell_res[i][[1]]
  res_snp <- bell_res_snp[i][[1]]
  bell_named_list_res[[i]] <- rbindlist(pbmclapply(1:length(res_vec), 
                                                   function(x) if (!is.na(res_vec[x])) rbind(bell_long_df, list(name_ngam=ngam, name_nsnp=nsnp, name_cov=cov, name_res=res_vec[x])), 
                                                   mc.cores = getOption("mc.cores", 4)))
  bell_named_list_res_snp[[i]] <- rbindlist(pbmclapply(1:length(res_snp),
                                                       function(x) if (!is.na(res_snp[x])) rbind(bell_long_snp, list(name_ngam=ngam, name_nsnp=nsnp, name_cov=cov, name_res=res_snp[x])),
                                                       mc.cores = getOption("mc.cores", 4)))
}
bell_long_df_final <- rbindlist(bell_named_list_res)
bell_long_snp_final <- rbindlist(bell_named_list_res_snp)

bell_long_df_final$log_res <- log2(bell_long_df_final$name_res)

medians <- merge(merge(ddply(bell_long_df_final, .(as.factor(name_cov)), summarise, approx_med = median(name_res)*30000) %>% `colnames<-`(c("cov","approx_med")), 
                       ddply(bell_long_df_final, .(as.factor(name_cov)), summarise, log_med = median(log_res)) %>% `colnames<-`(c("cov", "log_med"))),
                 ddply(bell_long_snp_final, .(as.factor(name_cov)), summarise, true_med = median(name_res)) %>% `colnames<-`(c("cov", "true_med")))

p <- ggplot(bell_long_df_final, aes(as.factor(name_cov), log_res)) + theme_bw() + theme(panel.grid = element_blank()) 
p2 <- p + geom_boxplot() + geom_hline(yintercept=log2(2/30000), linetype='dashed') + labs(x = 'Coverage (x)', y='log2(resolution of breakpoint/total SNPs)') + annotate("text", x = "0.001", y = log2(2/30000), label = "2 SNP resolution", vjust = -0.5, hjust = "center", size=2.5)
p3 <- p2 + geom_text(data= medians, aes(x=cov, y=log_med+0.2, label=as.integer(true_med)))
show(p3)
ggsave('bell_sim_breakpoint_res.png', width=14, height=6, units="in")


load('/Users/kateweaver/mccoyLab_withOthers/transmission-distortion/plotting/supfig3_rdata/supfig3_fn_locs.Rdata')
bell_fn <- fns_locs$starts[roi]
bell_named_list_fn <- list()
for (i in 1:length(bell_fn)){
  bell_long_df2 <- data.frame()
  ngam <- bell_df$name_gams[i]
  nsnp <- bell_df$name_snps[i]
  cov <- bell_df$name_covs[i]
  fn_vec <- bell_fn[i][[1]]
  bell_named_list_fn[[i]] <- rbindlist(pbmclapply(1:length(fn_vec), 
                                                  function(x) if (!is.na(fn_vec[x])) rbind(bell_long_df2, list(name_ngam=ngam, name_nsnp=nsnp, name_cov=cov, name_fn=fn_vec[x])), 
                                                  mc.cores = getOption("mc.cores", 4)))
}
bell_long_df2_final <- rbindlist(bell_named_list_fn)

p <- ggplot(bell_long_df2_final, aes(as.factor(name_cov), name_fn)) + theme_bw() + theme(panel.grid = element_blank()) 
p + geom_point() + labs(x = 'Coverage (x)', y='relative location of false negatives')
ggsave('bell_sim_breakpoint_fn.png')


load('/Users/kateweaver/mccoyLab_withOthers/transmission-distortion/plotting/supfig3_rdata/supfig3_fp_locs.Rdata')
bell_fp <- fps_locs$starts[roi]
bell_named_list_fp = list()
for (i in 1:length(bell_fp)){
  bell_long_df3 <- data.frame()
  ngam <- bell_df$name_gams[i]
  nsnp <- bell_df$name_snps[i]
  cov <- bell_df$name_covs[i]
  fp_vec <- bell_fp[i][[1]]
  bell_named_list_fp[[i]] <- rbindlist(pbmclapply(1:length(fp_vec), 
                                                  function(x) if (!is.na(fp_vec[x])) rbind(bell_long_df3, list(name_ngam=ngam, name_nsnp=nsnp, name_cov=cov, name_fp=fp_vec[x])), 
                                                  mc.cores = getOption("mc.cores", 4)))
}
bell_long_df3_final <- rbindlist(bell_named_list_fp)

p <- ggplot(bell_long_df3_final, aes(as.factor(name_cov), name_fp)) + theme_bw() + theme(panel.grid = element_blank()) 
p + geom_point() + labs(x = 'Coverage (x)', y='relative location of false positives')
ggsave('bell_sim_breakpoint_fp.png')

