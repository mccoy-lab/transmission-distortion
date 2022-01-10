library(rhapsodi)
library(data.table)
library(HMM)
library(dplyr)
library(tidyr)
library(tibble)
library(stats)
library(parallel)
library(readr)
library(utils)
library(magrittr)
library(pbapply)
library(pbmcapply)
library(ggplot2)
library(patchwork)
n_snp = 1000
#TD_gams <- sim_run_generative_model_with_TD(num_gametes = 1000, num_snps = n_snp, coverage = 0.01, p_kill = 0.3)

TD_gams <- sim_run_generative_model_with_TD(num_gametes = 1000, num_snps = n_snp, coverage = 1, conversion_lambda = 5, p_convert = 0.8, TD_type = 'convert')
a <- TD_gams[[1]]
killer_snp <- TD_gams[[2]]
b <- rhapsodi_autorun(input_dt = a[[3]], use_dt = TRUE)

td_test <- function(sperm_matrix, row_index) {
  test_row <- sperm_matrix[row_index,]
  gt_vector <- unlist(test_row)[-1]
  one_count <- sum(gt_vector == 'h1' | gt_vector == 1 , na.rm = TRUE)
  two_count <- sum(gt_vector == 'h2' | gt_vector == 0, na.rm = TRUE)
  p_value <- binom.test(c(one_count, two_count))$p.value
  return(c(p_value, one_count, two_count))
}

true_gametes <- a[[2]]
true_gametes_df <- as.data.frame(as.matrix(true_gametes))
rownames(true_gametes_df) <- true_gametes_df[,1]
true_gametes_df <- true_gametes_df[, 2:length(true_gametes_df)]
colnames(true_gametes_df) <- c(sprintf("g_%d", seq(1:ncol(true_gametes_df))))

df_counts_pvals <- do.call(rbind, pbmclapply(1:nrow(true_gametes_df),
                                             function(x) td_test(true_gametes_df, x),
                                             mc.cores=getOption("mc.cores", 2))) %>%
  as_tibble() %>%
  add_column(rownames(true_gametes_df)) #bind the positions vector to df_counts_pvals
colnames(df_counts_pvals) <- c("pval", "h1_count", "h2_count", "genomic_position")

imputed_gametes <- b[[4]]
imputed_gametes_df <- as.data.frame(as.matrix(imputed_gametes))
rownames(imputed_gametes_df) <- imputed_gametes_df$pos
imputed_gametes_df <- imputed_gametes_df[, 3:length(imputed_gametes_df)]
colnames(imputed_gametes_df) <- c(sprintf("g_%d", seq(1:ncol(imputed_gametes_df))))

df_counts_pvals_rh <- do.call(rbind, pbmclapply(1:nrow(imputed_gametes_df),
                                             function(x) td_test(imputed_gametes_df, x),
                                             mc.cores=getOption("mc.cores", 2))) %>%
  as_tibble() %>%
  add_column(rownames(imputed_gametes_df)) #bind the positions vector to df_counts_pvals
colnames(df_counts_pvals_rh) <- c("pval", "h1_count", "h2_count", "genomic_position")

# Bonferroni
p_cutoff = 0.05 / n_snp

#miny = log10(min(df_counts_pvals$pval, df_counts_pvals_rh$pval))

### Plotting
plot1 <- ggplot(df_counts_pvals, aes(x = as.numeric(genomic_position), y = pval)) + geom_line() +
  scale_y_continuous(trans='log10')
plot1

plot2 <- ggplot(df_counts_pvals_rh, aes(x = as.numeric(genomic_position), y = pval)) + geom_line() +
  scale_y_continuous(trans='log10')
plot2

p1 <- plot1 / plot2 
p1 

### get out imputed genotypes
imputed_geno <- b[[3]]
imputed_geno_df <- as.data.frame(as.matrix(imputed_geno))
rownames(imputed_geno_df) <- imputed_geno_df$pos
imputed_geno_df <- imputed_geno_df[, 3:length(imputed_geno_df)]
colnames(imputed_geno_df) <- c(sprintf("g_%d", seq(1:ncol(imputed_geno_df))))

# gamete out imputed haplotypes
gam_haps <- b[[2]]
gam_haps_df <- as.data.frame(as.matrix(gam_haps))
rownames(gam_haps_df) <- gam_haps_df$pos
gam_haps_df <- gam_haps_df[, 3:length(gam_haps_df)]
colnames(gam_haps_df) <- c(sprintf("g_%d", seq(1:ncol(gam_haps_df))))

# calculate transmission rate
transmission_rates <- c()
for (i in seq(1: length(gam_haps_df[,1]))){
  transmission_rates <- c(transmission_rates, (length(which(gam_haps_df[i, ] == "h1"))/length(gam_haps_df[i, ])))
}

# Plot transmission rate
df_TR <- data.frame(pos = seq(1: length(gam_haps_df[,1])), rate = transmission_rates)
ggplot(data = df_TR, aes(x = pos, y = rate)) + geom_point() + ylim(0,1)

# Add h1 and h2 counts to df
#df_counts_pvals <- cbind(df_counts_pvals, h1_count = h1_count, h2_count = h2_count)

plt_me_df_counts_pvals <- cbind(df_counts_pvals, method = 'true')
plt_me_df_counts_pvals_rh <- cbind(df_counts_pvals_rh, method = 'imputed')
plt_me <- rbind(plt_me_df_counts_pvals, plt_me_df_counts_pvals_rh) %>%
  mutate(method = factor(method, levels = c("true", "imputed"))) %>%
  mutate(negative_log10_pval = -1 * log(pval, base = 10))

# Plot
# p3 <-  ggplot(df_counts_pvals, aes(x = as.numeric(genomic_position), y = pval, color = (pmax(h1_count,h2_count) / (h2_count + h1_count))*100)) + 
#   geom_point(size = 0.5) +
#   scale_color_viridis_c() +
#   scale_y_continuous(trans='log10') + 
#   geom_hline(yintercept = p_cutoff, col = 'blue') +
#   xlab("Genomic Position") + ylab("log10(p-value)") + 
#   theme_bw() + theme(panel.grid = element_blank()) +
#   labs(color = "Transmission Rate")

p3 <-  ggplot(plt_me, aes(x = 1000*as.numeric(genomic_position), y = negative_log10_pval, color = (pmax(h1_count,h2_count) / (h2_count + h1_count))*100)) +
  geom_vline(xintercept = 1000*killer_snp) +
  geom_point(size = 0.5) +
  scale_color_viridis_c() +
  facet_wrap('method', nrow = 2, labeller = as_labeller(c(`true` = "Simulated ground truth", `imputed` = 'Imputed'))) + 
  geom_hline(yintercept = -1 * log(p_cutoff, base = 10), col = 'blue') +
  xlab("Genomic Position") + 
  ylab(bquote(-log[10](italic(p)))) +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(color = "Transmission Rate") +
  theme(legend.key.height = unit(1.5, "cm"))
  


p3

p4 <-  ggplot(plt_me, aes(x = 1000*as.numeric(genomic_position), y = negative_log10_pval, color = (pmax(h1_count,h2_count) / (h2_count + h1_count))*100)) +
  #geom_vline(xintercept = 1000*killer_snp) +
  geom_point(size = 0.5) +
  scale_color_viridis_c() +
  facet_wrap('method', nrow = 2, labeller = as_labeller(c(`true` = "Simulated ground truth", `imputed` = 'Imputed'))) + 
  geom_hline(yintercept = -1 * log(p_cutoff, base = 10), col = 'blue') +
  xlab("Genomic Position") + 
  ylab(bquote(-log[10](italic(p)))) +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(color = "Transmission Rate") +
  theme(legend.key.height = unit(1.5, "cm"))



p4
                    
#true_gametes_df_mapped <- data.frame(lapply(true_gametes_df, function(x) x + 1))
#true_gametes_df_mapped <- data.frame(lapply(true_gametes_df_mapped, function(x) paste0("h", x)))
                                            
#length(which(imputed_geno_df[killer_snp,] != true_gametes_df[killer_snp,]))
lens <- c()
for (i in 1:n_snp){
  lens <- c(lens, length(which(imputed_geno_df[i,] != true_gametes_df[i,])))  
}

ggplot(data = data.frame(pos = 1:n_snp, lens = lens), aes(x =pos , y= lens)) + geom_point()


