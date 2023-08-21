library(data.table)
library(tidyverse)
library(scattermore)
library(pbapply)
library(viridisLite)
library(ggplot2)

flist_1 <- list.files("/scratch/groups/rmccoy22/kweave23/sc_transmission_distortion/run_full_filter_20210826/", 
                      "_pval.csv",
                      full.names = TRUE,
                      recursive = TRUE)

flist_2 <- list.files("/scratch/groups/rmccoy22/kweave23/sc_transmission_distortion/infertile_donors/", 
                      "_pval.csv",
                      full.names = TRUE,
                      recursive = TRUE)


file_list <- c(flist_1, flist_2)

read_pvals <- function(file_name) {
  fread(file_name) %>%
    .[, file_name := basename(file_name)] %>%
    return()
}

dt <- rbindlist(pblapply(1:length(file_list), function(x) 
read_pvals(file_list[x])))

dt[, c("subject", "chrom", "drop") := tstrsplit(file_name, "_", fixed = TRUE)]
dt[, file_name := NULL]
dt[, drop := NULL]

dt$chrom <- factor(dt$chrom, levels = 1:22)

dt_sparse <- rbind(dt[pval < 0.01], dt[pval >= 0.01][sample(1:nrow(dt[pval >= 0.01]), 1000000)])

r <- ggplot(data = dt_sparse,
            aes(x = genomic_position, y = subject, color = -log10(pval))) +
  geom_tile() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        axis.text.x = element_blank()) +
  scale_color_viridis_c(limits = c(0, 8)) +
  facet_grid(. ~ chrom, space = "free_x", scale = "free_x") +
  xlab("Genomic Position") +
  ylab("Subject")

ggsave("./tile.pdf", units = 'in', width = 8, height = 3, r)

setorder(dt, pval)
dt[, expected_p := .I / .N]

print("good to before q")

q <- ggplot(data = rbind(dt[pval < 0.01], dt[pval >= 0.01][sample(1:nrow(dt[pval >= 0.01]), 50000)]), 
            aes(x = -log10(expected_p), y = -log10(pval))) + xlab(expression(paste("Hb", A[1][c]," (%)", sep=""))) +
  geom_point(size = 2) +
  xlab(expression(Expected -log[10](italic(p)))) + 
  ylab(expression(Observed -log[10](italic(p)))) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlim(0, 8) +
  ylim(0, 8) +
  geom_abline(slope = 1, lty = "dashed", color = "gray80")

ggsave('./qq.pdf', units = 'in', width = 3, height = 3, q)