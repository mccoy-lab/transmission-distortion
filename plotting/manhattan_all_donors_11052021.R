library(data.table)
library(tidyverse)
library(scattermore)
library(pbapply)
library(viridisLite)
library(ggrastr)
library(ggnewscale)
library(ggplot2)

#flist_1 <- list.files("/scratch/groups/rmccoy22/kweave23/sc_transmission_distortion/run_full_filter_20210826/", 
#                      "_pval.csv",
#                      full.names = TRUE,
#                     recursive = TRUE)

#flist_2 <- list.files("/scratch/groups/rmccoy22/kweave23/sc_transmission_distortion/infertile_donors/", 
#                      "_pval.csv",
#                      full.names = TRUE,
#                      recursive = TRUE)

#flist_3 <- "~/scratch/sperm-seq/nc9ab_12_pval.csv"

#file_list <- c(flist_1, flist_2, flist_3)


file_list <- list.files("/scratch/groups/rmccoy22/kweave23/sc_transmission_distortion/runworkflow_20211105/",
                      "_pval.csv",
                      full.names = TRUE,
                      recursive = TRUE)

read_pvals <- function(file_name) {
  fread(file_name) %>%
    .[, file_name := basename(file_name)] %>%
    return()
}

dt <- rbindlist(pblapply(1:length(file_list), function(x) read_pvals(file_list[x])))

dt[, c("subject", "chrom", "drop") := tstrsplit(file_name, "_", fixed = TRUE)]
dt[, file_name := NULL]
dt[, drop := NULL]

dt$chrom <- factor(dt$chrom, levels = 1:22)


dt_sparse <- rbind(dt[pval < 0.01], dt[pval >= 0.01][sample(1:nrow(dt[pval >= 0.01]), 1000000)])
#fwrite(dt_sparse, file = "~/Dropbox/papers/2021_spermseq/td_sparse.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
fwrite(dt_sparse, file = "~/work/scarios1/transmission-distortion/td_sparse_11052021.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
# plot on local 

#dt_sparse <- fread("~/Dropbox/papers/2021_spermseq/td_sparse.txt")
dt_sparse <- fread("/Users/saracarioscia/mccoy-lab/transmission-distortion/plotting/td_sparse_11052021.txt")

r <- ggplot(data = dt_sparse,
            aes(x = genomic_position, y = subject, color = -log10(pval))) +
  geom_tile() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        axis.text.x = element_blank()) +
  scale_color_viridis_c(limits = c(0, 8)) +
  geom_hline(yintercept = 0.000000177703221, linetype="dashed", color="grey") +
  facet_grid(. ~ chrom, space = "free_x", scale = "free_x") +
  xlab("Genomic Position") +
  ylab("Subject")

#ggsave("~/tile.pdf", units = 'in', width = 8, height = 3, r)
ggsave("~/mccoy-lab/transmission-distortion/plotting/tile_11052021.pdf", units = 'in', width = 8, height = 3, r)

dt <- fread("/Users/saracarioscia/mccoy-lab/transmission-distortion/plotting/td_dt_11052021.txt")

setorder(dt, pval)
dt[, expected_p := .I / .N]
q <- ggplot(data = rbind(dt[pval < 0.01], dt[pval >= 0.01][sample(1:nrow(dt[pval >= 0.01]), 50000)]), 
            aes(x = -log10(expected_p), y = -log10(pval))) +
  geom_point_rast() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab(expression(Expected -log[10](italic(p)))) + 
  ylab(expression(Observed -log[10](italic(p)))) +
  xlim(0, 8) +
  ylim(0, 8) +
  geom_abline(slope = 1, lty = "dashed", color = "gray80")

ggsave('~/mccoy-lab/transmission-distortion/plotting/qq_11052021.pdf', units = 'in', width = 3, height = 3, q)

###

palette <- c("#9cba3f",
             "#9353c9",
             "#52b54b",
             "#d373dc",
             "#66be87",
             "#d03a93",
             "#3a7f4b",
             "#5670de",
             "#c9aa46",
             "#7358a7",
             "#e08a3c",
             "#466fb0",
             "#d4502d",
             "#3fbfbc",
             "#d64157",
             "#55a5db",
             "#a65035",
             "#9997e2",
             "#748438",
             "#9d4d8a",
             "#987033",
             "#d28bc5",
             "#e1917b",
             "#a2485d",
             "#e46f96")

subjects <- unique(dt_sparse$subject)

p <- ggplot(data = dt_sparse) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none") +
  facet_grid(. ~ chrom, space = "free_x", scale = "free_x") +
  xlab("Genomic Position") +
  ylab(bquote(-log[10](italic(p)))) +
  ylim(0, 8)

for (i in 1:25) {
  p <- p +
    new_scale("color") +
    geom_point_rast(data = dt_sparse[subject == subjects[i]], 
                    aes(x = genomic_position, y = -log10(pval), color = -log10(pval)),
                    size = 1) + 
    scale_color_gradient(low = "gray50", high = palette[i])
}
p <- p + geom_hline(yintercept = -log10(1.777703221e-7), color = "gray50", lty = "dashed")

#ggsave("~/Dropbox/papers/2021_spermseq/manhattan.pdf", units = 'in', width = 8, height = 3, p)
ggsave("/Users/saracarioscia/mccoy-lab/transmission-distortion/plotting/manhattan_all_donors_11052021.pdf", units = 'in', width = 8, height = 3, p)
