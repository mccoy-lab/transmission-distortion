library(tidyverse)
library(data.table)
read_pvals <- function(sample_id, chrom) {
  message(paste(sample_id, chrom))
  tryCatch(
    {
      dt <- fread(paste0("/scratch/groups/rmccoy22/kweave23/sc_transmission_distortion/",
                         "table_with_uncorrected_pval_",
                         sample_id,
                         "_chr",
                         chrom,
                         ".csv"), 
                  header = TRUE) %>%
        .[, sample_id := sample_id] %>%
        .[, chrom := chrom]
    },
    error = function(e){
      dt <- data.table(V1 = NA, pval = NA, h1_count = NA, h2_count = NA, 
                       genomic_position = NA, sample_id = NA, chrom = NA)
    })
}
dt_pvals <- rbind(
  rbindlist(lapply(1:22, function(x) read_pvals("nc1abnov17", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc2absept17", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc3aboct17", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc4abnov17", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc6abcd", x)))
)
dt_pvals_downsample <- rbind(dt_pvals[pval < 0.01],
                             dt_pvals[pval >= 0.01] %>%
                               .[sample(1:nrow(.), 100000),])
pdf("~/td_manhattan.pdf", height = 8, width = 12)
ggplot(data = dt_pvals_downsample, 
       aes(x = genomic_position, y = -log10(pval))) +
  geom_point(size = 0.5) +
  theme_bw() +
  facet_grid(sample_id ~ factor(chrom), scales = "free_x", space = "free") + 
  theme(panel.spacing.x = unit(0, "lines"),
        axis.text.x = element_blank()) +
  xlab("Genomic Coordinate (bp)") +
  ylab(expression(-log[10](p)))
dev.off()