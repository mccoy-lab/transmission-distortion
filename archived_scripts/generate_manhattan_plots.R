## Script to generate manhattan plots of our td data with every donor and chromosome 

library(tidyverse)
library(data.table)

# Adjust filename process based on naming convention and script location
read_pvals <- function(sample_id, chrom) {
  message(paste(sample_id, chrom))
  tryCatch(
    {
      data <- fread(paste0("~/work/kweave23/sc_transmission_distortion/run_bell_20210302/", 
                           sample_id,
                           "/csv_out/",
                           sample_id, 
                           "_",
                           chrom,
                           "_pval.csv"), header = TRUE) %>%
        .[, sample_id := sample_id] %>%
        .[, chrom := chrom]
    },
    error = function(e){
      data <- data.table(pval = NA, h1_count = NA, h2_count = NA, 
                         genomic_position = NA, sample_id = NA, chrom = NA)
    })
}

# Apply function to all chromosomes in all donors 
dt_pvals <- rbind(
  rbindlist(lapply(1:22, function(x) read_pvals("nc1abnov17", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc2absept17", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc3aboct17", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc4abnov17", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc6abcd", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc8ab", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc9ab", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc10oldoil", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc11ab", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc12ab", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc13ab", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc14ab", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc15ab", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc16ab", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc17ab", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc18ab", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc22abcd", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc25abcd", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc26abcd", x))),
  rbindlist(lapply(1:22, function(x) read_pvals("nc27aboct17", x)))
)

# Take all p values below 0.01 and then sample from those above 0.01
dt_pvals_downsample <- rbind(dt_pvals[pval < 0.01],
                             dt_pvals[pval >= 0.01] %>%
                               .[sample(1:nrow(.), 100000),])

# Might want to change format of "pdf" generating if 
pdf("~/work/scarios1/transmission-distortion/td_manhattan03032.pdf", height = 8, width = 12)

# Have color bar be ratio of TD 
ggplot(data = dt_pvals_downsample, 
       aes(x = genomic_position, y = -log10(pval), color = (pmax(h1_count,h2_count) / (h2_count + h1_count))*100)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c() +
  theme_bw() +
  facet_grid(sample_id ~ factor(chrom), scales = "free_x", space = "free") + 
  theme(panel.spacing.x = unit(0, "lines"),
        axis.text.x = element_blank()) +
  xlab("Genomic Coordinate (bp)") +
  ylab(expression(-log[10](p)))
dev.off()