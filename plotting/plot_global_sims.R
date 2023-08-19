library(data.table)
library(tidyverse)

obs <- 0.499996

dt <- fread("~/Downloads/null_sim.txt") %>%
  setnames(., c("seed", "sim"))

# not all simulations finished - take first 500

ggplot(data = dt[1:500], aes(x = sim)) +
  geom_histogram(bins = 50, fill = "gray30") +
  xlim(0.5 - 0.00015, 0.5 + 0.00015) +
  geom_vline(xintercept = obs, color = "red") +
  geom_vline(xintercept = 0.5,
             color = "gray70",
             lty = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("Proportion shared alleles") +
  ylab("Number of simulations")

mean(dt[1:500]$sim >= obs)