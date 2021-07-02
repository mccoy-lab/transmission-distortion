# Plot peak with random sperm cells 

library(data.table)
library(tidyverse)
library(pbapply)
library(cowplot)
library(viridis)

# Input to this is the filled_sperm from the rhapsodi pipeline 
dt <- fread("~/Downloads/NC17chr10_filledsperm.csv", header = TRUE, stringsAsFactors = TRUE)
setnames(dt, "V1", "pos")
td <- function(data, row_index) {
  data_vector <- as.numeric(unname(data[row_index,]))[-1]
  data_vector <- data_vector[!is.na(data_vector)]
  test <- binom.test(x = sum(data_vector == 1), n = length(data_vector))
  return(data.table(est = test$estimate, ll = test$conf.int[1], ul = test$conf.int[2], p = test$p.value))
}
results <- pblapply(1:nrow(dt), function(x) td(dt, x))
results <- rbindlist(results)
results[, pos := dt$V1]
setorder(results, p)
order <- dt[pos == 78372560] %>%
  melt(id.vars = "pos") %>%
  setorder(., value) %>%
  .$variable
dt_peak <- dt[pos > 76700000 & pos < 80000000][, c(1, sample(2:ncol(dt), 100)), with = FALSE] %>%
  .[, pos_index := .I] %>%
  melt(id.vars = c("pos_index", "pos"))
dt_peak$variable <- factor(dt_peak$variable, levels = order)
setorder(results, pos)
results[, pos_index := .I]
m1 <- ggplot(data = results[pos > 76700000 & pos < 80000000], aes(x = pos_index, y = -log10(p), color = est)) +
  geom_point() +
  scale_color_viridis_c(name = "Transmission rate") +
  theme_minimal() +
  ylim(0, 12)
p1 <- ggplot(data = dt_peak, aes(x = pos_index, y = variable, fill = value)) +
  geom_tile() +
  theme(axis.text = element_blank(),
        panel.background = element_blank()) +
  ylab("Sperm cell") +
  xlab("SNP") +
  scale_fill_manual(values = c("#fb8072", "#80b1d3"), name = "")
plot_grid(m1, p1, nrow = 2, align = "hv", axis = "tblr")
m2 <- ggplot(data = results[pos > 78200000 & pos < 78600000], aes(x = pos_index, y = -log10(p), color = est)) +
  geom_point() +
  scale_color_viridis_c(name = "Transmission rate") +
  theme_minimal() +
  ylim(0, 12)
p2 <- ggplot(data = dt_peak[pos > 78200000 & pos < 78600000], aes(x = pos_index, y = variable, fill = value)) +
  geom_tile() +
  theme(axis.text = element_blank(),
        panel.background = element_blank()) +
  ylab("Sperm cell") +
  xlab("SNP") +
  scale_fill_manual(values = c("#fb8072", "#80b1d3"), name = "")
plot_grid(m2, p2, nrow = 2, align = "hv", axis = "tblr")