# calculating power to detect TD in single-cell sequencing 

library(tidyverse)
library(ggplot2)
library(data.table)
library(cowplot)
library(pbapply)

sim_binom_test <- function(n, bias, p = 0.5) {
  # simulate 1 coin flip n times with the specified bias
  coin <- rbinom(1, n, bias)
  # run a binomial test on the simulated data for the specified p
  btest <- binom.test(coin, n, p)
  # return the p-value of this test
  return(btest$p.value)
}

compute_power <- function(n_reps, n, bias, alpha, p = 0.5) {
  my_reps <- replicate(n_reps, sim_binom_test(n = n, bias = bias))
  power <- (mean(my_reps < alpha)) # this is power
  return(power)
}


rate_td <- seq(.5,.7,.01) # different rates of distortion
sample_size <- seq(500, 3500, 50)
alpha <- c(0.05, 1.78e-7)
parameter_grid <- data.table(expand.grid(rate_td, sample_size, alpha)) %>%
  setnames(., c("rate_td", "sample_size", "alpha"))

power_results <- pbmapply(compute_power, n_reps = 1e3, p = 0.5,
                          n = parameter_grid$sample_size, 
                          bias = parameter_grid$rate_td, 
                          alpha = parameter_grid$alpha)

parameter_grid[, power := power_results]

panel_a <- ggplot(data = parameter_grid[alpha == 0.05], aes(x = sample_size, y = rate_td, fill = power)) +
  geom_tile() +
  geom_vline(xintercept = 1711, colour = "red") + 
  #geom_segment(aes(x = 1711 , y = 0.5, xend = 1711, yend = 0.7), colour="red") +
  xlab("Number of gametes") +
  ylab("Transmission rate") +
  scale_fill_viridis_c(name = "Power") + 
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text=element_text(size=20), axis.title=element_text(size=26), 
        legend.key.size = unit(3, 'cm'), legend.title = element_text(size=26), legend.text = element_text(size=20)) 

panel_b <- ggplot(data = parameter_grid[alpha == 1.78e-7], aes(x = sample_size, y = rate_td, fill = power)) +
  geom_tile() +
  geom_vline(xintercept = 1711, colour = "red") + 
  #geom_segment(aes(x = 1711 , y = 0.5, xend = 1711, yend = 0.7), colour="red") +
  xlab("Number of gametes") +
  ylab("Transmission rate") +
  scale_fill_viridis_c(name = "Power") + 
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text=element_text(size=20), axis.title=element_text(size=26), 
        legend.key.size = unit(3, 'cm'), legend.title = element_text(size=26), legend.text = element_text(size=20))

#plot_grid(panel_a, panel_b, nrow = 1, labels = c("A.", "B."))


# try to simulate based on sperm count from a few individuals 
panel_c <- ggplot(data = parameter_grid[alpha == 0.05], aes(x = sample_size, y = power, color = rate_td)) +
  geom_vline(xintercept = 1711, colour = "red") + 
  #geom_segment(aes(x = 1711 , y = 0, xend = 1711, yend = 1), colour="red") +
  geom_point(aes(colour = rate_td)) +
  scale_y_continuous(breaks=seq(0,1.00,0.1)) +
  xlab("Number of gametes") +
  ylab("Power") +
  scale_colour_viridis_c(name = "Transmission rate") + 
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.text=element_text(size=20), axis.title=element_text(size=26), 
        legend.key.size = unit(3, 'cm'), legend.title = element_text(size=26), legend.text = element_text(size=20))


panel_d <- ggplot(data = parameter_grid[alpha == 1.78e-7], aes(x = sample_size, y = power, color = rate_td)) +
  geom_vline(xintercept = 1711, colour = "red") + 
  #geom_segment(aes(x = 1711 , y = 0, xend = 1711, yend = 1), colour="red") +
  geom_point(aes(colour = rate_td)) +
  scale_y_continuous(breaks=seq(0,1.00,0.1)) +
  xlab("Number of gametes") +
  ylab("Power") +
  scale_colour_viridis_c(name = "Transmission rate") + 
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text=element_text(size=20), axis.title=element_text(size=26), 
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        legend.key.size = unit(3, 'cm'), legend.title = element_text(size=26), legend.text = element_text(size=20))
