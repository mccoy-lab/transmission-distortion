# Compute number of samples necessary for a specified power using the TDT 

# load libraries
library(pbapply)

sim_mcnemar_test <- function(n, bias, p = 0.5) {
  # simulate 1 coin flip n times with the specified bias
  b <- rbinom(1, n, bias)
  # run a mcnemar test on the simulated data for the specified p
  c <- n - b 
  mcnemar_test <- pchisq((b-c)^2/(b+c), df = 1, lower.tail = FALSE)
  # return the p-value of this test
  return(mcnemar_test)
}

compute_power <- function(n_reps, n, bias, alpha, p = 0.5) {
  my_reps <- replicate(n_reps, sim_mcnemar_test(n = n, bias = bias))
  power <- (mean(my_reps < alpha)) # this is power
  return(power)
}


rate_td <- seq(.5,.7,.01) # different rates of distortion
sample_size <- seq(500, 3500, 50)
alpha <- c(0.05, 10e-7)
parameter_grid <- data.table(expand.grid(rate_td, sample_size, alpha)) %>%
  setnames(., c("rate_td", "sample_size", "alpha"))

power_results <- pbmapply(compute_power, n_reps = 1e3, p = 0.5,
                          n = parameter_grid$sample_size, 
                          bias = parameter_grid$rate_td, 
                          alpha = parameter_grid$alpha)



parameter_grid[, power := power_results]

panel_a <- ggplot(data = parameter_grid[alpha == 0.05], aes(x = sample_size, y = rate_td, fill = power)) +
  geom_tile() +
  #geom_vline(xintercept = 2362, colour = "red") + 
  #geom_segment(aes(x = 2362 , y = 0.5, xend = 2362, yend = 0.7), colour="red") +
  xlab("Number of informative transmissions") +
  ylab("Transmission rate") +
  scale_fill_viridis_c(name = "Power") + 
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text=element_text(size=20), axis.title=element_text(size=26), 
        legend.key.size = unit(3, 'cm'), legend.title = element_text(size=26), legend.text = element_text(size=20))

panel_b <- ggplot(data = parameter_grid[alpha == 10e-7], aes(x = sample_size, y = rate_td, fill = power)) +
  geom_tile() +
  #geom_vline(xintercept = 2362, colour = "red") + 
  #geom_segment(aes(x = 2362 , y = 0.5, xend = 2362, yend = 0.7), colour="red") +
  xlab("Number of informative transmissions") +
  ylab("Transmission rate") +
  scale_fill_viridis_c(name = "Power") + 
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text=element_text(size=20), axis.title=element_text(size=26), 
        legend.key.size = unit(3, 'cm'), legend.title = element_text(size=26), legend.text = element_text(size=20))
