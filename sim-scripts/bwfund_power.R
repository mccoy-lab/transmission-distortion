library(tidyverse)
library(ggplot2)
library(data.table)

sim_binom_test <- function(n, bias, p = 0.5) {
  # simulate 1 coin flip n times with the specified bias
  coin <- rbinom(1, n, bias)
  # run a binomial test on the simulated data for the specified p
  btest <- binom.test(coin, n, p)
  # return the p-value of this test
  btest$p.value
}
rate_td <- seq(.5,.6,.01) # different rates of distortion
alpha <- 0.05 # this does not always have to be 0.05
# loop through different sample sizes 
sample_size <- seq(500, 10000, 500)
# empty data frame to save the info to                         
full_list <- data.frame(matrix(ncol = 3, nrow = 0))
# iterate through the levels of TD, and run the simulation 1000 times 
for(j in rate_td) {
  for(i in sample_size) {
    my_reps <- replicate(1e3, sim_binom_test(n=i,bias=j)) # this holds the p values
    power <- (mean(my_reps < alpha)) # this is power
    df <- data.frame(i,j,power)
    #print(df)
    full_list <- rbind(full_list, df)  
  }
}
td_rate <- as.factor(full_list$j)

ggplot(data = full_list, aes(x = i, y = factor(j),  fill = power)) +
  geom_tile() +
  xlab("Number of sperm") +
  ylab("Transmission rate") +
  scale_fill_viridis_c(name = "Power") + 
  theme(panel.grid = element_blank(), panel.background = element_blank())
