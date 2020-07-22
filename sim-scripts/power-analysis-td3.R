# Power analysis for TD project: simulate different levels of transmission distortion, see how often we can detect it (# that way if we have a negative result (i.e. can’t find any TD), we know how big of effects we can rule out)

library(data.table)
library(ggplot2)

# simulate data using rbinom() with different amounts of transmission distortion, and then see if you can detect it using a binomial test, while varying the sample size of sperm
# e.g., https://psyteachr.github.io/msc-data-skills/sim.html#binomial, and the function `sim_binomial_test` for an example of how to do this

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

# loop through different sample sizes, we have 974-2,274 gametes per donor, 20 donors
sample_size <- c(1e2, 300, 500, 700, 900, 1e3, 1e4)

# empty data frame to save the info to                         
full_list <- data.frame(matrix(ncol = 7, nrow = 0))

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
ggplot(data = full_list, aes(x = i, y = power,  color = factor(j))) +
  geom_point() +
  geom_line() +
  theme_bw() +
  scale_x_log10() +
  xlab("Sample size") +
  ylab("Power") +
  scale_color_manual(values = unname(palette.colors(11, "Alphabet")), name = "Transmission bias")



# info about -- mean(my_reps < alpha)
# calculate the power of your analysis by checking the proportion of your simulated analyses that have a p-value 
# less than your alpha (the probability of rejecting the null hypothesis when the null hypothesis is true)
# if power prints out 0.70, it means you have 70% power to detect the difference
