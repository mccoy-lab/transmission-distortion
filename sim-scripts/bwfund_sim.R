library(data.table)
library(tidyverse)

# generate a heterozygous sample at 30,000 sites
hap <- data.frame(V1 = sample(c(0, 1), size = 30000, replace = TRUE))
hap$V2 <- 1 - hap$V1

# we want a new table with N_sperm columns and N_snp rows
rows <- nrow(hap)
print(rows)
M <- matrix(ncol = 10000, nrow = rows)
print(M)

# using x, assign recombination breakpoints randomly. Send columns to be columns in M 
# at some point, also import legend file so you can replace the 0 and 1 with the allele letters
recomb_spot = sample(1:(rows - 1), 10000, replace=T)
print(recomb_spot)

# the change starts / the breakpoint is in the i row, including the i row. (The i row already belongs to the second haplotype)
for(i in 1:10000) {
  haps <- c(1, 2)
  hap_choice = sample(haps, size = 1)
  if (hap_choice < 2) {
    hap_second <- 2 
  } else {
    hap_second <- 1
  }
  # take first column of hap; take the first recomb_spot[i] rows
  # randomly choose either haplotype 1 or haplotype 2 for the first and then the other for the second part 
  M[1:recomb_spot[i], i] = hap[1:recomb_spot[i], hap_choice]
  M[(recomb_spot[i] + 1):rows, i] = hap[(recomb_spot[i] + 1):rows, hap_second]
}
print(M)

p_vals <- apply(M, MARGIN = 1, FUN = function(x) binom.test(table(x))$p.value)
hist(p_vals, breaks = 100)
plot(-log10(p_vals), ylim = c(0, 8))

#

snp_of_interest <- 15000
soi_zero <- which(M[snp_of_interest,] == 0)
soi_one <- which(M[snp_of_interest,] == 1)
M_biased <- M[, c(sample(soi_zero, 300, replace = FALSE), sample(soi_one, 200, replace = FALSE))]
p_vals_biased <- apply(M_biased, MARGIN = 1, FUN = function(x) binom.test(table(x))$p.value)
hist(p_vals_biased, breaks = 100)
plot(-log10(p_vals_biased), ylim = c(0, 8))

dt <- data.table(p_val = p_vals_biased) %>%
  .[, index := .I]

ggplot(data = dt, aes(x = index, y = -log10(p_val))) +
  geom_point(size = 0.1) +
  geom_vline(xintercept = 15000, color = 'red', lty = "dashed") +
  theme_classic() +
  xlab("SNP position on chromosome") +
  theme(axis.text.x = element_blank()) +
  ylab(expression("-"*log[10]*"(p-value)"))



