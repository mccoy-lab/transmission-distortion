# Plot randomly sampled sperm, after assigning stretches to each haplotype

# filled_sperm is the output of other steps in assign_haplotypes 
outcome_new_5 <- as.data.frame(filled_sperm)
outcome_new_5 <- outcome_new_5[,sample(1:ncol(outcome_new_5),20, replace=FALSE)]
rownames(outcome_new_5) <- positions

# write filled sperm to its own file for further use
write.csv(outcome_new_5,"/Users/saracarioscia/mccoy-lab/transmission-distortion/raw_data_for_tests/NC17chr6_filledsperm.csv", row.names = TRUE)
write.csv(outcome_new_5,"/Users/saracarioscia/mccoy-lab/transmission-distortion/raw_data_for_tests/NC17chr6_filledsperm_20sample.csv", row.names = TRUE)

# reformat filled_sperm
outcome_new6 <- outcome_new_5 %>% 
  tibble::rownames_to_column('location') %>%
  transform(location = as.numeric(location)) %>%
  pivot_longer(-'location', names_to = 'chromosome', values_to = 'haplotype')

# grab the cell IDs of the chromosomes randomly chosen above 
inferred_cell_IDs <- unique(outcome_new6$chromosome)

library(ggplot2)
# Plot whole length of chromosome 
outcome_new6 %>% ggplot(aes(x=location, y = chromosome)) +
  geom_point(aes(color=haplotype)) +
  ggtitle("donor 17 chromosome 6") + 
  theme_minimal()

# plot only the region of the chromosome with the peak 
outcome_new6 %>% ggplot(aes(x=location, y = chromosome)) +
  geom_point(aes(color=haplotype)) +
  xlim(158135000,170612671) + 
  ggtitle("donor 17 chromosome 6") + 
  theme_minimal()