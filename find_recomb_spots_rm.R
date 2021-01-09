find_recomb_spots <- function(input_matrix, x, identities, genomic_positions){
  ident <- identities[x]
  input_tibble <- input_matrix[, x] %>%
    mutate(., index = row_number()) %>%
    mutate(., positions = genomic_positions)
  complete_cases_tibble <- input_tibble[complete.cases(input_tibble),]
  input_vec <- as.factor(complete_cases_tibble[[1]])
  switch_indices <- which(input_vec[-1] != input_vec[-length(input_vec)])
  switch_indices_input <- complete_cases_tibble[switch_indices,]$index
  crossover_start <- input_tibble[switch_indices_input,]$positions
  rev_input_tibble <- arrange(input_tibble, -index) %>%
    mutate(., index = row_number())
  complete_cases_rev_tibble <- rev_input_tibble[complete.cases(rev_input_tibble),]
  rev_input_vec <- as.factor(complete_cases_rev_tibble[[1]])
  rev_switch_indices <- which(rev_input_vec[-1] != rev_input_vec[-length(rev_input_vec)])
  rev_switch_indices_input <- complete_cases_rev_tibble[rev_switch_indices,]$index
  crossover_end <- rev(rev_input_tibble[rev_switch_indices_input,]$positions)
  recomb_spots <- tibble(Ident = ident, Genomic_start = crossover_start, Genomic_end = crossover_end)
  return(recomb_spots)
}
test_tibble <- tibble(sperm_1 = c(rep("h1", 5), rep(NA, 5), rep("h2", 5), rep(NA, 10), rep("h2", 50), rep("h1", 25)), sperm_2 = rep("h1", 100))
test_tibble2 <- tibble(sperm_1 = c("h1", "h1", "h1", "h1", "h1", NA,   "h2", NA,   "h1", NA,   NA,   "h2", "h2", "h2", "h2"), 
                       sperm_2 = c("h1", "h1", "h1", "h1", "h1", NA, NA, NA, NA, NA, NA, "h2", "h2", "h2", "h2"),
                       sperm_3 = c("h1", "h1", "h1", NA, NA, "h2", NA, "h1", NA, NA, "h2", NA, "h1", "h1", "h1"),
                       sperm_4 = c("h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1"),
                       sperm_5 = c("h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h2", "h2", "h2", "h2", "h2"),
                       sperm_6 = c("h1", "h1", "h1", "h1", "h1", "h2", "h2", "h2", "h2", "h2", "h2", "h1", "h1", "h1", "h1"),
                       sperm_7 = c("h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", NA, NA, "h2", "h1", NA, "h2", "h2"),
                       sperm_8 = c("h1", "h1", "h1", "h1", "h1", NA, "h1", "h1", NA, NA, "h1", "h1", NA, "h2", "h2"))

find_recomb_spots(test_tibble, 
                  1, 
                  "test_sperm", 
                  1:100)
find_recomb_spots(test_tibble, 
                  2, 
                  "test_sperm", 
                  1:100)

find_recomb_spots(test_tibble2, 1, "test_sperm", 1:15)
find_recomb_spots(test_tibble2, 2, "test_sperm", 1:15)
find_recomb_spots(test_tibble2, 3, "test_sperm", 1:15)
find_recomb_spots(test_tibble2, 4, "test_sperm", 1:15)
find_recomb_spots(test_tibble2, 5, "test_sperm", 1:15)
find_recomb_spots(test_tibble2, 6, "test_sperm", 1:15)
find_recomb_spots(test_tibble2, 7, "test_sperm", 1:15)
find_recomb_spots(test_tibble2, 8, "test_sperm", 1:15)
