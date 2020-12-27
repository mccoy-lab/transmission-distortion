library(tidyverse)

test <- c("h1", "h1", "h1", "h1", "h1", NA,   "h2", NA,   "h1", NA,   NA,   "h2", "h2", "h2", "h2")
test <- c("h1", "h1", "h1", "h1", "h1", NA, NA, NA, NA, NA, NA, "h2", "h2", "h2", "h2")
test <- c("h1", "h1", "h1", "h1", "h1", NA, NA, "h2", NA, "h1", NA, NA, "h2", NA, "h1", "h1", "h1")
test <- c("h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1")
test <- c("h1", "h1", "h1", "h1", "h1", "h2", "h2", "h2", "h2", "h2")
test <- c("h1", "h1", "h1", "h1", "h1", "h2", "h2", "h2", "h2", "h2", "h1", "h1", "h1")
test <- c("h1", "h1", "h1", "h1", "h1", NA, "h1", "h1", NA, NA, "h1", "h1", NA, "h2" ) # this case will introduce two recombination spots that don't exist


find_recomb_spots <- function(test_vec) {
single_rs <- FALSE

na_locs <- which(is.na(test_vec)) #which locations have NAs
num_na_locs <- sum(is.na(test_vec)) #how many NAs
rle_result <- rle(test_vec)

if ((num_na_locs == 0) & (length(rle_result$values) == 1)) { #no recombination spots
  return (as_tibble(c()))
}

if ((num_na_locs == 0) & (length(rle_result$values) > 1)) { #recombination spots, but no NA's buffering them
  recomb_spot_starts <- lapply(rle_result["lengths"], cumsum)$lengths[1:length(rle_result$values)-1]
  recomb_spot_ends <- recomb_spot_starts + 1
  recomb_spots <- as_tibble(cbind(recomb_spot_starts, recomb_spot_ends))
  colnames(recomb_spots) <- c("Start Index", "End Index")
  return (recomb_spots)
}

neighboring_difs <- diff(na_locs)

subsetted_na_locs <- c(na_locs[1])
for (i in 2:length(na_locs)){
  if (neighboring_difs[i-1] != 1){
    subsetted_na_locs <- cbind(subsetted_na_locs, na_locs[i])
  }
}


first_na_loc <- subsetted_na_locs[1] #first NA loc
first_na_loc_fi <- first_na_loc #copying variable to use it for indexing later
last_na_loc <- na_locs[length(na_locs)] #last NA loc
i <- 1 #start incremental counter

recomb_spots <- c() #initialize empty list

if ((last_na_loc - first_na_loc + 1) == num_na_locs) { #there's a single recombination spot
  single_rs <- TRUE
  recomb_spots <- rbind(recomb_spots, c(first_na_loc - 1, last_na_loc + 1)) 
} else {
  
  while ((last_na_loc - first_na_loc + 1) > num_na_locs) { #means that we must have some haplotypes in the middle or there must be more than one recombination spot
    first_hp_loc <- which(!is.na(test_vec[first_na_loc_fi : last_na_loc]))[i] #this will find the first haplotype location for each "recombination spot", unless it's the last, then it should return an NA
    if (!is.na(first_hp_loc)) { #for all but the last recombination spots
      recomb_spots <- rbind(recomb_spots, c(first_na_loc - 1, first_na_loc_fi  + first_hp_loc -1)) #add the indices surrounding the recombination spot
  
  
      i <- i + 1 #add to incremental counter
      first_na_loc <- subsetted_na_locs[i] #find new first na for next recombination spot
      num_na_locs <- sum(is.na(test_vec[(first_na_loc + first_hp_loc-1):length(test_vec)])) #find new number of remaining/unaccounted for NAs
            
      }
    else { #reached last recombination spots
      break
    }
  }
}
  
if (((last_na_loc - first_na_loc + 1) > 0) & !single_rs) { #add last recombination spots
  message('got here')
  recomb_spots <- rbind(recomb_spots, c(first_na_loc - 1, last_na_loc + 1))
}

recomb_spots <- as_tibble(recomb_spots)
colnames(recomb_spots) <- c("Start Index", "End Index")
}
