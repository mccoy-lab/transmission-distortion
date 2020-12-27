library(tidyverse)
library(pbapply)

test1 <- c("h1", "h1", "h1", "h1", "h1", NA,   "h2", NA,   "h1", NA,   NA,   "h2", "h2", "h2", "h2")
test2 <- c("h1", "h1", "h1", "h1", "h1", NA, NA, NA, NA, NA, NA, "h2", "h2", "h2", "h2")
test3 <- c("h1", "h1", "h1", NA, NA, "h2", NA, "h1", NA, NA, "h2", NA, "h1", "h1", "h1")
test4 <- c("h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1")
test5 <- c("h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", "h2", "h2", "h2", "h2", "h2")
test6 <- c("h1", "h1", "h1", "h1", "h1", "h2", "h2", "h2", "h2", "h2", "h2", "h1", "h1", "h1", "h1")
overall_test <- cbind(test1, test2, test3, test4, test5, test6)
#test <- c("h1", "h1", "h1", "h1", "h1", NA, "h1", "h1", NA, NA, "h1", "h1", NA, "h2" ) # this case will cause an error
#test <- c("h1", "h1", "h1", "h1", "h1", "h1", "h1", "h1", NA, NA, "h2", "h1", NA, "h2" ) #this will catch two of the 3 recombination spots


find_recomb_spots <- function(input_tib, x, identities) {
input_vec <- input_tib[,x]
ident <- identities[x]
message(ident)
single_rs <- FALSE

na_locs <- which(is.na(input_vec)) #which locations have NAs
num_na_locs <- sum(is.na(input_vec)) #how many NAs
rle_result <- rle(input_vec)

if ((num_na_locs == 0) & (length(rle_result$values) == 1)) { #no recombination spots
  recomb_spots <- as_tibble(cbind(ident, "None", "None"))
  colnames(recomb_spots) <- c("Ident", "Start Index", "End Index")
  return (recomb_spots)
}

if ((num_na_locs == 0) & (length(rle_result$values) > 1)) { #recombination spots, but no NA's buffering them
  recomb_spot_starts <- lapply(rle_result["lengths"], cumsum)$lengths[1:length(rle_result$values)-1]
  recomb_spot_ends <- recomb_spot_starts + 1
  recomb_spots <- as_tibble(cbind(ident, recomb_spot_starts, recomb_spot_ends))
  colnames(recomb_spots) <- c("Ident", "Start Index", "End Index")
  return (recomb_spots)
}

neighboring_difs <- diff(na_locs)

subsetted_na_locs <- c(na_locs[1])
for (i in 2:length(na_locs)){
  if (neighboring_difs[i-1] != 1){
    subsetted_na_locs <- cbind(subsetted_na_locs, na_locs[i]) #need this for when we have different number of NAs, especially more than 1 between the different haplotypes
  }
}


first_na_loc <- subsetted_na_locs[1] #first NA loc
first_na_loc_fi <- first_na_loc #copying variable to use it for indexing later
last_na_loc <- na_locs[length(na_locs)] #last NA loc
i <- 1 #start incremental counter

recomb_spots <- c() #initialize empty list

if ((last_na_loc - first_na_loc + 1) == num_na_locs) { #there's a single recombination spot
  single_rs <- TRUE
  if (input_vec[first_na_loc-1] != input_vec[last_na_loc + 1]) {
  recomb_spots <- rbind(recomb_spots, c(ident, first_na_loc - 1, last_na_loc + 1))
  } else {stop(paste0("Haplotypes surrounding first set of NAs match but appear to be a recombination spot at ", first_na_loc -1, " and ", last_na_loc + 1))}
} else {
  
  while ((last_na_loc - first_na_loc + 1) > num_na_locs) { #means that we must have some haplotypes in the middle or there must be more than one recombination spot
    first_hp_loc <- which(!is.na(input_vec[first_na_loc_fi : last_na_loc]))[i] #this will find the first haplotype location for each "recombination spot", unless it's the last, then it should return an NA
    if (!is.na(first_hp_loc)) { #for all but the last recombination spots
      if (input_vec[first_na_loc - 1] != input_vec[first_na_loc_fi + first_hp_loc - 1]) {#make sure haplotypes are different
      recomb_spots <- rbind(recomb_spots, c(ident, first_na_loc - 1, first_na_loc_fi  + first_hp_loc -1)) #add the indices surrounding the recombination spot
      } else {stop(paste0("Haplotypes surrounding NAs match but appear to be a recombination spot at ", first_na_loc -1, " and ", first_na_loc_fi + first_hp_loc - 1))}
  
      i <- i + 1 #add to incremental counter
      first_na_loc <- subsetted_na_locs[i] #find new first na for next recombination spot
      num_na_locs <- sum(is.na(input_vec[(first_na_loc + first_hp_loc-1):length(input_vec)])) #find new number of remaining/unaccounted for NAs
            
      }
    else { #reached last recombination spots
      break
    }
  }
}
  
if (((last_na_loc - first_na_loc + 1) > 0) & !single_rs) { #add last recombination spots
  if (input_vec[first_na_loc - 1] != input_vec[last_na_loc + 1]) { #make sure haplotypes are different
  recomb_spots <- rbind(recomb_spots, c(ident, first_na_loc - 1, last_na_loc + 1))
  } else {stop(paste0("Haplotypes surrounding last recombination spot NAs match but appear to be a recombination spot at ", first_na_loc -1, " and ", last_na_loc + 1))}
}

recomb_spots <- as_tibble(recomb_spots)
colnames(recomb_spots) <- c("Ident", "Start Index", "End Index")

return (recomb_spots)

}

recomb_spots_all <- do.call(rbind, pblapply(1:ncol(overall_test), function(x) find_recomb_spots(overall_test, x, c('test1', 'test2', 'test3', 'test4', 'test5', 'test6'))))