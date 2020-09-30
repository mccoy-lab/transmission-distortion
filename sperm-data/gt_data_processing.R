#install.packages("tidyr")
#install.packages("dplyr")
library(tidyr)
library(dplyr)

"
Sample information for benchmarking: 
pos gt
23366809 1 1
28759494 0 0 0 
28755346 0
"

#pwd -P --> /scratch/groups/rmccoy22/scarios1/transmission-distortion

test_table <- read.delim("~/work/scarios1/transmission-distortion/test_nc9_21.cellsbyrow.txt", header = TRUE, sep = "\t", dec = ".")
new_table <- test_table %>% pivot_wider(names_from = cell, values_from = gt) %>% arrange(pos)

write.csv(new_table, file="~/work/scarios1/transmission-distortion/sample_data.csv", row.names = FALSE, col.names = TRUE, quote = FALSE)