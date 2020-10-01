#!/usr/bin/env Rscript

library(tidyr)
library(dplyr)

args <-commandArgs(TRUE)

table <- read.table(args[1], header = TRUE, sep = "\t", dec = ".")
new_table <- table %>% pivot_wider(names_from = cell, values_from = gt) %>% arrange(pos)
write.csv(new_table, file=args[2], row.names=FALSE, col.names=TRUE, quote=FALSE)