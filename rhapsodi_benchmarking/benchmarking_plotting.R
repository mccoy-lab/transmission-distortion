library(ggplot2)
library(tidyverse)

files <- list.files(path="data/times_out", full.names=FALSE, recursive=FALSE)

times <- data.frame(matrix(ncol = 8, nrow = 0))
colNames_times <- c("num_gametes", "coverage", "rs", "user.self", "sys.self", "elapsed", "user.child", "sys.child")
colnames(times) <- colNames_times

for (f in files){
  fname_split <- strsplit(f, '_')
  n_gam <- as.integer(strsplit(fname_split[[1]][1], "g")[[1]][[2]]) 
  coverage <- as.numeric(strsplit(fname_split[[1]][2], "c")[[1]][[2]])
  rs <- as.numeric(strsplit(fname_split[[1]][3], "rs")[[1]][[2]])
  
  time_file <- read.csv(paste0("data/times_out/", f))
  
  meta <- data.frame(num_gametes = n_gam, coverage = coverage, rs = rs)
  t <- cbind(meta, time_file)
  
  times <- rbind(times, t)
}

pdf("rhapsodi_run_time.pdf", width = 10, height = 10)

ggplot(data=times, aes(x=coverage, y=(user.self + sys.self + user.child + sys.child), group = num_gametes, colour = factor(num_gametes))) + geom_point() + 
  scale_color_viridis_d(begin = 0,end = 0.96, name = "Number of Gametes") +
  theme_bw() + theme(panel.grid = element_blank(), panel.background = element_blank()) + xlab("Coverage") + ylab("Run Time (CPU Seconds)") + 
  theme(text = element_text(size = 15)) 

# Closing the graphical device
dev.off() 