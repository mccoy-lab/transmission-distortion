# Comparing depth of bam file on full chromosome, within peak, and outside of peak 
# Uses bam file that is the output of https://github.com/mccoy-lab/transmission-distortion/blob/master/shell-scripts/process_bam_files.txt

# Load libraries
library(tidyr)
library(ggplot2)

# Read in depth file 
SRR10140446_chr10_depth <- read.table("/Users/saracarioscia/mccoy-lab/transmission-distortion/SRR10140446_chr10_depth.txt", sep = "\t", header = FALSE)
SRR10140446_chr10_depth_only <- SRR10140446_chr10_depth[,3] %>% as.data.frame()

# Plot 
ggplot(SRR10140446_chr10_depth_only,aes(x=.) ) + geom_histogram(binwidth = 1) 

# Cut tail of plots
ggplot(SRR10140446_chr10_depth_only,aes(x=.) ) + geom_histogram(binwidth = 1) + xlim(0,150) + xlab("Depth") + geom_vline(mean) + geom_vline(sd) + geom_point(x=28.094)

mean <- mean(SRR10140446_chr10_depth_only[,1]) # 20.17671
sd <- sd(SRR10140446_chr10_depth_only[,1]) # 17.14197
cutoff <- mean + sd # 37.31868
