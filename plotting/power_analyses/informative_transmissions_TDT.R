# count informative transmissions
# using supplemental data from Meyer et al. 2012 

# remove headers on command line: sed -i -e 1,13d TableS1_no_headers.txt
library(scales) # to force y axis labels
library(dplyr)

supplemental_table_S1 <- read.table("/Users/saracarioscia/mccoy-lab/transmission-distortion/TableS1_no_headers.txt", header = TRUE, sep = " ", fill = TRUE)
names(supplemental_table_S1) <- c("SNP", "CHR", "BP", "AGRE_A1", "AGRE_A2", "AGRE_T", "AGRE_U", "AGRE_P", "AGRE_T_PAT", "AGRE_U_PAT", "AGRE_P_PAT", "AGRE_T_MAT", "AGRE_U_MAT", "AGRE_P_MAT", "HUTT_A1", "HUTT_A2", "HUTT_T", "HUTT_U", "HUTT_P", "HUTT_T_PAT", "HUTT_U_PAT", "HUTT_P_PAT", "HUTT_T_MAT", "HUTT_U_MAT", "HUTT_P_MAT")

# 

# get count of each SNP in AGRE 
supplemental_table_S1$AGRE_T <- supplemental_table_S1$AGRE_T %>% as.numeric() 
supplemental_table_S1$AGRE_U <- supplemental_table_S1$AGRE_U %>% as.numeric()
supplemental_table_S1$AGRE_transmissions <- supplemental_table_S1$AGRE_T + supplemental_table_S1$AGRE_U
# for X axes number of transmissions: min 200, max 1684

# get count of each SNP in HUTT 
supplemental_table_S1$HUTT_T <- supplemental_table_S1$HUTT_T %>% as.numeric() 
supplemental_table_S1$HUTT_U <- supplemental_table_S1$HUTT_U %>% as.numeric() 
supplemental_table_S1$HUTT_transmissions <- supplemental_table_S1$HUTT_T + supplemental_table_S1$HUTT_U
# for X axes number of transmissions: min 50, max 1319
# remove those with fewer than 200 transmissions, to match filtering on AGRE 
HUTT_transmissions_filtered <- supplemental_table_S1[supplemental_table_S1$HUTT_transmissions > 200,]

# plot the two histograms 
panel_AGRE <- ggplot(data = supplemental_table_S1, aes(x = AGRE_transmissions)) + 
  geom_histogram() + 
  theme_bw() + 
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text=element_text(size=20), axis.title=element_text(size=26), 
        legend.key.size = unit(3, 'cm'), legend.title = element_text(size=26)) + 
  xlab("Number of informative transmissions") +
  ylab("Number of SNPs") +
  scale_y_continuous(labels = comma) + 
  xlim(0, 1700) 

panel_HUTT <- ggplot(data = supplemental_table_S1, aes(x = HUTT_transmissions)) + 
  geom_histogram() + 
  theme_bw() + 
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text=element_text(size=20), axis.title=element_text(size=26), 
        legend.key.size = unit(3, 'cm'), legend.title = element_text(size=26)) + 
  xlab("Number of informative transmissions") +
  ylab("Number of SNPs") +
  scale_y_continuous(labels = comma) + 
  xlim(0, 1700)
