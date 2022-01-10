
library(ggplot2)
library(tidyverse)
library(patchwork)
library(svglite)
library(viridis)

# Prevent scientific notation (for file names that have higher snp counts)
options(scipen = 999)

#########PHASING###############################################################################################

# Initialize empty DF for phasing analysis results
phasing_results <- data.frame(matrix(ncol=10, nrow=0))
colnames(phasing_results) = c('n_gametes', 'n_snp', 'cov', 'seqerr', 'avgr', 'rs', 'acc', 'com', 'lhs', 'ser')

# Read in phasing data
# Loop through gamete numbers 
for (gam_number in c(3,15,50,150)){
  # Loop through files, reading them in
  for (fileName in list.files(paste0("~/work/abortvi2/rhapsodi/command_line_results/out/", gam_number, "_gmt/SAP"))){
    file <- read.csv(paste0("~/work/abortvi2/rhapsodi/command_line_results/out/", gam_number, "_gmt/SAP/", fileName), header=FALSE)
    # For multi-row files, only get the last line. 
    if (nrow(file) != 1){
      file = file[nrow(file),]
    } 
    # Add to data frame 
    phasing_results[nrow(phasing_results) + 1, ] <-  file
  }
}

# Add file name column for merging with rhapsodi data
phasing_results <- phasing_results %>%
  mutate(phasing_results, fname = paste('g', n_gametes, 's', n_snp, 'c', cov, 'se', seqerr, 'r', avgr, 'rs', rs, sep="_"))
#########RECOMBINATION###########################################################################################

# Initialize empty DF for recombination analysis results
recombination_results <- data.frame(matrix(ncol=19, nrow=0))
colnames(recombination_results) = c("n_gametes","n_snp","cov","seqerr","avgr","rs","precision","recall","accuracy","specificity","fdr","fpr","f1","true_n","pred_n","tn","fn","tp","fp")

# Read in recombination data
# Loop through gamete numbers 
## ADD 150 ONCE THOSE ARE DONE
for (gam_number in c(3,15,50,150)){
  # Loop through files, reading them in
  for (fileName in list.files(paste0("~/work/abortvi2/rhapsodi/command_line_results/out/", gam_number, "_gmt/SAR"))){
    file <- read.csv(paste0("~/work/abortvi2/rhapsodi/command_line_results/out/", gam_number, "_gmt/SAR/", fileName), 
                     header=TRUE,
                     row.names = NULL)
    # For multi-row files, only get the last line. 
    if (nrow(file) != 1){
      file = file[nrow(file),]
    } 
    file <- subset( file, select = -row.names)
    recombination_results[nrow(recombination_results) + 1, ] <-  file
  }
}

# Add file name column for merging with rhapsodi data
recombination_results <- recombination_results %>%
  mutate(recombination_results, fname = paste('g', n_gametes, 's', n_snp, 'c', cov, 'se', seqerr, 'r', avgr, 'rs', rs, sep="_"))

#########GAMETE#IMPUTATION######################################################################################

# Initialize empty DF for gamete imputation analysis results
gam_imputation_results <- data.frame(matrix(ncol=11, nrow=0))
# If I don't want to average across gametes, uncomment the following line
#colnames(gam_imputation_results) = c("row.names","n_gametes","n_snp","cov","seqerr","avgr","rs","acc","com","lhs","ser")
colnames(gam_imputation_results) = c("n_gametes","n_snp","cov","seqerr","avgr","rs","acc","com","lhs","ser")

# This loop takes averages of all variables, which is ultimately what we want.
for (gam_number in c(3,15,50,150)){
  # Loop through files, reading them in
  for (fileName in list.files(paste0("~/work/abortvi2/rhapsodi/command_line_results/out/", gam_number, "_gmt/SAGI"))){
    file <- read.csv(paste0("~/work/abortvi2/rhapsodi/command_line_results/out/", gam_number, "_gmt/SAGI/", fileName), 
                     header=TRUE,
                     row.names=NULL,
                     stringsAsFactors = FALSE)#%>%
    #tibble::rownames_to_column("gamete")
    # For multi-row files, only get the last line. 
    if (nrow(file) != gam_number){
      file = tail(file, n = gam_number) 
    } 
    # Remove gamete indeces
    file <- subset( file, select = -row.names)
    file[] <- lapply(file, as.numeric)
    # Average columns for data columns and add to overall dataframe 
    gam_imputation_results <-  rbind(gam_imputation_results, c(file[1, 1:6], colMeans(file[7:10])))
  }
}

colnames(gam_imputation_results) = c("n_gametes","n_snp","cov","seqerr","avgr","rs","acc","com","lhs","ser")

# Add file name column for merging with rhapsodi data
gam_imputation_results <- gam_imputation_results %>%
  mutate(gam_imputation_results, fname = paste('g', n_gametes, 's', n_snp, 'c', cov, 'se', seqerr, 'r', avgr, 'rs', rs, sep="_"))

###########RHAPSODI#RESULTS############################################################

# Initialize dataframe for analysis output
rhapsodi_phasing_results <- data.frame(matrix(ncol=5, nrow=0))
colnames(rhapsodi_phasing_results) = c('fname', 'acc_rh', 'com_rh', 'lhs_rh', 'ser_rh')

rhapsodi_recombination_results <- data.frame(matrix(ncol=14, nrow=0))
colnames(rhapsodi_recombination_results) = c("fname","precision_rh","recall_rh","accuracy_rh","specificity_rh","fdr_rh","fpr_rh","f1_rh","true_n_rh","pred_n_rh","tn_rh","fn_rh","tp_rh","fp_rh")

rhapsodi_gam_imputation_results <- data.frame(matrix(ncol=5, nrow=0))
colnames(rhapsodi_gam_imputation_results) = c("fname","mean_acc_rh","mean_com_rh","mean_lhs_rh","mean_ser_rh")

## Loading Kate's rhapsodi assess data

# directory with new data
rhapsodi_dir <- "~/work/kweave23/assess_rhapsodi/changing_model_params/"

# loop through the numbers of gametes analyzed by both Hapi and Rhapsodi
for (n_gam in c(3, 15, 50, 150)){
  # extract all the files from the rhapsodi directory corresponding to this gamete number and loop through them 
  gam_pattern <- paste0("^g", n_gam, "_")
  for (dir in list.files(rhapsodi_dir, pattern = gam_pattern)){
    g = str_split(str_split(dir, '_')[[1]][1], 'g')[[1]][2]
    s = str_split(str_split(dir, '_')[[1]][2], 's')[[1]][2]
    c = str_split(str_split(dir, '_')[[1]][3], 'c')[[1]][2]
    se = str_split(str_split(dir, '_')[[1]][4], 'se')[[1]][2]
    r = str_split(str_split(dir, '_')[[1]][5], 'r')[[1]][2]
    
    # Loop through all random seeds
    for (rs in c(1848, 357, 42)){
      #Check if file exists
      input_fname = paste0("assess_out_rs_", rs, ".Rdata")
      
      # Check that rhapsodi output exists and if so load it
      if (input_fname %in% list.files(paste0(rhapsodi_dir,
                                             dir))){
        load(paste0(rhapsodi_dir,
                    dir, 
                    "/assess_out_rs_", rs, ".Rdata"))
        
        # get name of the simulation. used to merge with hapi output
        sim_fname <- paste('g', g, 's', s, 'c', c, 'se', se, 'r', r, 'rs', rs, sep="_")
        
        if ("phasing" %in% names(assess_out)){
          # Add data to phasing dataset
          rhapsodi_phasing_results[nrow(rhapsodi_phasing_results) + 1, ] <- c(sim_fname, assess_out$phasing)
        }
       
        if ("recomb" %in% names(assess_out)){
          # Add data to recombination dataset
          rhapsodi_recombination_results[nrow(rhapsodi_recombination_results) + 1, ] <- c(sim_fname, assess_out$recomb)
        }
       
        if ("gam_imputation" %in% names(assess_out)){
          # Manipulate gamete imputation data
          mean_acc <- mean(unlist(assess_out$gam_imputation$acc))
          mean_com <- mean(unlist(assess_out$gam_imputation$com))
          mean_lhs <- mean(unlist(assess_out$gam_imputation$lhs))
          mean_ser <- mean(unlist(assess_out$gam_imputation$ser))
          
          # Add data to gamete imputation dataset
          rhapsodi_gam_imputation_results[nrow(rhapsodi_gam_imputation_results) + 1, ] <- c(sim_fname, mean_acc, mean_com, mean_lhs, mean_ser)
        }
        
      }
    }
    
    
  }
  
}

#########MERGING#HAPI#AND#RHAPSODI#DATA###########################################################################################

phasing_merged <- merge(phasing_results, rhapsodi_phasing_results, by="fname") %>% 
  mutate(delta_acc = (acc_rh - acc)/100) %>%
  #mutate(delta_com = (com_rh - com) * 100) %>%
  mutate(delta_com = com_rh - com) %>%
  mutate(gmt_facet = paste(n_gametes, "gametes")) %>%
  mutate(gmt_facet = factor(gmt_facet, levels= c("3 gametes", "15 gametes", "50 gametes", "150 gametes")))

recombination_merged <- merge(recombination_results, rhapsodi_recombination_results, by="fname") %>%
  mutate(accuracy_rh = as.numeric(accuracy_rh)) %>%
  mutate(delta_accuracy = accuracy_rh - accuracy) %>%
  mutate(f1_rh = as.numeric(f1_rh)) %>%
  mutate(delta_f1 = f1_rh - f1) %>%
  mutate(specificity_rh = as.numeric(specificity_rh)) %>%
  mutate(delta_specificity = specificity_rh - specificity) %>%
  mutate(gmt_facet = paste(n_gametes, "gametes")) %>%
  mutate(gmt_facet = factor(gmt_facet, levels= c("3 gametes", "15 gametes", "50 gametes", "150 gametes")))

gamete_imputation_merged <- merge(gam_imputation_results, rhapsodi_gam_imputation_results, by="fname") %>% 
  mutate(mean_acc_rh = as.numeric(mean_acc_rh)) %>%
  mutate(delta_acc = (mean_acc_rh - acc)/100) %>%
  mutate(mean_com_rh = as.numeric(mean_com_rh)) %>%
  #mutate(delta_com = (mean_com_rh - com) * 100)%>%
  mutate(delta_com = (mean_com_rh - com)) %>%
  mutate(gmt_facet = paste(n_gametes, "gametes")) %>%
  mutate(gmt_facet = factor(gmt_facet, levels= c("3 gametes", "15 gametes", "50 gametes", "150 gametes")))

# Set color scale 
# Fall Foliage
cc <- scales::seq_gradient_pal("#9e0000", "#ffe53b", "Lab")(seq(0,1,length.out=9))
# Viridis
#cc <- scales::seq_gradient_pal("#fde725", "#440154", "Lab")(seq(0,1,length.out=9))
cc <- rev(c("#fde725", "#addc30", "#5ec962", "#28ae80", "#21918c", "#2c728e", "#3b528b", "#472d7b", "#440154"))

phasing_plot <- ggplot(data = phasing_merged, aes(x=delta_acc, y=delta_com, color = factor(cov))) + 
  geom_point(size = 0.8) + 
  theme_bw() + theme(panel.grid = element_blank()) +
  geom_hline(yintercept = 0, size = 0.2) + 
  geom_vline(xintercept = 0, size = 0.2) + 
  facet_grid(. ~ gmt_facet) + 
  xlim(-1,1) + 
  ylim(-1,1) + 
  xlab("\u0394 Accuracy") + 
  ylab("\u0394 Completeness") + 
  ggtitle("A: Donor Haplotype Phasing") +
  #scale_color_viridis(trans = "log10", guide = guide_legend(override.aes = list(alpha = 0) )) +
  #scale_colour_brewer(palette = "YlOrRd", direction=-1, guide = guide_legend(override.aes = list(alpha = 0) )) +
  scale_colour_manual(values=cc, guide = guide_legend(override.aes = list(alpha = 0) )) +
  theme(legend.title = element_text(color = "transparent"),
        legend.text = element_text(color = "transparent"),
        text = element_text(size = 10)) + 
  theme(axis.title.x=element_blank()) 

imputation_plot <- ggplot(data = gamete_imputation_merged, aes(x=delta_acc, y=delta_com, color = factor(cov))) + 
  geom_point(size = 0.8) + 
  geom_hline(yintercept = 0, size = 0.2) + 
  geom_vline(xintercept = 0, size = 0.2) + 
  #scale_colour_brewer(palette = "YlOrRd", direction=-1) +
  #scale_color_gradient(low="#b39f1b", high="#b31b1b") +
  scale_colour_manual(values=cc) +
  #scale_color_viridis(trans="log10") +
  facet_grid(. ~ gmt_facet) + 
  ggtitle("B: Gamete Genotype Imputation") +
  labs(color="Coverage") +
  xlim(-1,1) + 
  ylim(-1,1) + 
  xlab("\u0394 Accuracy") + 
  ylab("\u0394 Completeness") + 
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(
    axis.title.x=element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    text = element_text(size = 10)) 

#recombination_plot <- ggplot(data = recombination_merged, aes(x=delta_accuracy, y=delta_f1, color = factor(cov))) + 
recombination_plot <- ggplot(data = recombination_merged, aes(x=delta_accuracy, y=delta_specificity, color = factor(cov))) +
  geom_point(size = 0.8) + 
  theme_bw() + theme(panel.grid = element_blank()) +
  ggtitle("C: Meiotic Recombination Discovery") +
  geom_hline(yintercept = 0, size = 0.2) + 
  geom_vline(xintercept = 0, size = 0.2) + 
  facet_grid(. ~ gmt_facet) + 
  theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  text = element_text(size = 10)) +
  xlim(-1,1) + 
  ylim(-1,1) + 
  xlab("\u0394 Accuracy") + 
  #ylab("\u0394 F1") +
  ylab("\u0394 Specificity") +
  #scale_color_viridis(trans="log10", guide = guide_legend(override.aes = list(alpha = 0) )) +
  #scale_colour_brewer(palette = "YlOrRd", direction=-1, guide = guide_legend(override.aes = list(alpha = 0) )) +
  scale_colour_manual(values=cc, guide = guide_legend(override.aes = list(alpha = 0) )) +
  theme(legend.title = element_text(color = "transparent"),
        legend.text = element_text(color = "transparent")) +
  theme(panel.grid = element_blank()) 

p1 <- phasing_plot / imputation_plot / recombination_plot

p1

####Plotting just Rhapsodi Results#########################################################################

# split file name into files, filter by gamete number
rhapsodi_phasing_results_plt <- rhapsodi_phasing_results %>% 
  separate(fname, c("first", "g", "second", "s", "third", "cov"), sep = "_", remove = FALSE) %>% select(-c('first', 'second', 'third')) %>%
  filter(as.numeric(g) < 151)

# Find files that are in the rhapsodi results and not in hapi results
sdiff_phasing <- setdiff(sort(rhapsodi_phasing_results$fname), sort(phasing_results$fname))

# split file name into files, filter by gamete number
rhapsodi_recombination_results_plt <- rhapsodi_recombination_results %>% 
  separate(fname, c("first", "g", "second", "s", "third", "cov"), sep = "_", remove = FALSE) %>% select(-c('first', 'second', 'third')) %>%
  filter(as.numeric(g) < 151)

# Find files that are in the rhapsodi results and not in hapi results
sdiff_recomb <- setdiff(sort(rhapsodi_recombination_results$fname), sort(recombination_results$fname))

# split file name into files, filter by gamete number
rhapsodi_gam_imputation_results_plt <- rhapsodi_gam_imputation_results %>% 
  separate(fname, c("first", "g", "second", "s", "third", "cov"), sep = "_", remove = FALSE) %>% select(-c('first', 'second', 'third')) %>%
  filter(as.numeric(g) < 151)

# Find files that are in the rhapsodi results and not in hapi results
sdiff_imputation <- setdiff(sort(rhapsodi_gam_imputation_results$fname), sort(gam_imputation_results$fname))

## Add a column for each input file in the rhapsodi dataframes indicating whether hapi was able to run that file correctly
rhapsodi_phasing_results_plt <- rhapsodi_phasing_results_plt %>%
  mutate(hapi_completed = ifelse(fname %in% sdiff_phasing, 0, 1)) %>%
  #mutate(com_rh = 100*com_rh) %>%
  mutate(acc_rh = as.numeric(acc_rh)/100) %>%
  mutate(gmt_facet = paste(g, "gametes")) %>%
  mutate(gmt_facet = factor(gmt_facet, levels= c("3 gametes", "15 gametes", "50 gametes", "150 gametes")))

rhapsodi_recombination_results_plt <- rhapsodi_recombination_results_plt %>%
  mutate(hapi_completed = ifelse(fname %in% sdiff_recomb, 0, 1)) %>%
  mutate(accuracy_rh = accuracy_rh) %>%
  mutate(gmt_facet = paste(g, "gametes")) %>%
  mutate(gmt_facet = factor(gmt_facet, levels= c("3 gametes", "15 gametes", "50 gametes", "150 gametes")))

rhapsodi_gam_imputation_results_plt <- rhapsodi_gam_imputation_results_plt %>%
  mutate(hapi_completed = ifelse(fname %in% sdiff_imputation, 0, 1))%>%
  #mutate(mean_com_rh = 100*as.numeric(mean_com_rh)) %>%
  mutate(mean_com_rh = as.numeric(mean_com_rh)) %>%
  mutate(mean_acc_rh = as.numeric(mean_acc_rh)/100) %>%
  mutate(gmt_facet = paste(g, "gametes")) %>%
  mutate(gmt_facet = factor(gmt_facet, levels= c("3 gametes", "15 gametes", "50 gametes", "150 gametes")))

#### Plotting

# Set color scale 
#cc <- scales::seq_gradient_pal("#9e0000", "#ffe53b", "Lab")(seq(0,1,length.out=2))
cc <- c("#28ae80", "#440154")

rhapsodi_phasing_plot <- ggplot(data = rhapsodi_phasing_results_plt, aes(x=acc_rh, y=com_rh, color = factor(hapi_completed))) + 
  geom_point(size = 0.8) + 
  theme_minimal() + 
  geom_hline(yintercept = 0, size = 0.2) + 
  geom_vline(xintercept = 0, size = 0.2) + 
  facet_grid(. ~ gmt_facet) + 
  xlim(0,1) + 
  ylim(0,1) + 
  ggtitle("A: Donor Haplotype Phasing") +
  xlab("Accuracy") + 
  ylab("Completeness") + 
  theme_bw() + theme(panel.grid = element_blank()) +
  #scale_colour_brewer(palette = "YlOrRd", direction=-1, guide = guide_legend(override.aes = list(alpha = 0) )) +
  theme(legend.title = element_text(color = "transparent"),
        legend.text = element_text(color = "transparent"), 
        text = element_text(size = 10)) + 
  scale_colour_manual(values=cc, guide = guide_legend(override.aes = list(alpha = 0) )) +
  theme(axis.title.x=element_blank())

rhapsodi_imputation_plot <- ggplot(data = rhapsodi_gam_imputation_results_plt, aes(x=as.numeric(mean_acc_rh), y=as.numeric(mean_com_rh), color = factor(hapi_completed))) + 
  geom_point(size = 0.8) + 
  theme_minimal() + 
  geom_hline(yintercept = 0, size = 0.2) + 
  geom_vline(xintercept = 0, size = 0.2) + 
  #scale_colour_brewer(palette = "YlOrRd", direction=-1) +
  #scale_color_gradient(low="#b39f1b", high="#b31b1b") +
  scale_colour_manual(values=cc, labels = c("False","True")) + 
  facet_grid(. ~ gmt_facet) + 
  theme_bw() + theme(panel.grid = element_blank()) +
  ggtitle("B: Gamete Genotype Imputation") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    text = element_text(size = 10)) + 
  labs(color="Hapi Completion") +
  xlim(0,1) + 
  ylim(0,1) + 
  xlab("Accuracy") + 
  ylab("Completeness") + 
  theme(axis.title.x=element_blank())

#rhapsodi_recombination_plot <- ggplot(data = rhapsodi_recombination_results_plt, aes(x=accuracy_rh, y=f1_rh, color = factor(hapi_completed))) + 
rhapsodi_recombination_plot <- ggplot(data = rhapsodi_recombination_results_plt, aes(x=accuracy_rh, y=specificity_rh, color = factor(hapi_completed))) + 
  geom_point(size = 0.8) + 
  #theme_minimal() + 
  geom_hline(yintercept = 0, size = 0.2) + 
  geom_vline(xintercept = 0, size = 0.2) + 
  theme_bw() + theme(panel.grid = element_blank()) +
  facet_grid(. ~ gmt_facet) + 
  ggtitle("C: Meiotic Recombination Discovery") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    text = element_text(size = 10) 
    #text = element_text(size = 14, face = "bold")
    ) + 
  xlim(0,1) + 
  ylim(0,1) + 
  xlab("Accuracy") + 
  #ylab("F1") +
  ylab("Specificity") +
  #scale_colour_brewer(palette = "YlOrRd", direction=-1, guide = guide_legend(override.aes = list(alpha = 0) )) +
  scale_colour_manual(values=cc, guide = guide_legend(override.aes = list(alpha = 0) )) +
  theme(legend.title = element_text(color = "transparent"),
        legend.text = element_text(color = "transparent")) 

p2 <- rhapsodi_phasing_plot/rhapsodi_imputation_plot/rhapsodi_recombination_plot

p2

##########SURVIVAL###CURVE########################################################################

# Initialize dataframe for log file read-in
hapi_logs <- data.frame(matrix(ncol=15, nrow=0))
colnames(hapi_logs) = c("n_gam","n_snp","n_het_snp","cov","seqerr","avgr","rs","gmt_import","filter_error","frame_select","impute_and_phase","final_draft","consenus","consensus_end","CO")

# Loop through gamete numbers
for (gam_number in c(3, 15, 50, 150)){
  # Loop through all files starting with the string "log"
  for (log_file in list.files(paste0("/home-1/abortvi2@jhu.edu/work/abortvi2/rhapsodi/command_line_results/logs/", gam_number, "_gmt/"),
                              pattern = "^log")){
    # Read in file
    file <- read.csv(paste0("/home-1/abortvi2@jhu.edu/work/abortvi2/rhapsodi/command_line_results/logs/", gam_number, "_gmt/", log_file),
                     header=TRUE,
                     stringsAsFactors = FALSE)
    # Get out final row
    if (nrow(file) != 1){
      file = file[nrow(file),]
    } 
    # Add to dataframe of log files
    hapi_logs[nrow(hapi_logs) + 1, ] <-  file
  }
}

## Plotting bulk 
hapi_survival_df <- data.frame(matrix(ncol=15, nrow=0))
colnames(hapi_survival_df) = c('stage', 'hapi_survival', 'percent_survival', 'n_gam')

#loop through gamete number 
for (n_gam in unique(hapi_logs$n_gam)){
  hapi_logs_subset <- hapi_logs[hapi_logs$n_gam == n_gam, c('gmt_import', 'impute_and_phase', 'consensus_end', 'CO')]
  colnames(hapi_logs_subset) <- c('gmt_import', 'phasing', 'imputation', 'CO')
  # Get the total number of simulations 
  num_hapi_sims <- nrow(hapi_logs_subset)
  # Sum columns
  hapi_survival <- colSums(hapi_logs_subset, na.rm = TRUE)
  # Manipulate data for plotting
  hapi_survival <- as.data.frame(hapi_survival) %>%
    tibble::rownames_to_column("stage") %>%
    mutate(stage = factor(stage, levels= c('gmt_import', 'phasing', 'imputation', 'CO'))) %>%
    mutate(percent_survival = hapi_survival/num_hapi_sims)
  to_bind <- cbind(hapi_survival, n_gam)
  hapi_survival_df <- rbind(hapi_survival_df, to_bind)
}

survival_plot <- ggplot(data = hapi_survival_df, aes(x = stage, y = percent_survival, group = 1)) + 
  geom_step() +
  facet_wrap('~n_gam', ncol = 1)+
  theme_classic() + 
  ylim(0,1)

survival_plot

###########RHAPSODI#SURVIVAL############################################################

# Initialize dataframe for analysis output
# rhapsodi_phasing_results <- data.frame(matrix(ncol=5, nrow=0))
# colnames(rhapsodi_phasing_results) = c('fname', 'acc_rh', 'com_rh', 'lhs_rh', 'ser_rh')
# 
# rhapsodi_recombination_results <- data.frame(matrix(ncol=14, nrow=0))
# colnames(rhapsodi_recombination_results) = c("fname","precision_rh","recall_rh","accuracy_rh","specificity_rh","fdr_rh","fpr_rh","f1_rh","true_n_rh","pred_n_rh","tn_rh","fn_rh","tp_rh","fp_rh")
# 
# rhapsodi_gam_imputation_results <- data.frame(matrix(ncol=5, nrow=0))
# colnames(rhapsodi_gam_imputation_results) = c("fname","mean_acc_rh","mean_com_rh","mean_lhs_rh","mean_ser_rh")

## Loading Kate's rhapsodi assess data

# directory with new data
rhapsodi_dir <- "~/work/kweave23/assess_rhapsodi/changing_model_params/"

# Initialize empty df 
rhapsodi_survival_df <- data.frame(matrix(ncol=15, nrow=0))
colnames(rhapsodi_survival_df) = c("gmt_import", "phasing", "imputation", "CO", "gam")

# loop through the numbers of gametes analyzed by both Hapi and Rhapsodi
for (n_gam in c(3, 15, 50, 150)){
  # extract all the files from the rhapsodi directory corresponding to this gamete number and loop through them 
  gam_pattern <- paste0("^g", n_gam, "_")
  for (dir in list.files(rhapsodi_dir, pattern = gam_pattern)){
    g = str_split(str_split(dir, '_')[[1]][1], 'g')[[1]][2]
    s = str_split(str_split(dir, '_')[[1]][2], 's')[[1]][2]
    c = str_split(str_split(dir, '_')[[1]][3], 'c')[[1]][2]
    se = str_split(str_split(dir, '_')[[1]][4], 'se')[[1]][2]
    r = str_split(str_split(dir, '_')[[1]][5], 'r')[[1]][2]
    
    # Loop through all random seeds
    for (rs in c(1848, 357, 42)){
      #Check if file exists
      input_fname = paste0("survivor_arr_rs_", rs, ".csv")
      
      # Check that rhapsodi output exists and if so load it
      if (input_fname %in% list.files(paste0(rhapsodi_dir,dir))){
         temp_survival_data <- read.csv((paste0(rhapsodi_dir,
                     dir, '/',
                     input_fname)))
         
         rownames(temp_survival_data) <- c('gmt_import', 'phasing', 'imputation', 'CO', 'assess_phasing', 'assess_imputation', 'assess_CO')
         colnames(temp_survival_data) <- paste(g, s, c, se, r, sep='_')
         temp_survival_data <- rbind(temp_survival_data, "gam" = g)
         temp_survival_data_T <- as.data.frame(t(as.matrix(temp_survival_data))) %>% 
           select(c("gmt_import", "phasing", "imputation", "CO", "gam"))
         
         rhapsodi_survival_df <- rbind(rhapsodi_survival_df, temp_survival_data_T)
         }
        
      }
    }
}
  
# Summed df of rhaspodi survival for plotting
summed_rhapsodi_survival <- data.frame(matrix(ncol=4, nrow=0))
colnames(summed_rhapsodi_survival) <- c('stage', 'rhapsodi_survival', 'percent_survival', 'n_gam')

for (n_gam in c(3, 15, 50, 150)){
  rh_logs_subset <- rhapsodi_survival_df[rhapsodi_survival_df$gam == n_gam, !(names(rhapsodi_survival_df) == "gam")]
  num_rhapsodi_sims <- nrow(rh_logs_subset)
  rh_logs_subset <- rh_logs_subset %>% 
    mutate(gmt_import = as.numeric(as.character(gmt_import))) %>%
    mutate(phasing = as.numeric(as.character(phasing))) %>%
    mutate(imputation = as.numeric(as.character(imputation))) %>%
    mutate(CO = as.numeric(as.character(CO)))
  temp_rhapsodi_survival <- colSums(rh_logs_subset, na.rm = TRUE)
  
  temp_rhapsodi_survival <- as.data.frame(temp_rhapsodi_survival) %>%
    tibble::rownames_to_column("stage") %>%
    mutate(stage = factor(stage, levels= c('gmt_import', 'phasing', 'imputation', 'CO'))) %>%
    mutate(temp_rhapsodi_survival = as.numeric(as.character(temp_rhapsodi_survival))) %>%
    mutate(percent_survival = temp_rhapsodi_survival/num_rhapsodi_sims)
  
  temp_rhapsodi_survival <- cbind(temp_rhapsodi_survival, "n_gam" = n_gam) %>%
    mutate(n_gam = as.numeric(as.character(n_gam)))
  
  colnames(temp_rhapsodi_survival) <- c('stage', 'rhapsodi_survival', 'percent_survival', 'n_gam')
  
  summed_rhapsodi_survival <- rbind(summed_rhapsodi_survival, temp_rhapsodi_survival)
}

summed_rhapsodi_survival <- cbind(summed_rhapsodi_survival, 'tool' = 'rhapsodi')
colnames(summed_rhapsodi_survival) <- c("stage", "number_survival", "percent_survival", "n_gam", "Tool")
hapi_survival_df <- cbind(hapi_survival_df, 'tool' = 'hapi')
colnames(hapi_survival_df) <- c("stage", "number_survival", "percent_survival", "n_gam", "Tool")

all_survival <- rbind(summed_rhapsodi_survival, hapi_survival_df)

all_survival_plot <- ggplot(data = all_survival, aes(x = stage, y = percent_survival, group = interaction(n_gam, Tool), linetype = Tool, color = as.factor(n_gam))) + 
  #facet_wrap('~n_gam', ncol = 1)+
  theme_bw() + theme(panel.grid = element_blank()) +
  geom_step(size = 1) +
  xlab('Stage') +
  ylab('Percent Completed') +
  ylim(0,1) + guides(color=guide_legend(title="Number of Gametes", order = 1),
                     linetype=guide_legend(order = 4))

all_survival_plot 






















































































































