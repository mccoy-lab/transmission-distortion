#!/bin/bash

ml R/4.0.2
ml gcc/5.5.0

echo "Plot Main Heatmap"
date; time Rscript plot_rhapsodi_main_gen_fig.R
echo "Plotted Main Heatmap and saved to plots/ directory"

echo "Plot Supplementary Heatmap"
date; time Rscript plot_rhapsodi_supp_gen_fig.R
echo "Plotted Supplementary Heatmap and saved to plots/ directory"

echo "Plot Robustness Heatmaps"
date; time Rscript plot_rhapsodi_robust_fig.R
echo "Plotted Robustness Heatmaps and saved to plots/ directory"

echo "Plot Recomb Resolution Boxplot"
date; time Rscript plot_bp_res.R
echo "Plotted Recomb Resolution Boxplot and saved to plots/ directory"

echo "Plot Recomb FP and FN Scatterplot"
date; time Rscript plot_supfig3_recomb_fn_fp.R
echo "Plotted Recomb FP and FN Scatterplot and saved to plots/ directory"

echo "Plot Recombination Map"
date; time Rscript plot_recombination_map.R
echo 'Plotted Recombination Map Barplot and saved to plots/ directory"
