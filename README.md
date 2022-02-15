# A method for investigating transmission distortion among human sperm

This repository includes the analyses performed for the study "Strict adherence to Mendel's First Law across a large sample of human sperm genomes". This study is currently [posted on bioRxiv](https://www.biorxiv.org/content/10.1101/2021.11.19.469261v2)

<!--
- filtering_bell_data: info for filtering our Bell data based on genome consortia studies
- full_donors: data from Bell et al. 2020
- plotting: scripts for figures  
- shell-scripts: command line scripts for processing raw data  
- sim-scripts: R files for simulations of TD (real and simulated chromosomes)
- sperm-data: steps for processing genotype data
-->

# Results: Evaluating performance on simulated data

## Data Generation

For each different combo of study designs:

* [`sim-scripts/generative_model_for_rhapsodi.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/generative_model_for_rhapsodi.R) run by [`sim-scripts/run_genModel.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/run_genModel.sh) or [`sim-scripts/run_genModel_3args.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/run_genModel_3args.sh)

followed by one of the following

* [`sim-scripts/assess_with_rhapsodi.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/assess_with_rhapsodi.R) run by [`sim-scripts/run_assess_with_rhapsodi.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/run_assess_with_rhapsodi.sh) or [`sim-scripts/run_assess_with_rhapsodi_2args.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/run_assess_with_rhapsodi_2args.sh)
* [`sim-scripts/assess_with_rhapsodi_thread_arg.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/assess_with_rhapsodi_thread_arg.R) run by [`sim-scripts/run_assess_with_rhapsodi_thread_arg.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/run_assess_with_rhapsodi_thread_arg.sh)
* [`sim-scripts/assess_with_rhapsodi_mparams.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/assess_with_rhapsodi_mparams.R) run by [`sim-scripts/run_assess_with_rhapsodi_mparams.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/run_assess_with_rhapsodi_mparams.sh)

## Stats reported in text

* [`get_sim_stats.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/get_sim_stats.R)

## Plotting

### Fig 2

* [`plotting/plot_rhapsodi_main_gen_fig.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/plot_rhapsodi_main_gen_fig.R)
  * Rdata files for streamlined plotting in [`plotting/main_heatmap_rdata`](https://github.com/mccoy-lab/transmission-distortion/tree/master/plotting/main_heatmap_Rdata)

### Fig 2 Supplemental Fig 2

* [`plotting/plot_rhapsodi_supp_gen_fig.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/plot_rhapsodi_supp_gen_fig.R)
  * Rdata files for streamlined plotting in [`plotting/supp_heatmap_rdata`](https://github.com/mccoy-lab/transmission-distortion/tree/master/plotting/supp_heatmap_Rdata)

### Fig 2 Supplemental Fig 3

Making of Rdata files from
* [`sim-scripts/save_true_nsnp_vals.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/save_true_nsnp_vals.R)
  * Rdata file in [`test_data_rhapsodi_gen`](https://github.com/mccoy-lab/transmission-distortion/tree/master/test_data_rhapsodi_gen)
* [`plotting/supp_recomb/plot_supfig3_recombination_info.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/supp_recomb/plot_supfig3_recombination_info.R)
  * Rdata files saved in [`plotting/supfig3_rdata`](https://github.com/mccoy-lab/transmission-distortion/tree/master/plotting/supfig3_rdata) and later used by actual plotting scripts

Plotting itself
* [`plotting/supp_recomb/plot_bp_res.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/supp_recomb/plot_bp_res.R) and [`plotting/supp_recomb/plot_supfig3_recomb_fn_fp.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/supp_recomb/plot_supfig3_recomb_fn_fp.R)
  * Rdata files for streamlined plotting in [`plotting/supp_recomb`](https://github.com/mccoy-lab/transmission-distortion/tree/master/plotting/supp_recomb)

### Fig 2 Supplemental Fig 4 - 7

* [`plotting/robustness/plot_rhapsodi_robust_fig.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/robustness/plot_rhapsodi_robust_fig.R)
  * Rdata files for streamlined plotting in [`plotting/robustness`](https://github.com/mccoy-lab/transmission-distortion/tree/master/plotting/robustness)

# Results: Benchmarking against existing methods

## Data Generation

* Same from data generation in "Evaluating performance on simulated data"
* [`hapi_benchmarking.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/hapi_benchmarking.R)

## Plotting

### Fig 3, Fig 3 Supplemental Fig 1, & Fig 3 Supplemental Fig 2

* [`plotting/rhapsodi_hapi_benchmarking_plot_scripts.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/rhapsodi_hapi_benchmarking_plot_scripts.R)

# Results: Application to data from human sperm

## Data Generation

* [`sperm-data/subset_data.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sperm-data/subset_data.R) run by [`sperm-data/subset_bash.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sperm-data/subset_bash.sh)
* [`filter_rhapsodi_TDscan.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/filter_rhapsodi_TDscan.R) run by [`slurm_filter_rhapsodi_TDscan.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/slurm_filter_rhapsodi_TDscan.sh) and [`submit_slurm_filter_rhapsodi_TDscan.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/submit_slurm_filter_rhapsodi_TDscan.sh)

## Stats reported in text

* [`get_stats.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/get_stats.R)

# Results: Strict adherence to Mendelian expectations across sperm genomes

## Data Generation

* Same as "Application to data from human sperm section"

## Plotting

### Fig 4

* [`plotting/manhattan_all_donors_11052021.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/manhattan_all_donors_11052021.R)

### Fig 4 Supplemental Fig 2

* [`plotting/TD_sim_plotting.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/TD_sim_plotting.R)

### Fig 4 Supplemental Fig 3

* [`plotting/plot_recombination_map.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/plot_recombination_map.R)
  * Rdata files for streamlined plotting in [`plotting/supp_recomb`](https://github.com/mccoy-lab/transmission-distortion/tree/master/plotting/supp_recomb)

### Fig 4 Supplemental Fig 4

* [`plotting/power_analysis.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/power_analysis.R)

### Fig 4 Supplemental Fig 5

* [`plotting/identify_segmental_dups.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/identify_segmental_dups.R)

# Results: No global signal of biased transmission in human sperm

## Data Generation

* [`run_plink/make_ped_file.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/run_plink/make_ped_file.R) run by [`run_plink/submit_slurm_make_ped.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/run_plink/submit_slurm_make_ped.sh) calling [`run_plink/slurm_make_ped.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/run_plink/slurm_make_ped.sh) and [`run_plink/slurm_inf_make_ped.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/run_plink/slurm_inf_make_ped.sh)
<!--* [`plink/run_plink.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plink/run_plink.sh) OR [`run_plink/run_plink.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/run_plink/run_plink.sh)-->
* [`get_snps_prune.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/get_snps_prune.R)
* [`sim-scripts/null_sim.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/null_sim.R) run by [`sim-scripts/sim_wrapper.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/sim_wrapper.sh)

## Plotting

### Fig 5
<!--* Panel A: [`plotting/TD_plot_qq_tile.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/TD_plot_qq_tile.R) OR [`plotting/qqplot.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/qqplot.R) NOTE both use the two original and infertile separate runs; I don't see an obvious difference between the two-->
* Panel B: [`plotting/plot_global_sims.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/plot_global_sims.R)
