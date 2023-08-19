# A method for investigating transmission distortion among human sperm

This repository contains the analyses performed for the study "A method for low-coverage single-gamete sequence analysis demonstrates adherence to Mendel’s first law across a large sample of human sperm". This paper has been publised in [eLife](https://elifesciences.org/articles/76383)

<!--
- filtering_bell_data: info for filtering our Bell data based on genome consortia studies
- full_donors: data from Bell et al. 2020
- plotting: scripts for figures  
- shell-scripts: command line scripts for processing raw data  
- sim-scripts: R files for simulations of TD (real and simulated chromosomes)
- sperm-data: steps for processing genotype data
-->

# Results

## Evaluating performance on simulated data

### Data generation

For each different combination of study design:

* [`sim-scripts/generative_model_for_rhapsodi.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/generative_model_for_rhapsodi.R) run by [`sim-scripts/run_genModel.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/run_genModel.sh) or [`sim-scripts/run_genModel_3args.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/run_genModel_3args.sh)

Followed by one of the following:

* [`sim-scripts/assess_with_rhapsodi.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/assess_with_rhapsodi.R) run by [`sim-scripts/run_assess_with_rhapsodi.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/run_assess_with_rhapsodi.sh) or [`sim-scripts/run_assess_with_rhapsodi_2args.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/run_assess_with_rhapsodi_2args.sh)
* [`sim-scripts/assess_with_rhapsodi_thread_arg.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/assess_with_rhapsodi_thread_arg.R) run by [`sim-scripts/run_assess_with_rhapsodi_thread_arg.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/run_assess_with_rhapsodi_thread_arg.sh)
* [`sim-scripts/assess_with_rhapsodi_mparams.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/assess_with_rhapsodi_mparams.R) run by [`sim-scripts/run_assess_with_rhapsodi_mparams.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/run_assess_with_rhapsodi_mparams.sh)

### Statstics reported in text

* [`analysis_scripts/get_sim_stats.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/analysis_scripts/get_sim_stats.R)

### Plotting

#### Fig 2

* [`plotting/plot_rhapsodi_main_gen_fig.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/plot_rhapsodi_main_gen_fig.R)
  * Rdata files for streamlined plotting in [`plotting/main_heatmap_rdata`](https://github.com/mccoy-lab/transmission-distortion/tree/master/plotting/main_heatmap_Rdata)

#### Fig 2 Supplemental Fig 2

* [`plotting/plot_rhapsodi_supp_gen_fig.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/plot_rhapsodi_supp_gen_fig.R)
  * Rdata files for streamlined plotting in [`plotting/supp_heatmap_rdata`](https://github.com/mccoy-lab/transmission-distortion/tree/master/plotting/supp_heatmap_Rdata)

#### Fig 2 Supplemental Fig 3

Make Rdata files 
* [`sim-scripts/save_true_nsnp_vals.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/save_true_nsnp_vals.R)
  * Rdata file in [`test_data_rhapsodi_gen`](https://github.com/mccoy-lab/transmission-distortion/tree/master/test_data_rhapsodi_gen)
* [`plotting/supp_recomb/plot_supfig3_recombination_info.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/supp_recomb/plot_supfig3_recombination_info.R)
  * Rdata files [`plotting/supfig3_rdata`](https://github.com/mccoy-lab/transmission-distortion/tree/master/plotting/supfig3_rdata) 

Plot
* [`plotting/supp_recomb/plot_bp_res.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/supp_recomb/plot_bp_res.R) and [`plotting/supp_recomb/plot_supfig3_recomb_fn_fp.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/supp_recomb/plot_supfig3_recomb_fn_fp.R)
  * Rdata files for streamlined plotting in [`plotting/supp_recomb`](https://github.com/mccoy-lab/transmission-distortion/tree/master/plotting/supp_recomb)

#### Fig 2 Supplemental Fig 4-7

* [`plotting/robustness/plot_rhapsodi_robust_fig.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/robustness/plot_rhapsodi_robust_fig.R)
  * Rdata files for streamlined plotting in [`plotting/robustness`](https://github.com/mccoy-lab/transmission-distortion/tree/master/plotting/robustness)


## Benchmarking against existing methods: Hapi and HapCUT2

### Data generation

* Same from data generation in "Evaluating performance on simulated data"
* [`analysis_scripts/hapi_benchmarking.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/analysis_scripts/hapi_benchmarking.R)

### Plotting

#### Fig 3, 3S1, 3S2

* [`plotting/rhapsodi_hapi_benchmarking_plot_scripts.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/rhapsodi_hapi_benchmarking_plot_scripts.R)


## Application to data from human sperm

### Data generation

* [`sperm-data/subset_data.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sperm-data/subset_data.R) run by [`sperm-data/subset_bash.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sperm-data/subset_bash.sh)
* [`filter_rhapsodi_TDscan.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/filter_rhapsodi_TDscan.R) run by [`slurm_filter_rhapsodi_TDscan.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/slurm_filter_rhapsodi_TDscan.sh) and [`submit_slurm_filter_rhapsodi_TDscan.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/submit_slurm_filter_rhapsodi_TDscan.sh)

### Statistics reported in text

* [`analysis_scripts/get_stats.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/analysis_scripts/get_stats.R)


## Statistical power to detect moderate and strong TD

### Plotting 

#### Fig 4, 4S1

* [`plotting/power_analyses/power_analysis_sperm.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/TD_sim_plotting.R)


#### Fig 4S2, 4S3

* [`plotting/TD_sim_plotting.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/TD_sim_plotting.R)


#### Fig 4S4

* [`plotting/power_analyses/informative_transmissions_TDT.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/TD_sim_plotting.R)



## Strict adherence to Mendelian expectations across sperm genomes

### Data generation

* Same from data generation in "Application to data from human sperm"

### Plotting

#### Fig 5

* [`plotting/manhattan_all_donors_11052021.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/manhattan_all_donors_11052021.R)

#### Fig 5S2

* ??

#### Fig 5S3

* [`plotting/identify_segmental_dups.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/identify_segmental_dups.R)

#### Fig 5S4, 5S5

* [`plotting/plot_recombination_map.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/plot_recombination_map.R)
* Rdata files for streamlined plotting in [`plotting/supp_recomb`](https://github.com/mccoy-lab/transmission-distortion/tree/master/plotting/supp_recomb)



## No global signal of biased transmission in human sperm

### Data generation

* [`run_plink/make_ped_file.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/run_plink/make_ped_file.R) run by [`run_plink/submit_slurm_make_ped.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/run_plink/submit_slurm_make_ped.sh) calling [`run_plink/slurm_make_ped.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/run_plink/slurm_make_ped.sh) and [`run_plink/slurm_inf_make_ped.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/run_plink/slurm_inf_make_ped.sh)
<!--* [`plink/run_plink.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plink/run_plink.sh) OR [`run_plink/run_plink.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/run_plink/run_plink.sh)-->
* [`analysis_scripts/get_snps_prune.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/analysis_scripts/get_snps_prune.R)
* [`sim-scripts/null_sim.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/null_sim.R) run by [`sim-scripts/sim_wrapper.sh`](https://github.com/mccoy-lab/transmission-distortion/blob/master/sim-scripts/sim_wrapper.sh)

### Plotting

#### Fig 6
* Panel A: [`plotting/qqplot.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/qqplot.R) 
* Panel B: [`plotting/plot_global_sims.R`](https://github.com/mccoy-lab/transmission-distortion/blob/master/plotting/plot_global_sims.R)
