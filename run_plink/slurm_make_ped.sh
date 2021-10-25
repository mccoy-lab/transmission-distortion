#!/bin/bash

#SBATCH --job-name=make_ped
#SBATCH -N 1
#SBATCH --partition=shared
#SBATCH --mem=100G
#SBATCH --time=72:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=kweave23@jhu.edu
#SBATCH --array=1-22%22

module load R/4.0.2
module load gcc/5.5.0

echo JobID_${SLURM_JOB_ID}_Sample_${SAMPLE_DIR}_chr${SLURM_ARRAY_TASK_ID} > ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_make_ped/${SAMPLE_DIR}/out_make_ped_${SAMPLE_DIR}_${SLURM_ARRAY_TASK_ID}.out
date; time Rscript make_ped_file.R ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_full_filter_20210826/${SAMPLE_DIR}/csv_out/${SAMPLE_DIR}_${SLURM_ARRAY_TASK_ID}_filled_sperm_unsmoothed.csv ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_full_filter_20210826/${SAMPLE_DIR}/csv_out/${SAMPLE_DIR}_${SLURM_ARRAY_TASK_ID}_parental_hap.csv ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_make_ped/${SAMPLE_DIR}/${SAMPLE_DIR}_${SLURM_ARRAY_TASK_ID}_full_sperm_converted.ped ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_make_ped/${SAMPLE_DIR}/${SAMPLE_DIR}_${SLURM_ARRAY_TASK_ID}_full_sperm_converted.map $SLURM_ARRAY_TASK_ID &>> ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_make_ped/${SAMPLE_DIR}/out_make_ped_${SAMPLE_DIR}_${SLURM_ARRAY_TASK_ID}.out
