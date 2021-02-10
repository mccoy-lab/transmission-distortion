#!/bin/bash

#SBATCH --job-name=no_smooth
#SBATCH -N 1
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=72:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=kweave23@jhu.edu
#SBATCH --array=1-22%22
#SBATCH --account=rmccoy22

module load R/4.0.2
module load gcc/5.5.0

echo JobID_${SLURM_JOB_ID}_Sample_${SAMPLE_DIR}_chr${SLURM_ARRAY_TASK_ID} > ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_wo_smooth_20210210/${SAMPLE_DIR}/out_kw_20210210_${SAMPLE_DIR}_${SLURM_ARRAY_TASK_ID}.out
Rscript assign_sperm_haplotypes_wo_smooth.R ~/work-rmccoy22/scarios1/bell_filtered_data/${SAMPLE_DIR}/${SAMPLE_DIR}_euploid_${SLURM_ARRAY_TASK_ID}.txt ${SAMPLE_DIR} chr${SLURM_ARRAY_TASK_ID} ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_wo_smooth_20210210/${SAMPLE_DIR}/csv_out/ 0.005 8 3000 FALSE &>> ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_wo_smooth_20210210/${SAMPLE_DIR}/out_kw_20210210_${SAMPLE_DIR}_${SLURM_ARRAY_TASK_ID}.out
