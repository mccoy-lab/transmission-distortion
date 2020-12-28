#!/bin/bash

#SBATCH --job-name=recomb_y_pval
#SBATCH -N 1
#SBATCH --partition=shared
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --time=72:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=kweave23@jhu.edu
#SBATCH --array=1-22%22
#SBATCH --account=rmccoy22

module load R/4.0.2
module load gcc/5.5.0

echo JobID_${SLURM_JOB_ID}_Sample_${SAMPLE_DIR}_chr${SLURM_ARRAY_TASK_ID} > ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_20201228/out_pval_y_rs_kw_20201228_${SAMPLE_DIR}_${SLURM_ARRAY_TASK_ID}.out
Rscript assign_sperm_haplotypes_rm_kw.R ~/work-rmccoy22/rmccoy22/mccarroll-sperm-seq-data/${SAMPLE_DIR}/data/${SAMPLE_DIR}_goodcellsreplicatebcs_filteredhetsnps_${SLURM_ARRAY_TASK_ID}.cellsbyrow.txt ${SAMPLE_DIR} chr${SLURM_ARRAY_TASK_ID} ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_20201228/${SAMPLE_DIR}/csv_out/ 0.005 2L 2500 &>> ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_20201228/out_pval_y_rs_kw_20201228_${SAMPLE_DIR}_${SLURM_ARRAY_TASK_ID}.out
