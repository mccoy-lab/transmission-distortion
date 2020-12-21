#!/bin/bash

#SBATCH --job-name=sperm_HMM
#SBATCH -N 1
#SBATCH --partition=shared
#SBATCH ---mem=10GB
#SBATCH --time=72:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=kweave23@jhu.edu
#SBATCH --array=1-22%22
#SBATCH --account=rmccoy22

source activate R_env

Rscript assign_sperm_haplotypes_kw.R ~/work-rmccoy22/rmccoy22/mccarroll-sperm-seq-data/${SAMPLE_DIR}/data/${SAMPLE_DIR}_goodcellsreplicatebcs_filteredhetsnps_${SLURM_ARRAY_TASK_ID}.cellsbyrow.txt ${SAMPLE_DIR}_chr${SLURM_ARRAY_TASK_ID} 0.005 &> out_assign_sperm_haplotypes_kw_20201220_${SAMPLE_DIR}_${SLURM_ARRAY_TASK_ID}.out

conda deactivate
