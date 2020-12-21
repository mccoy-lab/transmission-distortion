#!/bin/bash

#SBATCH --job-name=sperm_HMM
#SBATCH -N 1
#SBATCH --partition=lrgmem
#SBATCH ---mem=1008GB
#SBATCH --time=72:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=kweave23@jhu.edu
#SBATCH --array=1-22%22
#SBATCH --account=rmccoy22

module load R/4.0.2
module load gcc/5.5.0

echo JobID_${SLURM_JOB_ID}_Sample_${SAMPLE_DIR}_chr${SLURM_ARRAY_TASK_ID} > out_assign_sperm_haplotypes_kw_20201220_${SAMPLE_DIR}_${SLURM_ARRAY_TASK_ID}.out
Rscript assign_sperm_haplotypes_kw.R ~/work-rmccoy22/rmccoy22/mccarroll-sperm-seq-data/${SAMPLE_DIR}/data/${SAMPLE_DIR}_goodcellsreplicatebcs_filteredhetsnps_${SLURM_ARRAY_TASK_ID}.cellsbyrow.txt ${SAMPLE_DIR}_chr${SLURM_ARRAY_TASK_ID} 0.005 2500 &>> out_assign_sperm_haplotypes_kw_20201220_${SAMPLE_DIR}_${SLURM_ARRAY_TASK_ID}.out
echo DONE >> out_assign_sperm_haplotypes_kw_20201220_${SAMPLE_DIR}_${SLURM_ARRAY_TASK_ID}.out
