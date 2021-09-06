#!/bin/bash

#SBATCH --job-name=full_filter
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

echo JobID_${SLURM_JOB_ID}_Sample_${SAMPLE_DIR}_chr${SLURM_ARRAY_TASK_ID} > ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_full_filter_20210825/${SAMPLE_DIR}/out_kw_20210825_${SAMPLE_DIR}_${SLURM_ARRAY_TASK_ID}.out
date; time Rscript assign_sperm_haplotypes_cr_full_filter.R ~/work-rmccoy22/scarios1/bell_filtered_data/${SAMPLE_DIR}/${SAMPLE_DIR}_euploid_${SLURM_ARRAY_TASK_ID}.txt ~/work-rmccoy22/kweave23/sc_transmission_distortion/filter_input_files/encode_blacklist.bed.gz ~/work-rmccoy22/kweave23/sc_transmission_distortion/filter_input_files/giab_union_regions.bed ~/work-rmccoy22/kweave23/sc_transmission_distortion/filter_input_files/GRCh38_notinalldifficultregions.bed ~/work-rmccoy22/kweave23/sc_transmission_distortion/filter_input_files/1kgp_positions_autosomes.txt ${SAMPLE_DIR} ${SLURM_ARRAY_TASK_ID} ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_full_filter_20210825/${SAMPLE_DIR}/csv_out/ 8 &>> ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_full_filter_20210825/${SAMPLE_DIR}/out_kw_20210825_${SAMPLE_DIR}_${SLURM_ARRAY_TASK_ID}.out
