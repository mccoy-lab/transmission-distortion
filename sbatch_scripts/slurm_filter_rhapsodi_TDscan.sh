#!/bin/bash

#SBATCH --job-name=runworkflow
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

INPUTDIR=/home-3/kweave23@jhu.edu/work-rmccoy22/scarios1/bell_filtered_data/${SAMPLE_DIR}
BASEDIR=/home-3/kweave23@jhu.edu/work-rmccoy22/kweave23/sc_transmission_distortion
DIR=${BASEDIR}/runworkflow_20211105/${SAMPLE_DIR}
OUTDIR=${DIR}/csv_out/
FILTERDIR=${BASEDIR}/filter_input_files

echo JobID_${SLURM_JOB_ID}_Sample_${SAMPLE_DIR}_chr${SLURM_ARRAY_TASK_ID} > ${DIR}/out_kw_runworkflow_${SAMPLE_DIR}_${SLURM_ARRAY_TASK_ID}.out
date; time Rscript filter_rhapsodi_TDscan.R ${INPUTDIR}/${SAMPLE_DIR}_euploid_${SLURM_ARRAY_TASK_ID}.txt ${FILTERDIR}/encode_blacklist.bed.gz ${FILTERDIR}/giab_union_regions.bed ${FILTERDIR}/1kgp_positions_autosomes.txt ${SAMPLE_DIR} ${SLURM_ARRAY_TASK_ID} ${OUTDIR} 8 &>> ${DIR}/out_kw_runworkflow_${SAMPLE_DIR}_${SLURM_ARRAY_TASK_ID}.out
