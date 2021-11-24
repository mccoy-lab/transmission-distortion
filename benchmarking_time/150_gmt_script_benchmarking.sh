#!/bin/bash
#SBATCH --partition=shared
#SBATCH --time=72:0:0
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --account=rmccoy22
#SBATCH --mail-type=end
#SBATCH --mail-user=kweave23@jhu.edu


arg=$(sed "${SLURM_ARRAY_TASK_ID}q;d" 150_gamete_sims.txt)
time Rscript rhapsodi_benchmark_nothread.R ${arg} ${SLURM_ARRAY_TASK_ID}
