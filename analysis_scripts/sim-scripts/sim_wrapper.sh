#!/bin/bash

#SBATCH --job-name=sim
#SBATCH --nodes=1
#SBATCH -p shared
#SBATCH --ntasks-per-node=24
#SBATCH --time=2:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=rmccoy22@jhu.edu
#SBATCH --array=1-1000%20

ml gcc
ml R

iteration=${SLURM_ARRAY_TASK_ID}

cd /scratch/groups/rmccoy22/rmccoy22/spermseq_sim
Rscript null_sim.R ${iteration} /scratch/groups/rmccoy22/rmccoy22/spermseq_sim/output_slurm

# combine SLURM outputs and write to null_sim.txt