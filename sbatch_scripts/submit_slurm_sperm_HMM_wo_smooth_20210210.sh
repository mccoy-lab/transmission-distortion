#!/bin/bash

for SAMPLE_DIR in nc1abnov17 nc2absept17 nc3aboct17 nc4abnov17 nc6abcd nc8ab nc9ab nc10oldoil nc11ab nc12ab nc13ab nc14ab nc15ab nc16ab nc17ab nc18ab nc22abcd nc25abcd nc26abcd nc27aboct17
do
  sbatch --export=SAMPLE_DIR=$SAMPLE_DIR slurm_sperm_HMM_wo_smooth_20210210.sh
done
