#!/bin/bash

for SAMPLE_DIR in ff3a ff4a pb2a pb3a pb4a nc10oldoil nc11ab nc12ab nc13ab nc14ab nc15ab nc16ab nc17ab nc18ab nc1abnov17 nc22abcd nc25abcd nc26abcd nc27aboct17 nc2absept17 nc3aboct17 nc4abnov17 nc6abcd nc8ab nc9ab
do
  mkdir -p /home-3/kweave23@jhu.edu/work-rmccoy22/kweave23/sc_transmission_distortion/runworkflow_20211105/${SAMPLE_DIR}/csv_out
  sbatch --export=SAMPLE_DIR=$SAMPLE_DIR slurm_filter_rhapsodi_TDscan.sh
done
