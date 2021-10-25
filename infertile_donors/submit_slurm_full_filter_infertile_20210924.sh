#!/bin/bash

for SAMPLE_DIR in ff3a ff4a pb2a pb3a pb4a
do
  mkdir -p /home-3/kweave23@jhu.edu/work-rmccoy22/kweave23/sc_transmission_distortion/infertile_donors/${SAMPLE_DIR}/csv_out
  sbatch --export=SAMPLE_DIR=$SAMPLE_DIR slurm_full_filter_infertile_20210924.sh
done
