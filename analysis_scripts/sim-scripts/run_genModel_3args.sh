#!/bin/bash

for RSD in 42 357 1848
do
  for SEQE in 0.001 0.005 0.05
  do
    for AVGR in 0.6 1 3
    do
      mkdir -p /home/kweave23/gamete_data/gen_model_results_noDNM/g${1}_s${2}_c${3}_se${SEQE}_r${AVGR}
      date; time Rscript generative_model_for_rhapsodi.R  /home/kweave23/gamete_data/gen_model_results_noDNM/g${1}_s${2}_c${3}_se${SEQE}_r${AVGR}/ $1 $2 $3 $RSD $AVGR $SEQE &> /home/kweave23/gamete_data/gen_model_results_noDNM/g${1}_s${2}_c${3}_se${SEQE}_r${AVGR}/out_rs${RSD}.txt
    done
  done
done