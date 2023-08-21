#!/bin/bash

for RSD in 42 357 1848
do
  for SEQE in 0.001 0.005 0.05
  do
    for AVGR in 0.6 1 3
    do
      DIR=/home/kweave23/gamete_data/gen_model_results_noDNM/g${1}_s${2}_c${3}_se${SEQE}_r${AVGR}      
      date; time Rscript assess_with_rhapsodi_thread_arg.R $1 $2 $3 $SEQE $AVGR $RSD $4 &> $DIR/out_assess_${RSD}.txt
    done
  done
done
