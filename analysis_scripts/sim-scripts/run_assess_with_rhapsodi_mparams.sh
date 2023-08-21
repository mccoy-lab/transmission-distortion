#!/bin/bash
for COV in 0.001 0.01 0.1 0.223 0.357 0.511 0.693 1.204 2.303
do
  for RSD in 42 357 1848
  do
    for SEQE in 0.001 0.005 0.05
    do
      for AVGR in 0.6 1 3
      do
      date; time Rscript assess_with_rhapsodi_mparams.R $1 $2 $COV $SEQE $AVGR $RSD $3 &>> out_assess_${1}_${2}_${COV}_${SEQE}_${AVGR}_match_mparam_${RSD}.txt
      done
    done
  done
done
