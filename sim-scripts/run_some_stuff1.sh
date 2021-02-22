#!/bin/bash

mkdir run_sim3_20210218/$1/
for RSDVAL in 27 42 386 651 1059 2556 2563 2862 3417 4900
do
  date; time Rscript test_metrics.R 2 3000 1000 30000 0.01 $RSDVAL 1 &> run_sim3_20210218/$1/rsd${RSDVAL}.txt
  grep "Total number of recombination" $1/rsd${RSDVAL}.txt > $1/sim3_nrecomb_${RSDVAL}.txt
  grep "new number of snps" $1/rsd${RSDVAL}.txt > $1/sim3_nsnps_${RSDVAL}.txt
  grep "haplotype reconstruction accuracy" $1/rsd${RSDVAL}.txt > $1/sim3_hap_rec_${RSDVAL}.txt
  grep -A 7 "Conservative recombination spot identification metrics" $1/rsd${RSDVAL}.txt > $1/sim3_cons_recomb_${RSDVAL}.txt
  grep -A 7 "Liberal recombination spot identification metrics" $1/rsd${RSDVAL}.txt > $1/sim3_lib_recomb_${RSDVAL}.txt
done
