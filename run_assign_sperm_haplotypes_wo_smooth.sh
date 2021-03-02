#!/bin/bash

for NUM in {1..22}
do
  Rscript assign_sperm_haplotypes_wo_smooth.R /home/kweave23/gamete_data/filterhet/$1/${1}_euploid_${NUM}.txt /home/kweave23/gamete_data/misc_inputs/encode_blacklist.bed.gz /home/kweave23/gamete_data/misc_inputs/giab_union_regions.bed $1 $NUM /home/kweave23/gamete_data/run_bell_20210302/${1}/csv_out/ 5
done