#!/bin/bash

grep "The coverage of this simulation" $1 > $2/sim_coverage.txt

grep "Total number of recombination spots" $1 > $2/sim_recomb_spots_n.txt

grep "Number of de novo mutations" $1 > $2/sim_num_dnm.txt
grep "is filtered out" $1 >> $2/sim_num_dnm.txt

grep "new number of snps" $1 > $2/sim_num_snps_after_filter.txt

grep "Number of windows with overlap of" $1 >$2/sim_window_snp_tradeoff.txt

grep "Number of sperm affected for de novo mutation" $1 > $2/sim_num_sperm_affected.txt
grep "is filtered out" $1 >> $2/sim_num_sperm_affected.txt

grep "new rows vector" $1 | tail -n 1 > $2/sim_loc_dnms.txt
grep "is filtered out" $1 >> $2/sim_loc_dnms.txt

grep "affected sperm" $1 > $2/sim_affected_sperm.txt
grep "is filtered out" $1 >> $2/sim_affected_sperm.txt

grep -A 7 "Conservative recombination spot identification metrics" $1 > $2/sim_cons_recomb.txt
grep -A 7 "Liberal recombination spot identification metrics" $1 > $2/sim_lib_recomb.txt

grep "haplotype reconstruction accuracy" $1 > $2/sim_hap_reconstruction_acc.txt

grep "Time used for" $1 > $2/sim_inference_times.txt
