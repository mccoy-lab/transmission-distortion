#!/bin/bash


ml bedtools2/2.27.1

bedtools makewindows -g /home/kweave23/scr4_rmccoy22/kweave23/resources/hg38.chrom.sizes_main.txt -w 500000 > /home/kweave23/scr4_rmccoy22/kweave23/resources/windows.0pt5Mbp.bed

for CHROM in {1..22}
do
	for DONOR in ff3a ff4a nc10oldoil nc11ab nc12ab nc13ab nc14ab nc15ab nc16ab nc17ab nc18ab nc1abnov17 nc22abcd nc25abcd nc26abcd nc27aboct17 nc2absept17 nc3aboct17 nc4abnov17 nc6abcd nc8ab nc9ab pb2a pb3a pb4a
	do
		FILEOI=/home/kweave23/scr16_rmccoy22/kweave23/sc_transmission_distortion/runworkflow_20211105/${DONOR}/csv_out/assess_filtered/${DONOR}_chr${CHROM}_nonnapos.bed
		bedtools coverage -a /home/kweave23/scr4_rmccoy22/kweave23/resources/windows.0pt5Mbp.bed -b  $FILEOI -counts | grep chr${CHROM} > /home/kweave23/scr16_rmccoy22/kweave23/sc_transmission_distortion/runworkflow_20211105/${DONOR}/csv_out/assess_filtered/${DONOR}_${CHROM}_0pt5Mbp_coverage.txt 
	done
done

