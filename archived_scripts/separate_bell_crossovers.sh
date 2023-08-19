#!/bin/bash


for CHR in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22
	do
		awk '{if ($3 == "'$CHR'") {print}}' $1 > $2_$CHR.txt
	done


# call from command line with /Users/saracarioscia/mccoy-lab/transmission-distortion/donor_1$ ../separate_bell_crossovers.sh ./donor_1.txt donor1