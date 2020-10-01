#!/bin/bash
for file_name in ~/work/rmccoy22/mccarroll-sperm-seq-data/nc*/data/*cellsbyrow.txt; do
    individual=$(echo $file_name | cut -d"_" -f1)
    individual=$(echo $individual | cut -d"/" -f9)
    chromosome=$(echo $file_name | cut -d"_" -f4)
    chromosome=$(echo $chromosome | cut -d"." -f1)
    output="${individual}_output_${chromosome}.txt"
    ./genotype_tables.R $file_name $output
done