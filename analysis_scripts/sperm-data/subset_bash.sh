# This script accesses each file from the McCarroll lab and executes a filtering script to 
# keep only the sperm cells that were included in Bell's autosomalploidy file AND were euploid. `subset_data.R`


#!/bin/bash
for file_name in ~/work/rmccoy22/mccarroll-sperm-seq-data/nc9*/data/*cellsbyrow.txt; do
    individual=$(echo $file_name | cut -d"_" -f1)
    individual=$(echo $individual | cut -d"/" -f9)
    chromosome=$(echo $file_name | cut -d"_" -f4)
    chromosome=$(echo $chromosome | cut -d"." -f1)
    output="${individual}_euploid_${chromosome}.txt"
    ./subset_data.R $file_name $output $individual $chromosome
done