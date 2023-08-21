# Filtering Bell Data 

Filter data to only include: 
- Cells that passed Bell's filtering steps (i.e., were included in her `autosomalploidy.txt` file)
- Cells that were euploid for the chromosome of interest
- `filter_bell_cells.R`


Remove problematic regions of the genome based on: 
- Encode blacklist regions (exclude)
- Genome in a bottle "in union" regions (include) `GRCh38_notinalldifficultregions.bed`
- `filtering-td-data.R`

