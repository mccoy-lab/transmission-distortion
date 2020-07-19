# This pipeline will: 
# download vcf from 1000 genomes project 
# grab an individual (here we use bcftools view to get a vcf file of just one person; we chose one of European descent) 
# choose only (using biallelic SNPs bcftools view); this creates an IMPUTE2 format file. IMPUTE2 is simpler and removes much of the info from vcf files - useful for any problem with imputation
# In the IMPUTE2 format, we care about a) the legend (gives column descriptions) and b) the haplotype table (includes zeros and ones for each parental haplotype)
# 	In haplotype table, col1 is hap1, col2 is hap2, and the rows are the SNPs
# 	To create recombination, siply introduce a break at a given row and swap the top of hap1 with the top of hap2

# To download hts.lib, bcf tools: http://www.htslib.org/download/ (details incl. in github wiki)


#!/bin/bash
# use wget to download the file for the chromosome you'd like to simulate 
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr21.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz 
# paper: https://wellcomeopenresearch.org/articles/4-50/v2
# choose individual (here, NA06985) and output it to a single vcf file 
../bcftools-1.10.2/bin/bcftools view --output-type v --samples NA06985 ../vcf_data/ALL.chr21.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz --output-file NA06985.vcf 

# step necessary for use of hts-lib 
../htslib-1.10.2/bin/bgzip -c NA06985.vcf > NA06985.vcf.gz
../htslib-1.10.2/bin/tabix -f -p vcf NA06985.vcf.gz

# Options to include only biallelic SNPs; detailed info on each option available in github wiki 
../bcftools-1.10.2/bin/bcftools view \
NA06985.vcf.gz \
--samples NA06985 \
--exclude-types indels,mnps,ref,bnd,other \
--min-alleles 2 \
--max-alleles 2 \
--min-ac 1 \
--phased \
--exclude 'AN!=2*N_SAMPLES' \
--output-file NA06985_chr21_EUR_panel.bcf \
--output-type u
		
# Convert to donor file that offer legend and hap tables; details available in github wiki 		
../bcftools-1.10.2/bin/bcftools convert NA06985_chr21_EUR_panel.bcf --haplegendsample NA06985_donor

gunzip NA06985_donor.legend.gz
gunzip NA06985_donor.hap.gz


