# get number of SNPs from each donor for bootstrap

for SAMPLE_DIR in nc1abnov17 nc2absept17 nc3aboct17 nc4abnov17 nc6abcd nc8ab nc10oldoil nc9ab nc11ab nc12ab nc13ab nc14ab nc15ab nc16ab nc17ab nc18ab nc22abcd nc25abcd nc26abcd nc27aboct17 ff3a ff4a pb2a pb3a pb4a
do 
  echo ${SAMPLE_DIR} >> num_snps.txt
  cat ./${SAMPLE_DIR}/*prune.in | wc -l >> num_snps.txt
done 
