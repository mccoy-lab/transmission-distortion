for SAMPLE_DIR in nc1abnov17 nc2absept17 nc3aboct17 nc4abnov17 nc6abcd nc8ab nc10oldoil nc9ab nc11ab nc12ab nc13ab nc14ab nc15ab nc16ab nc17ab 
nc18ab nc22abcd nc25abcd nc26abcd nc27aboct17 ff3a ff4a pb2a pb3a pb4a
do
  for chr in {1..22}
  do
    mkdir -p ./${SAMPLE_DIR}
    plink --file ~/work/kweave23/sc_transmission_distortion/run_make_ped/${SAMPLE_DIR}/${SAMPLE_DIR}_${chr}_full_sperm_converted --no-fid --no-parents --no-sex --no-pheno --indep 50 5 2 --out ${SAMPlE_DIR}_${chr}
  done
done
