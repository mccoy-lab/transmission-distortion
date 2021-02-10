#!/bin/bash

function submit_it {
  sbatch slurm_extract_sim_info.sh $1 $2
}

BASE_SAMPLENAME=simulation
BASE_CHROM=chrT
BASE_SEQERROR=0.05
BASE_WINDOW=3000
BASE_NGAM=1000
BASE_NSNP=30000
BASE_NNNA=2
BASE_RANDSD=42
BASE_RLAM=1
BASE_NBP=2
BASE_NBPR=0.09
BASE_ASE=TRUE
BASE_ADNM=TRUE
BASE_SEA=0.05
BASE_DNL=5
BASE_DNA=7.5
BASE_DNB=10

BASE_OUTDIR=${BASE_SAMPLENAME}_${BASE_CHROM}_${BASE_RANDSD}_${BASE_SEQERROR}_${BASE_ASE}_${BASE_SEA}_${BASE_WINDOW}_${BASE_NGAM}_${BASE_NSNP}_${BASE_NNNA}_${BASE_RLAM}_${BASE_NBP}_${BASE_NBPR}_${BASE_ADNM}_${BASE_DNL}_${BASE_DNA}_${BASE_DNB}/

submit_it ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${BASE_OUTDIR}out_kw_sim.txt ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${BASE_OUTDIR}
for RANDSD in 38 567 1008 2921
do
  NEW_DIR=${BASE_SAMPLENAME}_${BASE_CHROM}_${RANDSD}_${BASE_SEQERROR}_${BASE_ASE}_${BASE_SEA}_${BASE_WINDOW}_${BASE_NGAM}_${BASE_NSNP}_${BASE_NNNA}_${BASE_RLAM}_${BASE_NBP}_${BASE_NBPR}_${BASE_ADNM}_${BASE_DNL}_${BASE_DNA}_${BASE_DNB}/
  submit_it ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}out_kw_sim.txt ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}
done

for NGAM in 5 15 50 200 500 2500 5000 10000 16000
do
  NEW_DIR=${BASE_SAMPLENAME}_${BASE_CHROM}_${BASE_RANDSD}_${BASE_SEQERROR}_${BASE_ASE}_${BASE_SEA}_${BASE_WINDOW}_${NGAM}_${BASE_NSNP}_${BASE_NNNA}_${BASE_RLAM}_${BASE_NBP}_${BASE_NBPR}_${BASE_ADNM}_${BASE_DNL}_${BASE_DNA}_${BASE_DNB}/
  submit_it ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}out_kw_sim.txt ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}
  for RANDSD in 38 567 1008 2921
  do
    NEW_DIR=${BASE_SAMPLENAME}_${BASE_CHROM}_${RANDSD}_${BASE_SEQERROR}_${BASE_ASE}_${BASE_SEA}_${BASE_WINDOW}_${NGAM}_${BASE_NSNP}_${BASE_NNNA}_${BASE_RLAM}_${BASE_NBP}_${BASE_NBPR}_${BASE_ADNM}_${BASE_DNL}_${BASE_DNA}_${BASE_DNB}/
    submit_it ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}out_kw_sim.txt ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}
  done
done

for NNNA in 3 4 5 10 25 50 100
do
  NEW_DIR=${BASE_SAMPLENAME}_${BASE_CHROM}_${BASE_RANDSD}_${BASE_SEQERROR}_${BASE_ASE}_${BASE_SEA}_${BASE_WINDOW}_${BASE_NGAM}_${BASE_NSNP}_${NNNA}_${BASE_RLAM}_${BASE_NBP}_${BASE_NBPR}_${BASE_ADNM}_${BASE_DNL}_${BASE_DNA}_${BASE_DNB}/
  submit_it ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}out_kw_sim.txt ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}
  for RANDSD in 38 567 1008 2921
  do
    NEW_DIR=${BASE_SAMPLENAME}_${BASE_CHROM}_${RANDSD}_${BASE_SEQERROR}_${BASE_ASE}_${BASE_SEA}_${BASE_WINDOW}_${BASE_NGAM}_${BASE_NSNP}_${NNNA}_${BASE_RLAM}_${BASE_NBP}_${BASE_NBPR}_${BASE_ADNM}_${BASE_DNL}_${BASE_DNA}_${BASE_DNB}/
    submit_it ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}out_kw_sim.txt ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}
  done
done

for NBP in 0 4 6 8 16
do
  NEW_DIR=${BASE_SAMPLENAME}_${BASE_CHROM}_${BASE_RANDSD}_${BASE_SEQERROR}_${BASE_ASE}_${BASE_SEA}_${BASE_WINDOW}_${BASE_NGAM}_${BASE_NSNP}_${BASE_NNNA}_${BASE_RLAM}_${NBP}_${BASE_NBPR}_${BASE_ADNM}_${BASE_DNL}_${BASE_DNA}_${BASE_DNB}/
  submit_it ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}out_kw_sim.txt ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}
  for RANDSD in 38 567 1008 2921
  do
    NEW_DIR=${BASE_SAMPLENAME}_${BASE_CHROM}_${RANDSD}_${BASE_SEQERROR}_${BASE_ASE}_${BASE_SEA}_${BASE_WINDOW}_${BASE_NGAM}_${BASE_NSNP}_${BASE_NNNA}_${BASE_RLAM}_${NBP}_${BASE_NBPR}_${BASE_ADNM}_${BASE_DNL}_${BASE_DNA}_${BASE_DNB}/
    submit_it ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}out_kw_sim.txt ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}
  done
done

for NBPR in 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2 0.4
do
  NEW_DIR=${BASE_SAMPLENAME}_${BASE_CHROM}_${BASE_RANDSD}_${BASE_SEQERROR}_${BASE_ASE}_${BASE_SEA}_${BASE_WINDOW}_${BASE_NGAM}_${BASE_NSNP}_${BASE_NNNA}_${BASE_RLAM}_${BASE_NBP}_${NBPR}_${BASE_ADNM}_${BASE_DNL}_${BASE_DNA}_${BASE_DNB}/
  submit_it ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}out_kw_sim.txt ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}
  for RANDSD in 38 567 1008 2921
  do
    NEW_DIR=${BASE_SAMPLENAME}_${BASE_CHROM}_${RANDSD}_${BASE_SEQERROR}_${BASE_ASE}_${BASE_SEA}_${BASE_WINDOW}_${BASE_NGAM}_${BASE_NSNP}_${BASE_NNNA}_${BASE_RLAM}_${BASE_NBP}_${NBPR}_${BASE_ADNM}_${BASE_DNL}_${BASE_DNA}_${BASE_DNB}/
    submit_it ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}out_kw_sim.txt ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}
  done
done

for WINDOW in 1000 1500 2000 2500 3500 4000 5000 5500
do
  NEW_DIR=${BASE_SAMPLENAME}_${BASE_CHROM}_${BASE_RANDSD}_${BASE_SEQERROR}_${BASE_ASE}_${BASE_SEA}_${WINDOW}_${BASE_NGAM}_${BASE_NSNP}_${BASE_NNNA}_${BASE_RLAM}_${BASE_NBP}_${BASE_NBPR}_${BASE_ADNM}_${BASE_DNL}_${BASE_DNA}_${BASE_DNB}/
  submit_it ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}out_kw_sim.txt ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}
  for RANDSD in 38 567 1008 2921
  do
    NEW_DIR=${BASE_SAMPLENAME}_${BASE_CHROM}_${RANDSD}_${BASE_SEQERROR}_${BASE_ASE}_${BASE_SEA}_${WINDOW}_${BASE_NGAM}_${BASE_NSNP}_${BASE_NNNA}_${BASE_RLAM}_${BASE_NBP}_${BASE_NBPR}_${BASE_ADNM}_${BASE_DNL}_${BASE_DNA}_${BASE_DNB}/
    submit_it ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}out_kw_sim.txt ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}
  done
done

for SEA in 0.01 0.025 0.075 0.1 0.2
do
  NEW_DIR=${BASE_SAMPLENAME}_${BASE_CHROM}_${BASE_RANDSD}_${BASE_SEQERROR}_${BASE_ASE}_${SEA}_${BASE_WINDOW}_${BASE_NGAM}_${BASE_NSNP}_${BASE_NNNA}_${BASE_RLAM}_${BASE_NBP}_${BASE_NBPR}_${BASE_ADNM}_${BASE_DNL}_${BASE_DNA}_${BASE_DNB}/
  submit_it ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}out_kw_sim.txt ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}
  for RANDSD in 38 567 1008 2921
  do
    NEW_DIR=${BASE_SAMPLENAME}_${BASE_CHROM}_${RANDSD}_${BASE_SEQERROR}_${BASE_ASE}_${SEA}_${BASE_WINDOW}_${BASE_NGAM}_${BASE_NSNP}_${BASE_NNNA}_${BASE_RLAM}_${BASE_NBP}_${BASE_NBPR}_${BASE_ADNM}_${BASE_DNL}_${BASE_DNA}_${BASE_DNB}/
    submit_it ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}out_kw_sim.txt ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim_20210202/${NEW_DIR}
  done
done
