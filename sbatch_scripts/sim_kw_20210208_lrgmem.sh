#!/bin/bash

#SBATCH --job-name=sim2
#SBATCH -N 1
#SBATCH --partition=lrgmem
#SBATCH --cpus-per-task=8
#SBATCH --mem=1008000M
#SBATCH --time=72:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=kweave23@jhu.edu
#SBATCH --account=rmccoy22

module load R/4.0.2
module load gcc/5.5.0
module load bedtools/2.27.0
module load htslib/1.9

export PATH=$PATH:/home-3/kweave23@jhu.edu/code/bin

echo JobID_${SLURM_JOB_ID} > ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim2_20210208/${OUTDIR}out_kw_sim.txt
Rscript sim_varying_coverage_num_sperm.R $SAMPLENAME $CHROM $OUTDIR 8 $SEQERROR $WINDOW $NGAM $NSNP $COVERAGE $RANDSD $RLAM $ASE $ADNM $SEA $DNL $DNA $DNB &>> ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim2_20210208/${OUTDIR}out_kw_sim.txt
bash extract_sim_info.sh ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim2_20210208/${OUTDIR}out_kw_sim.txt ~/work-rmccoy22/kweave23/sc_transmission_distortion/run_sim2_20210208/${OUTDIR}
