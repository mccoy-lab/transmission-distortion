ml r/4.0.2

arg=$(sed "${SLURM_ARRAY_TASK_ID}q;d" sim_params.txt)
ngam=$(echo $arg | cut -f1 -d' ')
cov=$(echo $arg | cut -f2 -d' ')
rs=$(echo $arg | cut -f3 -d' ')

Rscript run_gen_model_array.R $ngam $cov $rs
