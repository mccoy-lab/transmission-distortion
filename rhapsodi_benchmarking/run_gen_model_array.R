library(rhapsodi)

args = commandArgs(trailingOnly=TRUE)

num_gametes <-  as.integer(args[1])
coverage <-  as.numeric(args[2])
rs <-  as.integer(args[3])

sim_out <- sim_run_generative_model(num_gametes, 100000, coverage, 1, random_seed=rs)
fname <- paste0("simOut/g", num_gametes, "_c", coverage, "_rs", rs, ".RData")
save(sim_out, file=fname)
