library(rhapsodi)

args = commandArgs(trailingOnly=TRUE)

num_gametes <-  as.integer(args[1])
coverage <-  as.numeric(args[2])
rs <-  as.integer(args[3])
fname <- paste0("simOut/g", num_gametes, "_c", coverage, "_rs", rs, ".RData")

load(fname)

times <- system.time(
{
  rhapsodi_autorun(NULL, use_dt=TRUE, input_dt=sim_out$gam_na, threads=48)
})

outName <- paste0("times/g", num_gametes, "_c", coverage, "_rs", rs, "_out.csv")

write.csv(data.frame(t(data.matrix(times))), outName, row.names=FALSE)
