#!/usr/bin/env Rscript


library("IHWpaper")


register(MulticoreParam(workers=20, verbose=TRUE))

#----------------- General benchmark settings -------------------------------#
alphas <- seq(0.01, .1, length=10)
nreps <- 4000 # number of times each simulation should be repeated (monte carlo replicates)

#------------- Simulation function ------------------------------------------#
ms <- 20000

sim_funs <- list(null_sim_fun(ms))

#------------- Methods to be benchmarked ------------------------------------#

continuous_methods_list <- list(bh,
							    lsl_gbh,
							    tst_gbh,
							    stratified_bh,
							    clfdr,
							    IHWpaper:::scott_fdrreg,
							    ddhf,
							    ihw_5fold)
# ihw_naive we run on another script on its

fdr_methods <- lapply(continuous_methods_list, continuous_wrap)


###############################################################################
eval_table <- run_evals(sim_funs, fdr_methods, nreps, alphas, BiocParallel=TRUE)

saveRDS(eval_table, file="result_files/ihw_null_simulation_benchmark.Rds")
