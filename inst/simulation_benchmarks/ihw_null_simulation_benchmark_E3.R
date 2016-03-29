#!/usr/bin/env Rscript


library("IHWpaper")


register(MulticoreParam(workers=5, verbose=TRUE))

#----------------- General benchmark settings -------------------------------#
alphas <- seq(0.01, .1, length=10)
nreps <- 4000 # number of times each simulation should be repeated (monte carlo replicates)

#------------- Simulation function ------------------------------------------#
ms <- 20000

sim_funs <- list(null_sim_fun(ms))

#------------- Methods to be benchmarked ------------------------------------#

continuous_methods_list <- list(ihw_5fold_reg)


fdr_methods <- lapply(continuous_methods_list, continuous_wrap)


###############################################################################
eval_table <- run_evals(sim_funs, fdr_methods, nreps, alphas, BiocParallel=TRUE)

saveRDS(eval_table, file="result_files/ihw_null_simulation_benchmark_E3.Rds")
