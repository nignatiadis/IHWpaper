#!/usr/bin/env Rscript


library("ddhwPaper")


register(MulticoreParam(workers=2, verbose=TRUE))

#----------------- General benchmark settings -------------------------------#
alphas <- seq(0.02, .1, length=5)
nreps <- 10 # number of times each simulation should be repeated (monte carlo replicates)

#------------- Simulation function ------------------------------------------#
ms <- 20000

sim_funs <- list(null_sim_fun(ms))

#------------- Methods to be benchmarked ------------------------------------#
source("continuous_covariate_methods_collection.R")
fdr_methods <- wrapped_continuous_methods_list

###############################################################################3
eval_table <- run_evals(sim_funs, fdr_methods, nreps, alphas, BiocParallel=T)

saveRDS(eval_table, file="result_files/null_simulation_benchmark.Rds")
