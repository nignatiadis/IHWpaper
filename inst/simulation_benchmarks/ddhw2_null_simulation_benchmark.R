#!/usr/bin/env Rscript


library("ddhwPaper")


register(MulticoreParam(workers=20, verbose=TRUE))

#----------------- General benchmark settings -------------------------------#
alphas <- seq(0.01, .1, length=10)
nreps <- 4000 # number of times each simulation should be repeated (monte carlo replicates)

#------------- Simulation function ------------------------------------------#
ms <- 20000

sim_funs <- list(null_sim_fun(ms))

#------------- Methods to be benchmarked ------------------------------------#

continuous_methods_list <- list(bh,
							    storey_qvalue,
							    lsl_gbh,
							    tst_gbh,
							    stratified_bh,
							    clfdr,
							    scott_fdrreg,
							    ddhf,
							    ddhw_5fold,
							    ddhw_unregularized)

fdr_methods <- lapply(continuous_methods_list, continuous_wrap)


###############################################################################3
eval_table <- run_evals(sim_funs, fdr_methods, nreps, alphas, BiocParallel=T)

saveRDS(eval_table, file="result_files/ddhw2_null_simulation_benchmark.Rds")
