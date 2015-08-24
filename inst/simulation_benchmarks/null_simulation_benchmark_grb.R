#!/usr/bin/env Rscript
library("gurobi")

library("ddhwPaper")


register(MulticoreParam(workers=10, verbose=TRUE))

#----------------- General benchmark settings -------------------------------#
alphas <- seq(0.01, .1, length=10)
nreps <- 3000 # number of times each simulation should be repeated (monte carlo replicates)

#------------- Simulation function ------------------------------------------#
ms <- 20000

sim_funs <- list(null_sim_fun(ms))

#------------- Methods to be benchmarked ------------------------------------#

continuous_methods_list <- list(bh,
							    storey_qvalue,
							    lsl_gbh,
							    tst_gbh,
							    stratified_bh,
							    cai,
							    scott_fdrreg,
							    ddhf,
							    ddhw_lp_grb,
							    ddhw_lp_unregularized_grb)

fdr_methods <- lapply(continuous_methods_list, continuous_wrap)


###############################################################################3
#dir.create("result_files/null_simulation_benchmark_grb/")

eval_table <- run_evals(sim_funs, fdr_methods, nreps, alphas, BiocParallel=T, print_dir = "result_files/null_simulation_benchmark_grb/")

saveRDS(eval_table, file="result_files/null_simulation_benchmark_grb3000.Rds")
