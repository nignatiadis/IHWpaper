#!/usr/bin/env Rscript

library("ddhwPaper")


register(MulticoreParam(workers=10, verbose=TRUE))

#----------------- General benchmark settings -------------------------------#
alphas <- .1#c(0.1, 0.01)
nreps <- 500 # number of times each simulation should be repeated (monte carlo replicates)

#------------- Simulation function ------------------------------------------#
ms <- 20000
#sim_funs <- list(null_sim_fun(ms))
eff_sizes <- seq(1, 2.5, length=10)
sim_funs <- lapply(eff_sizes, function(x) du_ttest_sim_fun(20000,0.95,x, uninformative_filter=F))

#------------- Methods to be benchmarked ------------------------------------#
continuous_methods_list <- list(bh,
							    storey_qvalue,
							    lsl_gbh,
							    tst_gbh,
							    stratified_bh,
							    cai,
							    scott_fdrreg,
							    ddhf,
							    ddhw_lp,
							    ddhw_lp_unregularized)

fdr_methods <- lapply(continuous_methods_list, continuous_wrap)


#-----------------------------------------------------------------------------
dir.create("result_files/du_ttest_informative_simulations_benchmark/")
eval_table <- run_evals(sim_funs, fdr_methods, nreps, alphas, BiocParallel=T, print_dir = "result_files/du_ttest_informative_simulations_benchmark/")
eval_table <- mutate(eval_table,
		eff_size = sapply(strsplit(eval_table$sim_pars,"effect size:"),
				function(x) as.numeric(x[2])))

saveRDS(eval_table, file="result_files/du_ttest_informative_simulation_benchmark_2.Rds")
