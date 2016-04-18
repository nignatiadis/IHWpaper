#!/usr/bin/env Rscript

library("IHWpaper")


register(MulticoreParam(workers=5))

#----------------- General benchmark settings -------------------------------#
alphas <- 0.1
nreps <- 1000 # number of times each simulation should be repeated (monte carlo replicates)

#------------- Simulation function ------------------------------------------#
ms <- 20000
#sim_funs <- list(null_sim_fun(ms))
eff_sizes <- seq(1, 2.5, length=20)
sim_funs <- lapply(eff_sizes, function(x) du_ttest_sim_fun(ms,0.95,x, uninformative_filter=F))

#------------- Methods to be benchmarked ------------------------------------#
continuous_methods_list <- list(bh,
                                lsl_gbh,
                                tst_gbh,
                                stratified_bh,
                                clfdr,
                                IHWpaper:::scott_fdrreg,
                                ddhf,
                                ihw_5fold)


fdr_methods <- lapply(continuous_methods_list, continuous_wrap)


#-----------------------------------------------------------------------------

eval_table <- run_evals(sim_funs, fdr_methods, nreps, alphas, BiocParallel=TRUE)
eval_table <- mutate(eval_table,
                     eff_size = sapply(strsplit(eval_table$sim_pars,"effect size:"),
                                       function(x) as.numeric(x[2])))

saveRDS(eval_table, file="result_files/ihw_du_ttest_inform_simulation_benchmark.Rds")
