#!/usr/bin/env Rscript

# implement the below mainly for size-investing strategy illustration 
library("IHWpaper")

register(MulticoreParam(workers=5, verbose=TRUE))

#----------------- General benchmark settings -------------------------------#
alphas <- .1#c(0.1, 0.01)
nreps <- 500 # number of times each simulation should be repeated (monte carlo replicates)

#------------- Simulation function ------------------------------------------#
ms <- 20000

xi_maxs<- seq(3, 6, length=10)

sim_funs <- lapply(xi_maxs, function(x) wasserman_normal_sim_fun(20000,0.9,1, x))



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
eval_table <- run_evals(sim_funs, fdr_methods, nreps, alphas, BiocParallel=T)
eval_table$xi_max = sapply(strsplit(eval_table$sim_pars,"xi_max:"),
                           function(x) as.numeric(x[2]))

saveRDS(eval_table, file="result_files/ihw_wasserman_normal_simulation_benchmark.Rds")
