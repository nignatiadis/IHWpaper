#!/usr/bin/env Rscript

library("IHWpaper")

args <- commandArgs(trailingOnly=TRUE)
idx <- as.integer(args[1])

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

register(MulticoreParam(workers=12))

#----------------- General benchmark settings -------------------------------#
alphas <- 0.1
nreps <- 500  #number of times each simulation should be repeated (monte carlo replicates)

#------------- Simulation function ------------------------------------------#
ms <- 20000

xi_maxs<- seq(3, 6, length=10)
xi_maxs <- xi_maxs[idx]

sim_funs <- lapply(xi_maxs, function(x) wasserman_normal_sim_fun(20000,0.9,1, x))



#------------- Methods to be benchmarked ------------------------------------#
continuous_methods_list <- list(ihw_bonf_5fold_reg, bonf)

fdr_methods <- lapply(continuous_methods_list, continuous_wrap)


#-----------------------------------------------------------------------------
eval_table <- run_evals(sim_funs, fdr_methods, nreps, alphas, BiocParallel=T)
eval_table$xi_max = sapply(strsplit(eval_table$sim_pars,"xi_max:"),
                function(x) as.numeric(x[2]))

file_name <- file.path("result_files","ihw_bonf_wasserman_normal",
                        paste0("par_setting_",idx,".Rds"))
saveRDS(eval_table, file=file_name)
