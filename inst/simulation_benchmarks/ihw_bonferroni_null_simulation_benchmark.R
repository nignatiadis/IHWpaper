#!/usr/bin/env Rscript

library("IHWpaper")

args <- commandArgs(trailingOnly=TRUE)
idx <- as.integer(args[1])

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

register(MulticoreParam(workers=8))

#----------------- General benchmark settings -------------------------------#

alphas <- seq(0.01, .1, length=10)
nreps <- 4000 # number of times each simulation should be repeated (monte carlo replicates)

#------------- Simulation function ------------------------------------------#
ms <- 20000

sim_funs <- list(null_sim_fun(ms))


# grab corresponding alpha
alphas <- alphas[idx]


#------------- Methods to be benchmarked ------------------------------------#
continuous_methods_list <- list(ihw_bonf_5fold_reg, bonf)

fdr_methods <- lapply(continuous_methods_list, continuous_wrap)


#-----------------------------------------------------------------------------
eval_table <- run_evals(sim_funs, fdr_methods, nreps, alphas, BiocParallel=TRUE)

file_name <- file.path("result_files","ihw_bonf_null",
                        paste0("par_setting_",idx,".Rds"))

saveRDS(eval_table, file=file_name)











