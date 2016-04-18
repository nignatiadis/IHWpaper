#!/usr/bin/env Rscript

library("IHWpaper")

args <- commandArgs(trailingOnly=TRUE)
idx <- as.integer(args[1])

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

register(MulticoreParam(workers=8))

#----------------- General benchmark settings -------------------------------#
alphas <- 0.1
nreps <- 1000 # number of times each simulation should be repeated (monte carlo replicates)

#------------- Simulation function ------------------------------------------#
ms <- 20000
#sim_funs <- list(null_sim_fun(ms))
eff_sizes <- seq(1, 2.5, length=20)


# grab first effect size
eff_sizes <- eff_sizes[idx]

sim_funs <- lapply(eff_sizes, function(x) du_ttest_sim_fun(20000,0.95,x, uninformative_filter=F))

#------------- Methods to be benchmarked ------------------------------------#
continuous_methods_list <- list(ihw_bonf_5fold_reg, bonf)

fdr_methods <- lapply(continuous_methods_list, continuous_wrap)


#-----------------------------------------------------------------------------

eval_table <- run_evals(sim_funs, fdr_methods, nreps, alphas, BiocParallel=TRUE)
eval_table <- mutate(eval_table,
        eff_size = sapply(strsplit(eval_table$sim_pars,"effect size:"),
                function(x) as.numeric(x[2])))
eval_table <- mutate(eval_table, idx=idx)

file_name <- file.path("result_files","ihw_bonf_du_ttest_informative",
                        paste0("par_setting_",idx,".Rds"))
saveRDS(eval_table, file=file_name)
