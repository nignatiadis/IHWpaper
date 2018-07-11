context("simulation_check")


# Exemplary simulation with low number of Monte Carlo replicates


#----------------- General benchmark settings -------------------------------#
alphas <- 0.1
nreps <- 2 # number of times each simulation should be repeated (monte carlo replicates)

#------------- Simulation function ------------------------------------------#
ms <- 20000
eff_sizes <- seq(1, 2.5, length=2)
sim_funs  <- lapply(eff_sizes, function(x) du_ttest_sim_fun(ms, 0.95, x, uninformative_filter = FALSE))

#------------- Methods to be benchmarked ------------------------------------#
continuous_methods_list <- list(bh,
                                lsl_gbh,
                                tst_gbh,
                                stratified_bh,
                                clfdr,
                                ddhf,
                                ihw_5fold,
                                ihw_bonf_5fold_reg,
                                bonf,
                                storey_qvalue)


fdr_methods <- lapply(continuous_methods_list, continuous_wrap)


#-----------------------------------------------------------------------------

eval_table <- run_evals(sim_funs, fdr_methods, nreps, alphas, BiocParallel = FALSE)
eval_table <- dplyr::mutate(eval_table,
                     eff_size = sapply(strsplit(sim_pars, "effect size:"), function(x) as.numeric(x[2])))

testthat::expect_equal(eval_table$nsuccessful, rep(nreps, nrow(eval_table)))
