#!/usr/bin/env Rscript

# implement the below mainly for size-investing strategy illustration 
library("IHWpaper")
library("BiocParallel")
register(MulticoreParam(workers=10))

taus <- signif(10^seq(-4,-1, length=20),5)

methods_list <- list()

tau_wrap_ihwc <- function(tau){
  ihwc_method <- function(unadj_p, filterstat, alpha){
    obj <- ihw(unadj_p, filterstat, alpha, nbins=20, nfolds=5, lambdas=Inf,
               distrib_estimator="grenander", lp_solver="lpsymphony",
               censoring=TRUE, censoring_level=tau, null_proportion = FALSE)
    obj
  }
  
  
  attr(ihwc_method, "testing covariate") <- "continuous"
  attr(ihwc_method, "fdr_pars") <- tau
  attr(ihwc_method, "fdr_method") <- "IHWc"
  
  ihwc_method
}

tau_wrap_ihwc_storey <- function(tau){
  ihwc_storey_method <- function(unadj_p, filterstat, alpha){
    obj <- ihw(unadj_p, filterstat, alpha, nbins=20, nfolds=5, lambdas=Inf,
               distrib_estimator="grenander", lp_solver="lpsymphony",
               censoring=TRUE, censoring_level=tau, null_proportion = TRUE,
               null_proportion_level = 0.5)
    obj
  }
  
  attr(ihwc_storey_method, "testing covariate") <- "continuous"
  attr(ihwc_storey_method, "fdr_pars") <- tau
  attr(ihwc_storey_method, "fdr_method") <- "IHWc_Storey"
  ihwc_storey_method
}
  
for (i in 1:length(taus)){
  tau <- taus[i]
  methods_list <- c(methods_list, tau_wrap_ihwc(tau), tau_wrap_ihwc_storey(tau))
}




methods_list <- c(methods_list, ihw_5fold, ihw_storey_5fold, bh)

methods_list <- lapply(methods_list, continuous_wrap)

alphas <- .1
nreps <- 200# number of times each simulation should be repeated (monte carlo replicates)

#------------- Simulation function ------------------------------------------#
ms <- 20000

pi0s <- c(0.7,0.9)
xi_max <- 2.5

sim_funs <- lapply(pi0s, function(pi0) wasserman_normal_sim_fun(20000,pi0,0 , 2.5))

eval_table <- run_evals(sim_funs, methods_list, nreps, alphas, BiocParallel=TRUE)

saveRDS(eval_table, file="ihwc_wasserman_normal_sim.Rds")
