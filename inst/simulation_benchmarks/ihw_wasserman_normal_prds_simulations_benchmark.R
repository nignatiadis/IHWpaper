#!/usr/bin/env Rscript

library("IHWpaper")

args <- commandArgs(trailingOnly=TRUE)
idx <- as.integer(args[1])

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

register(MulticoreParam(workers=9))

#----------------- General benchmark settings -------------------------------#
alphas <- 0.1
nreps <- 400 # number of times each simulation should be repeated (monte carlo replicates)

#------------- Simulation function ------------------------------------------#
ms <- 20000

param_df <- expand.grid(rho=seq(0,0.5,length=5),
                        latent_factors=c(1,10),
                        pi0=c(0.9,0.95,1),
                        xi_max = c(3,4))

# grab parameters for current run
rho <- param_df$rho[idx]
latent_factors <- param_df$latent_factors[idx]
pi0 <- param_df$pi0[idx]
xi_max <- param_df$xi_max[idx]

sim_funs <- list(wasserman_normal_prds_sim_fun(ms, pi0,
                    rho=rho, latent_factors=latent_factors, xi_min=0, xi_max=xi_max))

#------------- Methods to be benchmarked ------------------------------------#
continuous_methods_list <- list(bh,
                                ihw_5fold_reg)

fdr_methods <- lapply(continuous_methods_list, continuous_wrap)


#-----------------------------------------------------------------------------

eval_table <- run_evals(sim_funs, fdr_methods, nreps, alphas, BiocParallel=TRUE)

extract_params <- function(x){
  comma_split <- strsplit(x,", ")[[1]]
  par_split <- lapply(comma_split, function(presplit) strsplit(presplit,":")[[1]])
  df <- do.call(cbind, lapply(par_split, function(tmp) {tmp_df <-data.frame(as.numeric(tmp[2]));
                              names(tmp_df) <- tmp[1]; tmp_df}))
  df
}

eval_table <- eval_table %>%
            do(cbind(.,extract_params(.$sim_pars))) %>%
            mutate(idx=idx)

file_name <- file.path("result_files","ihw_wasserman_normal_prds",
                        paste0("par_setting_",idx,".Rds"))
saveRDS(eval_table, file=file_name)
