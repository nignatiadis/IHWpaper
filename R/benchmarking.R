# general functions to make it easy to benchmark FDR methods on given simulations..


run_evals <- function(sim_funs, fdr_methods, nreps, alphas,...){
	rbind_all(lapply(sim_funs, function(x) sim_fun_eval(x, fdr_methods, nreps, alphas, ...)))
}

sim_fun_eval <- function(sim_fun, fdr_methods, nreps, alphas, BiocParallel=T, ...){
	sim_seeds <- 1:nreps
	if (BiocParallel){
		evaluated_sims <- bplapply(sim_seeds, function(x) sim_eval(sim_fun, x, fdr_methods, alphas,...))
	} else {
		evaluated_sims <- lapply(sim_seeds, function(x) sim_eval(sim_fun, x, fdr_methods, alphas,...))
	}
	df <- dplyr::rbind_all(evaluated_sims)
	df <- dplyr::summarize(group_by(df, fdr_method, fdr_pars, alpha), FDR = mean(FDP),
		 power= mean(power), rj_ratio = mean(rj_ratio), FPR = mean(FPR), FDR_sd = sd(FDP), 
		 nsuccessful = n())
	m  <- attr(sim_fun, "m")
	sim_method <- attr(sim_fun, "sim_method")
	sim_pars <- attr(sim_fun, "sim_pars")
	df <- mutate(df, sim_method = sim_method, m=  m, sim_pars = sim_pars)
	df
}

sim_eval <- function(sim_fun, seed, fdr_methods, alphas, print_dir = NULL){
	print(paste("@@@@@ seed is ", seed))
	if (!is.null(print_dir)){
		sim_pars <- attr(sim_fun, "sim_pars")
		file_name <- paste0(print_dir, "seed", seed,"_", sim_pars, ".Rds")
		if (file.exists(file_name)){
			df <- readRDS(file_name)
		} else {
			sim <- sim_fun(seed) 
			df <- rbind_all(lapply(fdr_methods, function(fdr_method) sim_fdrmethod_eval(sim, fdr_method, alphas)))
			saveRDS(df, file=file_name)
		}
	} else {
		sim <- sim_fun(seed) 
		df <- rbind_all(lapply(fdr_methods, function(fdr_method) sim_fdrmethod_eval(sim, fdr_method, alphas)))
	}
	df
}

sim_fdrmethod_eval <- function(sim, fdr_method, alphas){
	rbind_all(lapply(alphas, function(a) sim_alpha_eval(sim, fdr_method, a)))
}

sim_alpha_eval <- function(sim, fdr_method, alpha){

	# tryCatch block mainly because many of the local fdr methods will occasionally crash
	df <- tryCatch({

    	df <- calculate_test_stats(sim, fdr_method, alpha)

		if (is.null(attr(fdr_method,"fdr_pars"))){
			fdr_pars <- NA
		} else{
			fdr_pars <- attr(fdr_method, "fdr_pars")
		}	
		df <- mutate(df, alpha= alpha, 
			fdr_pars = fdr_pars,
			fdr_method=attr(fdr_method,"fdr_method"))
		df
		
	}, error = function(e) {
    	 df <- data.frame()
    	 df		
		})
	df
}

calculate_test_stats <- function(sim, fdr_method, alpha){
	fdr_method_result <- fdr_method(sim, alpha)
	rejected <- rejected_hypotheses(fdr_method_result)
	rjs <- sum(rejected)
	rj_ratio <- rjs/nrow(sim)
	FDP <- ifelse(rjs == 0, 0, sum(sim$H == 0 & rejected)/rjs)
	power <- ifelse(sum(sim$H) == 0, NA, sum(sim$H == 1 & rejected)/sum(sim$H ==1))
	# in principle I should take special care in case only alternatives, but we are not
	# interested in this scenario...
	FPR <-  sum(sim$H == 0 & rejected)/sum(sim$H == 0) 

	df <- data.frame(rj_ratio=rj_ratio, FDP=FDP, power=power, FPR=FPR)
	df
}

