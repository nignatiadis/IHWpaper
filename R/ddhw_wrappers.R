#diverse DDHW wrappers (e.g. LP relaxation, adaptive, SD=0 or SD=1)

# also look again again at DDHW MILP unreg, DDHW MILP adaptive etc...? <- probably not

ddhw_lp_adaptive  <- function(unadj_p, filterstat, alpha){
	obj <- ddhw_auto(unadj_p, filterstat, alpha , 
				solver="Rsymphony" , nbins=20, reg_par_number=10, pi0_adaptive=T, lp_relaxation=T, threads=1, 
				optim_pval_threshold = alpha, test_distrib="both", nreps=20, sd_penalty=1)
	class(obj) <- "ddhw"
	obj
}
attr(ddhw_lp_adaptive, "testing covariate") <- "continuous" 
attr(ddhw_lp_adaptive, "fdr_method")        <- "DDHW lp adaptive"     

ddhw_lp  <- function(unadj_p, filterstat, alpha){
	obj <- ddhw_auto(unadj_p, filterstat, alpha , 
				solver="Rsymphony" , nbins=20, reg_par_number=10, pi0_adaptive=F, lp_relaxation=T, threads=1,
				optim_pval_threshold = alpha, test_distrib="both", nreps=20, sd_penalty=0)
	class(obj) <- "ddhw"
	obj
}

attr(ddhw_lp, "testing covariate") <- "continuous" 
attr(ddhw_lp, "fdr_method")        <- "DDHW lp"     



ddhw_lp_unregularized  <- function(unadj_p, filterstat, alpha){
	obj <- ddhw(unadj_p, filterstat, 20, alpha , 
				solver="Rsymphony" ,  pi0_adaptive=F, lp_relaxation=T, threads=1,
				optim_pval_threshold = alpha)
	class(obj) <- "ddhw"
	obj
}

attr(ddhw_lp_unregularized, "testing covariate") <- "continuous" 
attr(ddhw_lp_unregularized, "fdr_method")        <- "DDHW lp unregularized"     


# gurobi versions of above funs:

ddhw_lp_grb  <- function(unadj_p, filterstat, alpha){
	obj <- ddhw_auto(unadj_p, filterstat, alpha , 
				solver="gurobi" , nbins=20, reg_par_number=10, pi0_adaptive=F, lp_relaxation=T, threads=1,
				 test_distrib="both", nreps=20, sd_penalty=0, optim_pval_threshold=alpha)
	class(obj) <- "ddhw"
	obj
}

attr(ddhw_lp_grb, "testing covariate") <- "continuous" 
attr(ddhw_lp_grb, "fdr_method")        <- "DDHW lp"     



ddhw_lp_unregularized_grb <- function(unadj_p, filterstat, alpha){
	obj <- ddhw(unadj_p, filterstat, 20, alpha , 
				solver="gurobi" , pi0_adaptive=F, lp_relaxation=T, threads=1,
				 optim_pval_threshold=alpha)
	class(obj) <- "ddhw"
	obj
}

attr(ddhw_lp_unregularized_grb, "testing covariate") <- "continuous" 
attr(ddhw_lp_unregularized_grb, "fdr_method")        <- "DDHW lp unregularized"     

###########################

rejected_hypotheses.ddhw <- function(obj, alpha= obj@alpha){
	adj_pvalues(obj) <= alpha
}


