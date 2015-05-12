#diverse DDHW wrappers (e.g. LP relaxation, adaptive, SD=0 or SD=1)

# also look again again at DDHW MILP unreg, DDHW MILP adaptive etc...? <- probably not

ddhw_lp_adaptive  <- function(unadj_p, filterstat, alpha){
	obj <- ddhw_auto2(unadj_p, filterstat, alpha , 
				solver="Rsymphony" , reg_par_number=20, pi0_adaptive=T, lp_relaxation=T, threads=1)
	class(obj) <- "DDHWwrapper"
	obj
}
attr(ddhw_lp_adaptive, "testing covariate") <- "continuous" 
attr(ddhw_lp_adaptive, "fdr_method")        <- "DDHW lp adaptive"     

ddhw_lp  <- function(unadj_p, filterstat, alpha){
	obj <- ddhw_auto2(unadj_p, filterstat, alpha , 
				solver="Rsymphony" , reg_par_number=20, pi0_adaptive=F, lp_relaxation=T, threads=1)
	class(obj) <- "DDHWwrapper"
	obj
}

attr(ddhw_lp, "testing covariate") <- "continuous" 
attr(ddhw_lp, "fdr_method")        <- "DDHW lp"     


ddhw_lp_adaptive_nosd  <- function(unadj_p, filterstat, alpha){
	obj <- ddhw_auto2(unadj_p, filterstat, alpha , 
				solver="Rsymphony" , reg_par_number=20, pi0_adaptive=T, lp_relaxation=T, threads=1, sd_penalty=0)
	class(obj) <- "DDHWwrapper"
	obj
}
attr(ddhw_lp_adaptive_nosd, "testing covariate") <- "continuous" 
attr(ddhw_lp_adaptive_nosd, "fdr_method")        <- "DDHW lp adaptive nosd"     



ddhw_lp_unregularized  <- function(unadj_p, filterstat, alpha){
	obj <- ddhw(unadj_p, filterstat, 20, alpha , 
				solver="Rsymphony" , reg_par_number=20, pi0_adaptive=F, lp_relaxation=T, threads=1)
	class(obj) <- "DDHWwrapper"
	obj
}

attr(ddhw_lp_unregularized, "testing covariate") <- "continuous" 
attr(ddhw_lp_unregularized, "fdr_method")        <- "DDHW lp unregularized"     


rejected_hypotheses.DDHWwrapper <- function(obj, alpha= alpha.ddhw(obj)){
	adj_pvalues.ddhw(obj) <= alpha
}


