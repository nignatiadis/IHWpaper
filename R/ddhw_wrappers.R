
# ddhw_unregularized <- function(unadj_p, filterstat, alpha){
# 	obj <- ddhw(unadj_p, filterstat, alpha, nbins=20, nfolds=1, lambdas=Inf)
# 	obj
# }

# ddhw_5fold <- function(unadj_p, filterstat, alpha){
# 	obj <- ddhw(unadj_p, filterstat, alpha, nbins=20, nfolds=5, lambdas=Inf)
# 	obj
# }

# ddhw_5fold_reg <- function(unadj_p, filterstat, alpha){
# 	obj <- ddhw(unadj_p, filterstat, alpha, nbins=20, nfolds=5)
# 	obj
# }

# attr(ddhw_unregularized, "testing covariate") <- "continuous"
# attr(ddhw_unregularized, "fdr_method")        <- "DDHW unregularized"

# attr(ddhw_5fold, "testing covariate") <- "continuous"
# attr(ddhw_5fold, "fdr_method")        <- "DDHW 5fold"


# attr(ddhw_5fold_reg, "testing covariate") <- "continuous"
# attr(ddhw_5fold_reg, "fdr_method")        <- "DDHW 5fold reg"

ihw_naive <- function(unadj_p, filterstat, alpha){
 	obj <- ihw(unadj_p, filterstat, alpha, nbins=20, nfolds=1, lambdas=Inf,
 		distrib_estimator="ECDF", lp_solver="gurobi")
 	obj
}

attr(ihw_naive, "testing covariate") <- "continuous"
attr(ihw_naive, "fdr_method")        <- "IHW naive"

ihw_ecdf_5fold <- function(unadj_p, filterstat, alpha){
	obj <- ihw(unadj_p, filterstat, alpha, nbins=20, nfolds=5, lambdas=Inf,
 		distrib_estimator="ECDF", lp_solver="gurobi")
 	obj
}

attr(ihw_ecdf_5fold, "testing covariate") <- "continuous"
attr(ihw_ecdf_5fold, "fdr_method")        <- "IHW ECDF"


ihw_5fold <- function(unadj_p, filterstat, alpha){
 	obj <- ihw(unadj_p, filterstat, alpha, nbins=20, nfolds=5, lambdas=Inf,
 		distrib_estimator="grenander", lp_solver="lpsymphony")
 	obj
}

attr(ihw_5fold, "testing covariate") <- "continuous"
attr(ihw_5fold, "fdr_method")        <- "IHW"


ihw_5fold_reg <- function(unadj_p, filterstat, alpha){
  obj <- ihw(unadj_p, filterstat, alpha, nbins=20, nfolds=5,
             distrib_estimator="grenander", lp_solver="lpsymphony")
  obj
}

attr(ihw_5fold_reg, "testing covariate") <- "continuous"
attr(ihw_5fold_reg, "fdr_method")        <- "IHW E3"



rejected_hypotheses.ihwResult <- function(object, alpha= object@alpha){
  adj_pvalues(object) <= alpha
}


###########################

