
ddhw_unregularized <- function(unadj_p, filterstat, alpha){
	obj <- ddhw(unadj_p, filterstat, alpha, nbins=20, nfolds=1, lambdas=Inf)
	obj
}

ddhw_5fold <- function(unadj_p, filterstat, alpha){
	obj <- ddhw(unadj_p, filterstat, alpha, nbins=20, nfolds=5, lambdas=Inf)
	obj
}

attr(ddhw_unregularized, "testing covariate") <- "continuous"
attr(ddhw_unregularized, "fdr_method")        <- "DDHW unregularized"

attr(ddhw_5fold, "testing covariate") <- "continuous"
attr(ddhw_5fold, "fdr_method")        <- "DDHW 5fold"


rejected_hypotheses.ddhwResult <- function(object, alpha= object@alpha){
  adj_pvalues(object) <= alpha
}

###########################

