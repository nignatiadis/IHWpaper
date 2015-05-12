#' bh: Wrapper for Benjamini Hochberg
#'
#' @param unadj_p  Numeric vector of unadjusted p-values.
#' @param alpha    Significance level at which to apply method

bh <- function(unadj_p, alpha){
  	adj_p <- p.adjust(unadj_p, method="BH")
  	obj <- list(adj_p = adj_p, alpha = alpha)
  	class(obj) <- "BH"
	obj
}
attr(bh, "testing covariate") <- "simple" #i.e. covariate not considered
attr(bh, "fdr_method")        <- "BH"     

rejected_hypotheses.BH <- function(obj, alpha= obj$alpha){
	obj$adj_p <= alpha
}


#' storey_qvalue: Wrapper for Storey's qvalue package
#'
#' @param unadj_p  Numeric vector of unadjusted p-values.
#' @param alpha    Significance level at which to apply method

storey_qvalue <- function(unadj_p, alpha){
	qval_res <- qvalue(unadj_p, alpha)
	obj <- list(adj_p = qval_res$qvalues, pi0=qval_res$pi0, alpha=alpha)
	class(obj) <- "StoreyQValue"
	obj
}
attr(storey_qvalue, "testing covariate") <- "simple" #i.e. covariate not considered
attr(storey_qvalue, "fdr_method")        <- "qvalue"     

rejected_hypotheses.StoreyQValue <- function(obj, alpha= obj$alpha){
	obj$adj_p <= alpha
}