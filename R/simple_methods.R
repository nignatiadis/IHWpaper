#' bh: Wrapper for Benjamini Hochberg
#'
#' @param unadj_p  Numeric vector of unadjusted p-values.
#' @param alpha    Significance level at which to apply method
#'
#' @return BH multiple testing object
#'
#' @examples
#'      sim_df <- du_ttest_sim(20000,0.95, 1.5)
#'      obj <- bh(sim_df$pvalue, .1)
#'      sum(rejected_hypotheses(obj))
#'
#' @export
#' @importFrom stats p.adjust
bh <- function(unadj_p, alpha){
  	adj_p <- p.adjust(unadj_p, method="BH")
  	obj <- list(adj_p = adj_p, alpha = alpha)
  	class(obj) <- "BH"
	obj
}
attr(bh, "testing covariate") <- "simple" #i.e. covariate not considered
attr(bh, "fdr_method")        <- "BH"     

rejected_hypotheses.BH <- function(object, alpha= object$alpha){
	object$adj_p <= alpha
}

setOldClass("BH")
setMethod("rejected_hypotheses", signature("BH"), rejected_hypotheses.BH)


#' storey_qvalue: Wrapper for Storey's qvalue package
#'
#' @param unadj_p  Numeric vector of unadjusted p-values.
#' @param alpha    Significance level at which to apply method
#'
#' @return StoreyQValue multiple testing object
#'
#' @examples
#'      sim_df <- du_ttest_sim(20000,0.95, 1.5)
#'      obj <- storey_qvalue(sim_df$pvalue, .1)
#'      sum(rejected_hypotheses(obj))
#'
#' @export
#' @importFrom qvalue qvalue
storey_qvalue <- function(unadj_p, alpha){
	qval_res <- qvalue(unadj_p, alpha)
	obj <- list(adj_p = qval_res$qvalues, pi0=qval_res$pi0, alpha=alpha)
	class(obj) <- "StoreyQValue"
	obj
}
attr(storey_qvalue, "testing covariate") <- "simple" #i.e. covariate not considered
attr(storey_qvalue, "fdr_method")        <- "qvalue"     

rejected_hypotheses.StoreyQValue <- function(object, alpha= object$alpha){
	object$adj_p <= alpha
}

setOldClass("StoreyQValue")
setMethod("rejected_hypotheses", signature("StoreyQValue"), rejected_hypotheses.StoreyQValue)


#' bonf: Wrapper for Bonferroni
#'
#' @param unadj_p  Numeric vector of unadjusted p-values.
#' @param alpha    Significance level at which to apply method
#'
#' @return Bonferroni multiple testing object
#'
#' @examples
#'      sim_df <- du_ttest_sim(20000,0.95, 1.5)
#'      obj <- bonf(sim_df$pvalue, .1)
#'      sum(rejected_hypotheses(obj))
#'
#' @export
bonf <- function(unadj_p, alpha){
    adj_p <- p.adjust(unadj_p, method="bonferroni")
    obj <- list(adj_p = adj_p, alpha = alpha)
    class(obj) <- "Bonferroni"
    obj
}
attr(bonf, "testing covariate") <- "simple" #i.e. covariate not considered
attr(bonf, "fdr_method")        <- "Bonferroni"

rejected_hypotheses.Bonferroni <- function(object, alpha= object$alpha){
    object$adj_p <= alpha
}

setOldClass("Bonferroni")
setMethod("rejected_hypotheses", signature("Bonferroni"), rejected_hypotheses.Bonferroni)
