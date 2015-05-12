#' scott_fdrreg: Wrapper for FDR regression
#'
#'
#' @param unadj_p  Numeric vector of unadjusted p-values.
#' @param filterstat   Factor to which different hypotheses belong
#' @param alpha    Significance level at which to apply method

scott_fdrreg <- function(unadj_p, filterstat, alpha, df=3, lambda=0.01){
	# no automated way to choose function space over which we optimize
	# so we just use bs(df=3) as done in their analysis
	# also no automated way of choosing ridge regularization parameter

	b <- bs(filterstat, df=df)
	Xs <- model.matrix( ~  b - 1)

	fdrreg_res <- FDRreg(qnorm(unadj_p), Xs, control=list(lambda = 0.01))

	FDR <- getFDR(fdrreg_res$postprob)$FDR
  	
  	obj <- list(adj_p = FDR,
              pi0_continuous = fdrreg_res$priorprob, pi0=fdrreg_res$p0,
              method="FDRreg", alpha=alpha)
  	class(obj) <- "FDRreg"
	obj
}
attr(scott_fdrreg, "testing covariate") <- "continuous" 
attr(scott_fdrreg, "fdr_method")        <- "FDRreg"     

rejected_hypotheses.FDRreg <- function(obj, alpha= obj$alpha){
	obj$adj_p <= alpha
}



#' ddhf: Greedy independent filtering
#'
#'
#' @param unadj_p  Numeric vector of unadjusted p-values.
#' @param filterstat   Factor to which different hypotheses belong
#' @param alpha    Significance level at which to apply method


ddhf <- function(unadj_p, filterstat, alpha){
  sortedp <- sort(unadj_p)
  ranksp <- rank(unadj_p, ties.method="first")
  filterorder <- order(filterstat)
  filterrank <- rank(filterstat)

  sortedp_filterorder <- ranksp[filterorder]
  rjs <- ddhf_order(sortedp, sortedp_filterorder, alpha)
  j <- which.max(rjs)[1]

  m <- length(unadj_p)
  adj_p <- rep(1, m)
  adj_p[filterrank >= j] <- p.adjust(unadj_p[filterrank >= j], method="BH")

  ddhf_obj <- list(adj_p = adj_p, rejections = max(rjs),
                   cutoff_quantile=j/m,
                   cutoff_value = filterstat[filterorder][j],
                   alpha=alpha)
  class(ddhf_obj) <- "DDHF"
  ddhf_obj
}
attr(ddhf, "testing covariate") <- "continuous" 
attr(ddhf, "fdr_method")        <- "DDHF"     
# actually greedy independent filtering does not need continuous covariates,
# it is sufficient if covariates can be ordered!

rejected_hypotheses.DDHF <- function(obj, alpha= obj$alpha){
	obj$adj_p <= alpha
}

