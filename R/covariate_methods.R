
groups_by_filter <- function(filter_statistic, nbins){
  rfs <- rank(filter_statistic, ties.method="first")/length(filter_statistic)
  as.factor(ceiling( rfs* nbins))
}

### general wrapper
continuous_wrap <- function(mt_method, nbins=20){
  print(attr(mt_method,"fdr_method"))
  if (attr(mt_method, "testing covariate") == "simple"){

    f <- function(sim, alpha){
      mt_method(sim$pvalue, alpha)
    }
    attr(f, "fdr_method") <- attr(mt_method, "fdr_method")

  } else if (attr(mt_method, "testing covariate") == "stratified"){

    f <- function(sim, alpha){
      groups <- groups_by_filter(sim$filterstat, nbins)
      mt_method(sim$pvalue, groups, alpha)
    }
    attr(f, "fdr_method") <- paste(attr(mt_method, "fdr_method"), nbins, "bins")

  } else if (attr(mt_method, "testing covariate") == "continuous"){
    f <- function(sim, alpha){
      mt_method(sim$pvalue, sim$filterstat, alpha)
    }

    attr(f, "fdr_method") <- attr(mt_method, "fdr_method")

  } else {
    stop("unknown covariate handling")
  }

  f
}



#' scott_fdrreg: Wrapper for FDR regression (https://github.com/jgscott/FDRreg)
#'
#'
#' @param unadj_p  Numeric vector of unadjusted p-values.
#' @param filterstat   Factor to which different hypotheses belong
#' @param alpha    Significance level at which to apply method
#' @param df       Degrees of freedom for B-slines
#' @param lambda Ridge regularization parameter
#'
#' @return FDRreg multiple testing object
#'
#' @references  James G. Scott, Ryan C. Kelly, Matthew A. Smith, Pengcheng Zhou, and Robert E. Kass. 
#'         "False discovery rate regression: application to neural synchrony detection in primary visual cortex." 
#'         Journal of the American Statistical Association (2015).
#' @export
scott_fdrreg <- function(unadj_p, filterstat, alpha, df=3, lambda=0.01){
	# no automated way to choose function space over which we optimize
	# so we just use bs(df=3) as done in their analysis
	# also no automated way of choosing ridge regularization parameter

	b <- bs(filterstat, df=df)
	Xs <- model.matrix( ~  b - 1)

	fdrreg_res <- FDRreg(qnorm(unadj_p), Xs, control=list(lambda = 0.01))

  # modification to make this method applicable to p-values
  # FDRreg tests in a two-sided way, this means that we might get rejections
  # due to the right tail (high p-values) which have a z-score > 0 !
  # for this reason, just set local-fdr to 1 or equivalently the posterior prob to 0

  posterior_probs <- ifelse(unadj_p <= 0.5, fdrreg_res$postprob, 1)

	FDR <- getFDR(fdrreg_res$postprob)$FDR

  obj <- list(adj_p = FDR,
              pi0_continuous = fdrreg_res$priorprob, pi0=fdrreg_res$p0,
              method="FDRreg", alpha=alpha)
  class(obj) <- "FDRreg"
	obj
}

attr(scott_fdrreg, "testing covariate") <- "continuous"
attr(scott_fdrreg, "fdr_method")        <- "FDRreg"

#' @importFrom IHW rejected_hypotheses
rejected_hypotheses.FDRreg <- function(object, alpha= object$alpha){
	object$adj_p <= alpha
}

setOldClass("FDRreg")
setMethod("rejected_hypotheses", signature("FDRreg"), rejected_hypotheses.FDRreg)


#' ddhf: Greedy independent filtering
#'
#' @param unadj_p  Numeric vector of unadjusted p-values.
#' @param filterstat   Factor to which different hypotheses belong
#' @param alpha    Significance level at which to apply method
#'
#' @return DDHF multiple testing object
#'
#' @export
#' @useDynLib IHWpaper
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
attr(ddhf, "fdr_method")        <- "Greedy Indep. Filt."     
# actually greedy independent filtering does not need continuous covariates,
# it is sufficient if covariates can be ordered!

rejected_hypotheses.DDHF <- function(object, alpha= object$alpha){
	object$adj_p <= alpha
}

setOldClass("DDHF")
setMethod("rejected_hypotheses", signature("DDHF"), rejected_hypotheses.DDHF)

