### general wrapper

#' Benchmarking wrapper: Given a multiple testing method, convert it so that it takes
#'   a simulation object (see simulation function) and a nominal level alpha
#'   as inputs
#'
#' @param mt_method Multiple testing method (e.g. a function such as gbh or ddhf)
#' @param nbins Integer, number of equally sized bins into which to stratify hypotheses
#'
#' @return A new multiple testing function which has an
#'           interface of the form f(sim_data_frame, alpha)
#' @examples
#'      sim_df <- du_ttest_sim(20000,0.95, 1.5)
#'      sim_df$group <- groups_by_filter(sim_df$filterstat, 20)
#'      obj <- tst_gbh(sim_df$pvalue, sim_df$group, .1)
#'      sum(rejected_hypotheses(obj))
#'
#'      tst_gbh_continuous <- continuous_wrap(tst_gbh)
#'      obj2 <- tst_gbh_continuous(sim_df, .1)
#'      sum(rejected_hypotheses(obj2))
#'
#' @export
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

scott_fdrreg <- function(unadj_p, filterstat, alpha, df=3, lambda=0.01){
	# no automated way to choose function space over which we optimize
	# so we just use bs(df=3) as done in their analysis
	# also no automated way of choosing ridge regularization parameter

  if (!requireNamespace("FDRreg", quietly=TRUE)){
     stop("FDRreg package required for this function to work.")
  }

  if (! packageVersion("FDRreg") %in% c("0.2.1","0.2")){
     stop(paste("Benchmarks were run against version 0.2 of FDRreg",
                "available on github via:",
                "devtools::install_github(repo= 'jgscott/FDRreg', subdir='R_pkg/', ref = 'a63cebae6faecb1fb0ebee634195296f39faa11b')"
                ))
  }
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
#' @examples
#'      sim_df <- du_ttest_sim(20000,0.95, 1.5)
#'      obj <- ddhf(sim_df$pvalue, sim_df$filterstat, .1)
#'      sum(rejected_hypotheses(obj))
#'
#' @export
#' @importFrom Rcpp evalCpp
#' @importFrom Rcpp sourceCpp
#' @useDynLib IHWpaper, .registration=TRUE
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

