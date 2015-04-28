# this file includes all methods which operate in a stratified fashion,
# i.e. as input they expect the raw p--values and a corresponding factor.

#' gbh: Grouped Benjamini Hochberg
#'
#'
#' @param unadj_p  Numeric vector of unadjusted p-values.
#' @param groups   Factor to which different hypotheses belong
#' @param alpha    Significance level at which to apply method
#' @param method   What pi0 estimator should be used (available "TST","LSL")
#' @param pi0_global GBH requires also a pi0 estimate for the marginal p-value distribution. Can either apply pi0 estimation method to
#'    all p-values (pi0_global="global") or use a weigted averarage (pi0_global="weighted_average") of the pi0 estimates within each stratum.
#'    This is not explicitly stated in the paper, but based on a reproduction of their paper figures it seems to be the weighted_average.

gbh <- function(unadj_p, groups, alpha, method="TST", pi0_global="weighted_average"){

  # special care has to be taken for TST GBH, uses alpha/(1+alpha) instead of alpha

  groups <- as.factor(groups)
  pv_list <- split(unadj_p, groups)

  if (method == "TST"){
    pi0_fun <- function(pv) tst_pi0_est(pv, alpha/(1+alpha)) # as described in GBH paper for this case
  } else if (method == "LSL") {
    pi0_fun <- lsl_pi0_est
  }

  m        <- length(unadj_p)
  m_groups <- sapply(pv_list, length)
  pi0_groups <- sapply(pv_list, pi0_fun)

  # in this case we reject no hypothesis
  if (all(pi0_groups ==1 )){
    pi0 <- 1
    adj_p <- rep(1,m)

  # do actual work
  } else {
    # get pi0 estimate (unstratified)
    if (pi0_global=="weighted_average"){
      pi0 <- sum(m_groups*pi0_groups)/m #average pi0 across all strata
    } else if (pi0_global=="global"){
      pi0 <- pi0_fun(unadj_p)
    }

    # helper functions
    my_mult <- function(a,b) ifelse(b==Inf, Inf, a*b)
    my_bh_adjust <- function(p, adjust_factor){      #standard step-up procedure
        i <- length(p):1L
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        pmin(1, cummin( my_mult(p[o], length(p)/i*adjust_factor))[ro])
    }

    weighted_pvals <- unsplit( mapply(function(pi0_g, pv) my_mult(pi0_g/(1-pi0_g), pv), pi0_groups, pv_list,  SIMPLIFY=FALSE), groups)

    adjust_factor <- if (method == "TST"){
                         (1-pi0)*(1+alpha) # again to make sure that TST GBH is applied at alpha/(1+alpha)
                      } else if (method=="LSL") {
                         (1-pi0)
                      }

    adj_p <- my_bh_adjust(weighted_pvals, adjust_factor)
  }

  obj <- list(weighted_pvals=weighted_pvals, adj_p = adj_p,
              pi0_groups = pi0_groups, pi0=pi0,
              method=paste(method,"GBH"), alpha=alpha)
  class(obj) <- "GBH"
  obj
}

#' tst_gbh: wrapper for gbh with method="TST"
#' lsl_gbh: wrapper for gbh with method="LSL"
#'

tst_gbh <- function(unadj_p, groups, alpha, ...) gbh(unadj_p, groups, alpha, method="TST", ...)
lsl_gbh <- function(unadj_p, groups, alpha, ...) gbh(unadj_p, groups, alpha, method="LSL", ...)

rejected_hypotheses.GBH <- function(gbh_object, alpha= gbh_object$alpha){
  gbh_object$adj_p <= alpha
}


#' stratified_bh: Stratified Benjamini Hochberg
#'
#'
#' @param unadj_p  Numeric vector of unadjusted p-values.
#' @param groups   Factor to which different hypotheses belong
#' @param alpha    Significance level at which to apply method
#'
#' @references Sun, Lei, et al. "Stratified false discovery control for large-scale hypothesis testing with application to genome-wide
#'    association studies." Genetic epidemiology 30.6 (2006): 519-530.
#' @references Yoo, Yun J., et al. "Were genome-wide linkage studies a waste of time? Exploiting candidate regions within genome-wide
#'    association studies." Genetic epidemiology 34.2 (2010): 107-118.

stratified_bh <- function(unadj_p, groups, alpha){
    groups <- as.factor(groups)
    pv_list <- split(unadj_p, groups)

    m        <- length(unadj_p)
    m_groups <- sapply(pv_list, length)

    adj_pv_list <- lapply(pv_list, function(pv) p.adjust(pv, method="BH"))
    adj_p <- unsplit(adj_pv_list, groups)

    obj <- list(adj_p = adj_p, alpha=alpha)
    class(obj) <- "SBH"
    obj
}

rejected_hypotheses.SBH <- function(obj, alpha= obj$alpha){
  obj$adj_p <= alpha
}
