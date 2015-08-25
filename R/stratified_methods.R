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

#' @references  Hu, James X., Hongyu Zhao, and Harrison H. Zhou. "False discovery rate control with groups." 
#'         Journal of the American Statistical Association 105.491 (2010).

gbh <- function(unadj_p, groups, alpha, method="TST", pi0_global="weighted_average"){

  # special care has to be taken for TST GBH, uses alpha/(1+alpha) instead of alpha

  groups <- as.factor(groups)
  pv_list <- split(unadj_p, groups)

  if (method == "TST"){
    pi0_fun <- function(pv) tst_pi0_est(pv, alpha) # as described in GBH paper for this case
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
    weighted_pvals <- unadj_p
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
                         (1-pi0)*(1+alpha) # to make sure that TST GBH is applied at alpha/(1+alpha)
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
attr(gbh, "testing covariate") <- "stratified" # i.e. covariates can be considered by stratifying based on them
attr(gbh, "fdr_method")        <- "GBH"

#' tst_gbh: wrapper for gbh with method="TST"
#' lsl_gbh: wrapper for gbh with method="LSL"
#'

tst_gbh <- function(unadj_p, groups, alpha, ...) gbh(unadj_p, groups, alpha, method="TST", ...)
attr(tst_gbh, "testing covariate") <- "stratified" # i.e. covariates can be considered by stratifying based on them
attr(tst_gbh, "fdr_method")        <- "TST GBH"

lsl_gbh <- function(unadj_p, groups, alpha, ...) gbh(unadj_p, groups, alpha, method="LSL", ...)
attr(lsl_gbh, "testing covariate") <- "stratified" # i.e. covariates can be considered by stratifying based on them
attr(lsl_gbh, "fdr_method")        <- "LSL GBH"

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

attr(stratified_bh, "testing covariate") <- "stratified" # i.e. covariates can be considered by stratifying based on them
attr(stratified_bh, "fdr_method")        <- "SBH"

rejected_hypotheses.SBH <- function(object, alpha= object$alpha){
    object$adj_p <= alpha
}


#' clfdr: Cai's local fdr based method
#'
#'
#' @param unadj_p  Numeric vector of unadjusted p-values.
#' @param groups   Factor to which different hypotheses belong
#' @param alpha    Significance level at which to apply method
#' @param lfdr_estimation  Method used to estimate the loca fdr, defaults to fdrtool

#' @references Cai, T. Tony, and Wenguang Sun. "Simultaneous testing of grouped hypotheses: Finding needles in multiple haystacks." 
#'           Journal of the American Statistical Association 104.488 (2009).

clfdr <- function(unadj_p, groups, alpha, lfdr_estimation="fdrtool"){

  # estimate local fdr within each stratum first

  lfdr_res <- lfdr_fit(unadj_p, groups, lfdr_estimation=lfdr_estimation)
  lfdrs <- lfdr_res$lfdr

  # now use the rejection rule described in Cai's paper

  # Remark:
  # When sorting lfdrs, we break ties by pvalues so that in the end within each stratum
  # we get monotonic adjusted p-values as a function of the p-values
  # This is mainly needed for grenander based lfdrs, with most other
  # lfdr estimation methods lfdr ties are not a problem usually

  o <- order(lfdrs, unadj_p)
  lfdrs_sorted <- lfdrs[o]
  fdr_estimate <- cumsum(lfdrs_sorted)/(1:length(unadj_p))
  adj_p <- rev(cummin(rev(fdr_estimate)))
  adj_p <- adj_p[order(o)]


  obj <- list(adj_p = adj_p, alpha=alpha, lfdr_estimation=lfdr_estimation)
  class(obj) <- "Clfdr"
  obj
}

attr(clfdr, "testing covariate") <- "stratified" # i.e. covariates can be considered by stratifying based on them
attr(clfdr, "fdr_method")        <- "Clfdr"

rejected_hypotheses.Clfdr <- function(obj, alpha= obj$alpha){
    obj$adj_p <= alpha
}

# helper function for cai
lfdr_fit <- function(unadj_p, group, lfdr_estimation="fdrtool"){
  if (lfdr_estimation == "covmod"){

    lfdrs <- covmod_grouped(unadj_p, group)$cm_fdr

  } else {

      pvals_list <- split(unadj_p, group)

      if (lfdr_estimation == "fdrtool"){
        lfdr_fun <- function(pv) fdrtool(pv, statistic="pvalue",plot=F,verbose=F)$lfdr
      } else if (lfdr_estimation == "ConcaveFDR"){
        lfdr_fun <- function(pv) ConcaveFDR(pv, statistic="pvalue",plot=F,verbose=F)$lfdr.log
      } else if (lfdr_estimation == "locfdr"){
        lfdr_fun <- function(pv) locfdr(qnorm(pv), nulltype=0, plot=0)$fdr
      } else if (lfdr_estimation == "mixfdr"){
        lfdr_fun <- function(pv) mixFdr(qnorm(pv), theonull=T, plots=F)$fdr
      }

      lfdr_list <- lapply(pvals_list, lfdr_fun)
      lfdrs <- unsplit(lfdr_list, group)
    }
    fit_obj <- data.frame(pvalue=unadj_p, lfdr=lfdrs, group=group)
    fit_obj
}
