# implement LSL and TST pi0 estimators which we need in order to benchmark the GBH method


#' TST (Two-Step) pi0 estimator
#' @param pvalue Numeric vector of unadjusted p-values.
#' @param alpha Nominal level for applying the TST procedure
#'
#' @return estimated proportion of null hypotheses (pi0)
#'
#' @export
tst_pi0_est <- function(pvalue, alpha){
    nrjs <- sum(p.adjust(pvalue, method="BH") <= alpha/(1+alpha))
    n <- length(pvalue)
    (n-nrjs)/n
}

#' LSL (Least-Slope) pi0 estimator
#' @param pvalue Numeric vector of unadjusted p-values.
#'
#' @return estimated proportion of null hypotheses (pi0)
#'
#' @export
lsl_pi0_est <- function(pvalue){
  n <- length(pvalue)
  ls <- (n:1)/(1-sort(pvalue))
  ls_diff <- ls[-c(1,2)] - ls[-c(1,n)]
  index <- min(which(ls_diff > 0))+2
  if (index == Inf) {
    pi0 <- 1
  } else {
    pi0 <- min(1, (1+floor(ls[index]))/n)
  }
  pi0
}