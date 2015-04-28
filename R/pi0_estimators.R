# implement LSL and TST pi0 estimators which we need in order to benchmark the GBH method
# implementation of GBH in structSSI R package is unfortunately too buggy

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

tst_pi0_est <- function(pvalue, alpha){
    nrjs <- sum(p.adjust(pvalue, method="BH") <= alpha/(1+alpha))
    n <- length(pvalue)
    (n-nrjs)/n
}
