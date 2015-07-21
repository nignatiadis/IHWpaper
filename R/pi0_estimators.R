# implement LSL and TST pi0 estimators which we need in order to benchmark the GBH method
# implementation of GBH in structSSI R package is unfortunately too buggy

tst_pi0_est <- function(pvalue, alpha){
    nrjs <- sum(p.adjust(pvalue, method="BH") <= alpha/(1+alpha))
    n <- length(pvalue)
    (n-nrjs)/n
}
