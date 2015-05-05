du_ttest_simulation <- function(m, pi0, effect_size, n_samples=10){
  m0 <- ceiling(m*pi0)
  false_nulls <- sample(1:m, m-m0)
  z_table <- matrix(rnorm(n_samples*m), ncol=n_samples)
  z_table[false_nulls,(n_samples/2+1):n_samples] <- matrix(rnorm(n_samples/2*(m-m0), effect_size), ncol=n_samples/2)
  H <- rep(0,m)
  H[false_nulls] <- 1
  ttests <- rowttests(z_table, factor(rep(1:2,each=n_samples/2)))
  sds <- rowSds(z_table)
  filter_pvalue <- 1-pchisq((n_samples-1)*sds^2,n_samples-1)
  simDf <- data.frame(H=H, pvalue= ttests$p.value, sd= sds, filterstat=sds, 
    filter_pvalue=filter_pvalue, random_filterstat=runif(m))
  simDf
}

du_ttest <- function(m, pi0, effect_size, n_samples=10){
  f <- function() du_ttest_simulation(m, pi0, effect_size, n_samples=n_samples)
  attr(f, "sim_method") <- "t-test"
  attr(f, "sim_pars") <- paste0("pi0:",pi0, ", effect size:", effect_size)
  attr(f, "m") <- m
  f
}