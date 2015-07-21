du_ttest_sim <- function(m, pi0, effect_size, n_samples=10, uninformative_filter=F, seed=NULL){
  if (!is.null(seed)) set.seed(seed)
  m0 <- ceiling(m*pi0)
  false_nulls <- sample(1:m, m-m0)
  z_table <- matrix(rnorm(n_samples*m), ncol=n_samples)
  z_table[false_nulls,(n_samples/2+1):n_samples] <- matrix(rnorm(n_samples/2*(m-m0), effect_size), ncol=n_samples/2)
  H <- rep(0,m)
  H[false_nulls] <- 1
  ttests <- rowttests(z_table, factor(rep(1:2,each=n_samples/2)))
  sds <- rowSds(z_table)
  filter_pvalue <- 1-pchisq((n_samples-1)*sds^2,n_samples-1)
  if (uninformative_filter){
    filterstats <- runif(m)
  } else {
    filterstats <- sds
  }

  simDf <- data.frame(H=H, pvalue= ttests$p.value, sd= sds, filterstat=filterstats, 
    filter_pvalue=filter_pvalue)
  simDf
}

du_ttest_sim_fun <- function(m, pi0, effect_size, n_samples=10, uninformative_filter=F){
  f <- function(seed) du_ttest_sim(m, pi0, effect_size, n_samples=n_samples, uninformative_filter=uninformative_filter, seed=seed)
  attr(f, "sim_method") <- "t-test"
  attr(f, "sim_pars") <- paste0("pi0:",pi0, ", effect size:", effect_size)
  attr(f, "m") <- m
  f
}

null_sim <- function(m, seed=NULL){
  if (!is.null(seed)) set.seed(seed)

  sim_df <- data.frame(pvalue = runif(m), filterstat = runif(m), H=0)
  return(sim_df)
}

null_sim_fun <- function(m){
  f <- function(seed) null_sim(m=m, seed=seed)
  attr(f, "sim_method") <- "only nulls"
  attr(f, "sim_pars") <- paste0("")
  attr(f, "m") <- m
  f
}


wasserman_normal_sim <- function(m, pi0, xi_min, xi_max, seed=NULL){
    if (!is.null(seed)) set.seed(seed)

    X   <- runif(m, min=xi_min, max=xi_max)
    H   <- rbinom(m,1,1-pi0)
    Z   <- rnorm(m, H*X)
    pvalue <- 1-pnorm(Z)
    simDf <- data.frame(pvalue=pvalue, filterstat=X,H=H, Z=Z)
}

wasserman_normal_sim_fun <- function(m, pi0, xi_min, xi_max){
  f <- function(seed) wasserman_normal_sim(m, pi0, xi_min, xi_max, seed=seed)
  attr(f, "sim_method") <- "wasserman sim"
  attr(f, "sim_pars") <- paste0("pi0:",pi0, ", xi_max:", xi_max)
  attr(f, "m") <- m
  f
}