#' t-test simulation: Simulate rowwise t-tests
#'
#' @param m Integer, total number of hypotheses 
#' @param pi0 Numeric, proportion of null hypotheses
#' @param effect_size Numeric, the alternative hypotheses will be 
#       distributed according to N(effect_size, 1), while nulls according to N(0,1).
#' @param n_samples Integer, number of samples for t-test, i.e.
#'      the comparison will be n_samples/2 vs n_samples/2
#' @param uninformative_filter Boolean, if TRUE will generate uniformly distributed filter statistic
#'     Otherwise will use the pooled standard deviations
#' @param seed Integer, Random seed to be used for simulation 
#'        (default: NULL, i.e. RNG state will be used as is)
#'
#' @return A data frame containing all information about the simulation experiment
#'
#' @examples sim_df <- du_ttest_sim(20000,0.95, 1.5)
#'
#' @export
du_ttest_sim <- function(m, pi0, effect_size, n_samples=10, uninformative_filter=FALSE, seed=NULL){
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

#' @describeIn du_ttest_sim Creates a closure function for a given seed
#' @export
du_ttest_sim_fun <- function(m, pi0, effect_size, n_samples=10, uninformative_filter=FALSE){
  f <- function(seed) du_ttest_sim(m, pi0, effect_size, n_samples=n_samples, uninformative_filter=uninformative_filter, seed=seed)
  attr(f, "sim_method") <- "t-test"
  attr(f, "sim_pars") <- paste0("pi0:",pi0, ", effect size:", effect_size)
  attr(f, "m") <- m
  f
}

#' Null simulation: Generate uniformly distributed p-values and covariates
#'
#' @param m Integer, total number of hypotheses
#' @param seed Integer, Random seed to be used for simulation 
#'        (default: NULL, i.e. RNG state will be used as is)
#'
#' @return A data frame containing all information about the simulation experiment
#'
#' @examples sim_df <- null_sim(20000)
#'
#' @export
null_sim <- function(m, seed=NULL){
  if (!is.null(seed)) set.seed(seed)

  sim_df <- data.frame(pvalue = runif(m), filterstat = runif(m), H=0)
  return(sim_df)
}

#' @describeIn null_sim Creates a closure function for a given seed
#' @export
null_sim_fun <- function(m){
  f <- function(seed) null_sim(m=m, seed=seed)
  attr(f, "sim_method") <- "only nulls"
  attr(f, "sim_pars") <- paste0("")
  attr(f, "m") <- m
  f
}

#' Normal simulation: Covariate is effect size under alternative
#'
#' @param m Integer, total number of hypotheses 
#' @param pi0 Numeric, proportion of null hypotheses
#' @param xi_min,xi_max  Numeric, covariates are drawn as uniform on xi_min, xi_max
#' @param seed Integer, Random seed to be used for simulation 
#'        (default: NULL, i.e. RNG state will be used as is)
#'
#' @return A data frame containing all information about the simulation experiment
#'
#' @examples sim_df <- wasserman_normal_sim(20000,0.9, 1, 5)
#'
#' @export
wasserman_normal_sim <- function(m, pi0, xi_min, xi_max, seed=NULL){
    if (!is.null(seed)) set.seed(seed)

    X   <- runif(m, min=xi_min, max=xi_max)
    H   <- rbinom(m,1,1-pi0)
    Z   <- rnorm(m, H*X)
    pvalue <- 1-pnorm(Z)
    simDf <- data.frame(pvalue=pvalue, filterstat=X,H=H, Z=Z)
}

#' @describeIn wasserman_normal_sim Creates a closure function for a given seed
#' @export
wasserman_normal_sim_fun <- function(m, pi0, xi_min, xi_max){
  f <- function(seed) wasserman_normal_sim(m, pi0, xi_min, xi_max, seed=seed)
  attr(f, "sim_method") <- "wasserman sim"
  attr(f, "sim_pars") <- paste0("pi0:",pi0, ", xi_max:", xi_max)
  attr(f, "m") <- m
  f
}



# Z_i correlated + PRDS
# Uninformative covariate


#' Normal PRDS simulation: Covariate is effect size under alternative, there are latent factors
#'      driving PRDS correlations among hypotheses
#'
#' @param m Integer, total number of hypotheses 
#' @param pi0 Numeric, proportion of null hypotheses
#' @param rho Numeric, correlation between z-scores of hypotheses driven by same latent factor
#' @param latent_factors Integer, number of latent factors driving the correlations
#' @param xi_min,xi_max  Numeric, covariates are drawn as uniform on xi_min, xi_max
#' @param seed Integer, Random seed to be used for simulation 
#'        (default: NULL, i.e. RNG state will be used as is)
#'
#' @return A data frame containing all information about the simulation experiment
#'
#' @examples sim_df <- wasserman_normal_prds_sim(20000,0.9, rho=0.1)
#'
#' @export
wasserman_normal_prds_sim <- function(m, pi0, rho=0.0, latent_factors=1, xi_min=0, xi_max=2.5, seed=NULL){
    if (!is.null(seed)) set.seed(seed)
    latent_idx <- sample(1:latent_factors, m, replace=TRUE)
    latent_Z <- rnorm(latent_factors)
    X   <- runif(m, min=xi_min, max=xi_max)
    H   <- rbinom(m,1,1-pi0)
    Z   <- sqrt(1-rho)*rnorm(m) + rho*latent_Z[latent_idx] + H*X
    pvalue <- 1-pnorm(Z)
    simDf <- data.frame(pvalue=pvalue, filterstat=X,H=H, Z=Z, latent_idx=latent_idx)
}

#' @describeIn wasserman_normal_prds_sim Creates a closure function for a given seed
#' @export
wasserman_normal_prds_sim_fun <- function(m, pi0, rho=0.0, latent_factors=1, xi_min=0, xi_max=2.5){
  f <- function(seed) wasserman_normal_prds_sim(m, pi0, rho=rho, latent_factors=latent_factors, 
                    xi_min=xi_min, xi_max=xi_max, seed=seed)
  attr(f, "sim_method") <- "wasserman normal prds sim"
  attr(f, "sim_pars") <- paste0("pi0:",pi0, ", xi_max:", xi_max, 
                                ", rho:", rho,
                                ", latent_factors:",latent_factors)
  attr(f, "m") <- m
  f
}

# TODO: Allow both covariate effects + dependence effects..
# pi0 not group dependent for now?
dependence_by_groups_sim <- function(m, pi0, rho=0.0, latent_factors=1,
                                     xi_min=0, xi_max=2.5, seed=NULL){
  if (!is.null(seed)) set.seed(seed)

  group_idx <- 1:latent_factors
  mgroups <- rep(floor(m/latent_factors), latent_factors)
  mgroups[1] <- mgroups[1] + m - sum(mgroups)

  xi_min_range <- xi_min + (0:(latent_factors-1))*(xi_max-xi_min)/latent_factors
  xi_max_range <- c(xi_min_range[-1], xi_max)

  sim_df <- rbind_all(mapply(wasserman_normal_prds_sim,
                             mgroups, pi0, rho, 1, xi_min_range,
                             xi_max_range,
                             SIMPLIFY=FALSE))
  sim_df$X <- sim_df$filterstat
  sim_df$filterstat <- unlist(mapply(function(idx, m) rep(idx,m), 
                                     group_idx, mgroups,
                                     SIMPLIFY=FALSE ))
  sim_df <- sim_df[sample(m,m),]
  sim_df
}

