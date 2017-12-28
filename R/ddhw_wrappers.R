#' IHW wrappers
#' @param unadj_p  Numeric vector of unadjusted p-values.
#' @param filterstat   Factor to which different hypotheses belong
#' @param alpha    Significance level at which to apply method
#'
#' @return ihwResult multiple testing object
#'
#' @examples
#'      sim_df <- du_ttest_sim(20000,0.95, 1.5)
#'      obj <- ihw_5fold(sim_df$pvalue, sim_df$filterstat, .1)
#'      sum(rejected_hypotheses(obj))
#'
#' @details
#'     These are closures, which apply IHW with custom
#'   prespecified parameters. These correspond to interesting
#'   settings, for which it is convenient to be able to immediately
#'   call the corresponding functions, rather than having to specify
#'   parameters each time. Thus they make it easier to benchmark.
#'   All of these wrappers are defined in 2 lines of code, 
#'   so the settings pertaining
#'   to each one can be inspected by typing the functions name into the
#'   console.
#'
#'
#' @import IHW
#' @export
#' @describeIn ihw_naive IHW naive
ihw_naive <- function(unadj_p, filterstat, alpha){
 	obj <- ihw(unadj_p, filterstat, alpha, nbins=20, nfolds=1, lambdas=Inf,
 		distrib_estimator="ECDF", lp_solver="gurobi")
 	obj
}

attr(ihw_naive, "testing covariate") <- "continuous"
attr(ihw_naive, "fdr_method")        <- "IHW naive"


#' @describeIn ihw_naive  IHW (E2) with 5 folds
#' @export
ihw_ecdf_5fold <- function(unadj_p, filterstat, alpha){
	obj <- ihw(unadj_p, filterstat, alpha, nbins=20, nfolds=5, lambdas=Inf,
 		distrib_estimator="ECDF", lp_solver="gurobi")
 	obj
}

attr(ihw_ecdf_5fold, "testing covariate") <- "continuous"
attr(ihw_ecdf_5fold, "fdr_method")        <- "IHW ECDF"


#' @describeIn ihw_naive  IHW (E1-E2) with 5 folds
#' @export
ihw_5fold <- function(unadj_p, filterstat, alpha){
 	obj <- ihw(unadj_p, filterstat, alpha, nbins=20, nfolds=5, lambdas=Inf,
 		distrib_estimator="grenander", lp_solver="lpsymphony")
 	obj
}

attr(ihw_5fold, "testing covariate") <- "continuous"
attr(ihw_5fold, "fdr_method")        <- "IHW"

#' @describeIn ihw_naive  IHW (E1-E2-E3) with 5 folds
#' @export
ihw_5fold_reg <- function(unadj_p, filterstat, alpha){
  obj <- ihw(unadj_p, filterstat, alpha, nbins=20, nfolds=5,
             distrib_estimator="grenander", lp_solver="lpsymphony")
  obj
}

attr(ihw_5fold_reg, "testing covariate") <- "continuous"
attr(ihw_5fold_reg, "fdr_method")        <- "IHW E3"



#--------------- IHW-Bonferroni----------------------------------------

#' @describeIn ihw_naive  IHW-Bonferroni (E1-E2-E3) with 5 folds
#' @export
ihw_bonf_5fold_reg <- function(unadj_p, filterstat, alpha){
  obj <- ihw(unadj_p, filterstat, alpha, nbins=20, nfolds=5,
             distrib_estimator="grenander", lp_solver="lpsymphony",
             adjustment_type="bonferroni")
  obj
}

attr(ihw_bonf_5fold_reg, "testing covariate") <- "continuous"
attr(ihw_bonf_5fold_reg, "fdr_method")        <- "IHW-Bonferroni E3"


#------------------ IHW-Storey ------------------------------------------
#' @describeIn ihw_naive  IHW (E1-E2) with 5 folds and Storey's pi0 estimator
#' @export
ihw_storey_5fold <- function(unadj_p, filterstat, alpha){
    obj <- ihw(unadj_p, filterstat, alpha, nbins=20, nfolds=5, lambdas=Inf,
        distrib_estimator="grenander", lp_solver="lpsymphony",
        null_proportion=TRUE, null_proportion_level = 0.5)
    obj
}

attr(ihw_storey_5fold, "testing covariate") <- "continuous"
attr(ihw_storey_5fold, "fdr_method")        <- "IHW-Storey"


