# here we gather again all methods which have been defined 
# in R/covariate_methods.R, R/stratified_methods.R , R/simple_methods.R
# so that they can be easily used in the benchmarking scripts

# In particular they are modified so as to take "simulations" as inputs.
# These are dataframe so that their p-values can be extracted by 
# sim$pvalue and the covariate by sim$filterstat


continuous_wrap <- function(mt_method, nbins=20){
	print(attr(mt_method,"fdr_method"))
	if (attr(mt_method, "testing covariate") == "simple"){

		f <- function(sim, alpha){
			mt_method(sim$pvalue, alpha)
		}
		attr(f, "fdr_method") <- attr(mt_method, "fdr_method")

	} else if (attr(mt_method, "testing covariate") == "stratified"){

		f <- function(sim, alpha){
			groups <- ddhw::groups_by_filter(sim$filterstat, nbins)
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

continuous_methods_list <- list(bh,
							    storey_qvalue,
							    lsl_gbh,
							    tst_gbh,
							    stratified_bh,
							    cai,
							    scott_fdrreg,#,
							    ddhf,
							    ddhw_lp_adaptive,
							    ddhw_lp,
							    ddhw_lp_adaptive_nosd,
							    ddhw_lp_unregularized)

wrapped_continuous_methods_list <- lapply(continuous_methods_list, continuous_wrap)