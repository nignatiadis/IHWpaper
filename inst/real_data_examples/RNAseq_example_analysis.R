library("IHW")
library("IHWpaper")

bottomly <- analyze_dataset("bottomly")

ihw_res <- ihw(bottomly$pvalue, bottomly$baseMean, 0.1, nbins=13, 
               nfolds_internal=4L, nsplits_internal=5L)

breaks <- IHW:::stratification_breaks(ihw_res)

get_alpha_df <- function(alpha, pvalue, filter_statistic, filter_thresholds,...){
  print(paste0("alpha:",alpha))
  x <- ihw(pvalue, filter_statistic, alpha,...)
  print("IHW finished")
  ws <- weights(x, levels_only=TRUE)
  group_levels <- levels(groups_factor(x))
  folds  <- factor(1:x@nfolds)
  df <- expand.grid(stratum=1:nlevels(groups_factor(x)),
                    fold=folds)
  df$group <- group_levels[df$stratum]
  df$weight <- mapply(function(x,y) ws[x,y], df$stratum, df$fold)
  df$alpha <- alpha
  df$rejections <-rejections(x)
  df$bh_rejections <- sum(p.adjust(pvalue, method="BH") <= alpha, na.rm=T)
  for (filter_t in filter_thresholds){
    filt_pvalue <- pvalue[filter_statistic <= filter_t]
    df[paste0("threshold:", filter_t)] <- sum(p.adjust(filt_pvalue, method="BH") <= alpha, na.rm=T)
  }
  df
}

alpha_df <- bind_rows(lapply(seq(0.05,0.1,length=5), get_alpha_df, 
                             bottomly$pvalue, bottomly$baseMean, c(100,1000), nbins=13,
                             nsplits_internal=5))

res <- list(alpha_df = alpha_df,
            breaks   = breaks,
            break_min = min(bottomly$baseMean),
            ihw_res  = ihw_res
            )

saveRDS(res, file="result_files/RNAseq_benchmark.Rds")
