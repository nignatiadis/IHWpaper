library("IHWpaper")
library("IHW")

proteomics_file <- system.file("extdata/real_data",
                                "science_signaling.csv", package = "IHWpaper")

proteomics_df <- read.csv(proteomics_file, stringsAsFactors = F)

proteomics_df$pvalue <- rank(proteomics_df$p1, ties.method="first")*proteomics_df$p1/nrow(proteomics_df) 

ihw_res <- ihw(proteomics_df$pvalue,proteomics_df$X..peptides, .1, nbins=4, 
                 nsplits_internal=5, lambdas=seq(0,3,length=20))


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


alpha_df <- bind_rows(lapply(seq(0.05, 0.1, length = 5), get_alpha_df, proteomics_df$pvalue, 
                            proteomics_df$X..peptides,
                            c(), nbins=4, 
                            nsplits_internal=5, 
                            lambdas=seq(0,3,length=20)))

breaks <- IHW:::stratification_breaks(ihw_res)
break_min <- min(proteomics_df$X..peptides)

res <- list(alpha_df = alpha_df,
            breaks   = breaks,
            break_min = break_min,
            ihw_res  = ihw_res
            )

saveRDS(res, file="result_files/proteomics_benchmark.Rds")
