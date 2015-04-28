
ddhf <- function(unadj_p, filterstat, alpha){
  sortedp <- sort(unadj_p)
  ranksp <- rank(unadj_p, ties.method="first")
  filterorder <- order(filterstat)
  filterrank <- rank(filterstat)

  sortedp_filterorder <- ranksp[filterorder]
  rjs <- ddhf_order(sortedp, sortedp_filterorder, alpha)
  j <- which.max(rjs)[1]

  m <- length(unadj_p)
  adj_p <- rep(1, m)
  adj_p[filterrank >= j] <- p.adjust(unadj_p[filterrank >= j], method="BH")

  ddhf_obj <- list(adj_pvalues = adj_p, rejections = max(rjs),
                   cutoff_quantile=j/m,
                   cutoff_value = filterstat[filterorder][j])
  class(ddhf_obj) <- "ddhf"
  ddhf_obj
}

 #pv <- runif(1000)
 #ord <- sample(1:1000,1000)
 # tmp <- ddhf_order(pv,ord, .1)
