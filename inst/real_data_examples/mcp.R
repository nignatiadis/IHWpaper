tmp <- read.csv("mcp.csv", stringsAsFactors = F)
pv <- as.numeric(tmp$p.value.10min.0min.rep1)
cov <- tmp$number.of.ratios.10min.0min.rep1

sum(p.adjust(pv, method="BH") <= .05, na.rm=T)
tmp <- ddhw(pv,cov,.05, nbins=3)

tmp <- read.csv("mcp_phospho.csv", stringsAsFactors=FALSE)
pv2 <- as.numeric(tmp$p.value.10min.0min.rep1)
cov2 <- as.factor(tmp$numPhosphoSites)
sum(p.adjust(pv2,method="BH") <=.1 ,na.rm=T)

ddhw_res <- ddhw(pv2,cov2,.1)


nmeth <- read.csv("nmeth.2518-S2.csv", stringsAsFactors=FALSE)
#nmeth <- mutate(nmeth, count=Ratio.H.L.count.Rep1 + Ratio.L.H.count.Rep2 + Ratio.H.L.count.Rep3)
nmeth <- mutate(nmeth, count=Ratio.H.L.count.Rep1 )

adj_pv <- nmeth$adj.P.Val
non_na <- !is.na(adj_pv)
adj_pv <- adj_pv[non_na]
order_pv <- order(adj_pv)
sorted_adj_pv <- adj_pv[order_pv]
sorted_pv <- pmin(sorted_adj_pv /length(sorted_adj_pv)*(1:length(sorted_adj_pv)),1)

cov <- nmeth$count[non_na][order_pv]
sum(adj_pv <= 0.1, na.rm=T)

sum(p.adjust(sorted_pv, method="BH")<=.1, na.rm=T)

ddhw_res <- ddhw(sorted_pv,cov,.1)

ggplot(nmeth, aes(x=rank(Proteins), y=log(adj.P.Val))) + geom_point() 

library("NBPSeq")
data(arab)

countData <- arab
colData <- data.frame(condition= gsub("\\d", "",colnames(arab)))
rownames(colData) <- colnames(countData)

dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
sum(p.adjust(res$pvalue, method="BH")<.05,na.rm=T)

tmp <- ddhw(res$pvalue,res$baseMean,.05, lambdas = seq(0.1,5, length=20))


## maxquant
tmp <- read.csv("maxquant.csv", stringsAsFactors = F)
pvalue <- tmp$Ratio.H.L.Significance.B.
pvalue2 <- tmp$Ratio.H.L.Significance.A.
intensity <- tmp$Intensity
sum(p.adjust(pvalue2, method="BH") <= .00005 & tmp$Ratio.H.L.Count>=3, na.rm=T)

sum(p.adjust(pvalue, method="BH") <= 0.1,na.rm=T)
tmp_ddhw <- ddhw(pvalue, intensity, .1, nbins=4)

tmp_ddhw3 <- ddhw(pvalue, tmp$Ratio.H.L.Count, .1, nbins=4)


sum(p.adjust(pvalue2, method="BH") <= 0.05,na.rm=T)
tmp_ddhw2 <- ddhw(pvalue2, intensity, .05, nbins=4)
plot_ddhw(tmp_ddhw2)

#tmp2 <- arrange(tmp, Ratio.H.L.Significance.B.)
tmp2 <- filter(tmp, Contaminant != "+", Reverse != "+", !is.na(Ratio.H.L.Normalized))
tmp2$group <- DDHW:::groups_by_filter(tmp2$Intensity,3)  
tmp2 <- group_by(tmp2,group) %>% mutate(pvalue=locfdr_wrapper(log(Ratio.H.L.Normalized)))
sum(p.adjust(tmp2$pvalue, method="BH") <= .1)
sum(p.adjust(tmp2$Ratio.H.L.Significance.B., method="BH") <= .1)

a <- ddhw(tmp2$pvalue, tmp2$Intensity, .1, nbins=3)
rejections(a)


locfdr_wrapper <- function(x){
  res <- locfdr(x, pct=0.02)
  mu0 <- res$fp0[5,1]
  sd0 <- res$fp0[5,2]
  pv <- 2*pmin(1 - pnorm(x, mu0, sd0),
               pnorm(x, mu0, sd0))
 # pv <- 1-pnorm(x,mu0,sd0)
  pv
}
fdrtool_wrapper <- function(x){
  df<-fdrtool(x, statistic="normal", plot=F)
  df$pval
}
fdrtool_wrapper_pi0 <- function(x){
  df<-fdrtool(x, statistic="normal", plot=F)
  df$param[3]
}

fdrreg_wrapper <- function(x){
    e1 = efron(x, nulltype='empirical', df=7)
    print(e1$mu0)
    print(e1$sig0)
    pv <- 2*pmin(1 - pnorm(x, e1$mu0, e1$sig0),
            pnorm(x, e1$mu0, e1$sig0))
    pv
}

fdrreg_wrapper <- function(x){
  e1 = efron(x, nulltype='empirical', df=7)
  print(e1$mu0)
  print(e1$sig0)
  #pv <- 2*pmin(1 - pnorm(x, e1$mu0, e1$sig0),
  #             pnorm(x, e1$mu0, e1$sig0))
  pv <- 1-pnorm(x,e1$mu, e1$sig0)
  pv
}

fdrreg_wrapper_pi0 <- function(x){
  e1 = efron(x, nulltype='empirical', df=7)
  print(e1$mu0)
  print(e1$sig0)
  #pv <- 2*pmin(1 - pnorm(x, e1$mu0, e1$sig0),
  #             pnorm(x, e1$mu0, e1$sig0))
 #data.frame(pi0=e1$p0, mu0=e1$mu0, sig0=e1$sig0)
  e1$mu0
}


group_by(tmp2,group) %>% summarize(pi0=fdrtool_wrapper_pi0(log(Ratio.H.L.Normalized)))

sum(p.adjust(tmp2$Ratio.H.L.Significance.A., method="BH")<=0.1, na.rm=T)
tmp_ddhw4 <- ddhw(tmp2$Ratio.H.L.Significance.A., tmp2$Intensity, .1, nbins=3)
rejections(tmp_ddhw4)
plot_ddhw(tmp_ddhw4)


###############

tmp <- read.csv("hsp.csv", stringsAsFactors = F)
qv <- qvalue(tmp$p.value)
get_bh_threshold(tmp$p.value, .1)
get_bh_threshold(tmp$p.value, .1/qv$pi0)

hist(tmp$p.value)
sum(p.adjust(tmp$p.value, method="BH")<=.0005)
sum(qv$qvalues <= .0005)
tst <- ddhw(tmp$p.value, tmp$Intensity, .0005, nbins=4, lambdas=seq(0,3, length=20))
plot_ddhw(tst) + scale_x_discrete(breaks=c(1,2,3,4)) + xlab("stratum by increasing intensity")
ggsave("HSP_weights.pdf", width=7,height=5)
tst_adaptive <- ddhw(tmp$p.value, tmp$Intensity, .0005/qv$pi0, nbins=4, lambdas=seq(0,3, length=20))
rejections(tst_adaptive)

rejections(tst)
ggplot(tst@df,aes(x=pvalue)) + geom_histogram() + facet_grid(.~group)

hist(tmp$p.value)


### SCIENCE signaling
tmp <- read.csv("science_signaling.csv", stringsAsFactors = F)
pvalue <- 2*(1-pt(abs(tmp$Avg.log2.ratio)/tmp$SD*sqrt(6), 5)) #maybe fix this later
pvalue <- rank(tmp$p1, ties.method="first")*tmp$p1/nrow(tmp) #use this for now
tmp$pvalue <- pvalue
sum(p.adjust(pvalue, method="BH") <= .05) #246
ddhw_res <- ddhw(tmp$pvalue,tmp$X..peptides, .1, nbins=5, 
                 nsplits_internal=50, lambdas=seq(0,3,length=20))
#ddhw_res2 <- ddhw(tmp$pvalue,tmp$X..peptides, .1, nbins=4, nsplits_internal=10,
 #                 lambdas=seq(0,2.5,length=20))

get_alpha_df <- function(alpha, pvalue, filter_statistic, filter_thresholds,...){
  print(paste0("alpha:",alpha))
  x <- ddhw(pvalue, filter_statistic, alpha,...)
  print("DDHW finished")
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


mydfs2 <- rbind_all(lapply(seq(0.05,0.1,length=5), get_alpha_df, tmp$pvalue, tmp$X..peptides, c(), nbins=5, 
                           nsplits_internal=10, lambdas=seq(0,3,length=20)))
breaks <- DDHW:::stratification_breaks(ddhw_res)
breaks <- c(0,breaks[-length(breaks)])
mydfs2 <- mutate(mydfs2, peptideCount = breaks[stratum])
summ2 <- group_by(mydfs2, alpha) %>% summarize(rjs_bh = max(bh_rejections), rjs_ddhw=max(rejections)) %>% 
  gather(method, rejections, rjs_bh, rjs_ddhw)
  
ggplot(summ2, aes(x=alpha,y=rejections,col=method)) +  geom_line() 
ggsave("gygi_rejections.pdf", width=7,height=5)

ggplot(filter(mydfs2, alpha==0.05), aes(x=peptideCount,y=weight, col=fold)) +
  geom_point()+ 
  geom_line()+
  xlab("peptideCount (left break)")
ggsave("gygi_weights.pdf", width=7,height=5)
