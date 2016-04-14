## ----warning=FALSE, message=FALSE----------------------------------------
library("ggplot2")
library("cowplot")

## ------------------------------------------------------------------------
set.seed(1)

m <- 10000
binwidth <- 0.025

# generic function to add some properties to gg histograms
gg_hist_aesthetic <-  function(gg_obj) {
    gg_obj + scale_x_continuous(expand = c(0.02, 0)) + 
             scale_y_continuous(expand = c(0.02, 0), limits=c(0,650)) + 
             xlab("p-value")+
             theme(axis.title = element_text(face="bold") )
}

## ------------------------------------------------------------------------
pv_unif <- data.frame(pvalue=runif(m))
gg_unif <- ggplot(pv_unif, aes(x=pvalue)) + 
            geom_histogram(binwidth = binwidth, boundary = 0, colour="lightgrey", fill="#939598") +
            geom_hline(yintercept=m*binwidth, size=2, col="darkblue")
gg_unif <- gg_hist_aesthetic(gg_unif)
 

gg_unif

## ----eval=FALSE----------------------------------------------------------
#  ggsave(plot=gg_unif, file="stratified_histograms_unif.pdf", width=4, height=3)

## ------------------------------------------------------------------------
pv_beta_1 <- data.frame(pvalue=c(runif(9000), rbeta(1000,0.5,7)))

gg_beta_1 <- ggplot(pv_beta_1, aes(x=pvalue)) + 
              geom_histogram(binwidth = binwidth, boundary = 0, colour="lightgrey", fill="#939598") +
              geom_hline(yintercept=9000*binwidth, size=2, col="darkblue")
gg_beta_1 <- gg_hist_aesthetic(gg_beta_1)
gg_beta_1

## ----eval=FALSE----------------------------------------------------------
#  ggsave(plot=gg_beta_1, file="stratified_histograms_beta1.pdf", width=4, height=3)

## ------------------------------------------------------------------------
pv_beta_2 <- data.frame(pvalue=c(runif(5500), rbeta(4500,1,4)))

gg_beta_2 <- ggplot(pv_beta_2, aes(x=pvalue)) + 
              geom_histogram(binwidth = binwidth, boundary = 0, colour="lightgrey", fill="#939598") +
              geom_hline(yintercept=5500*binwidth, size=2, col="darkblue")
gg_beta_2 <- gg_hist_aesthetic(gg_beta_2)
gg_beta_2

## ----eval=FALSE----------------------------------------------------------
#  ggsave(plot=gg_beta_2, file="stratified_histograms_beta2.pdf", width=4, height=3)

## ----fig.width=12, fig.height=4------------------------------------------
gg_stratified <- plot_grid(gg_unif, gg_beta_1, gg_beta_2,
                         nrow=1,
                         labels=c("a)", "b)", "c)"))
gg_stratified

## ----eval=FALSE----------------------------------------------------------
#  ggsave(plot=gg_stratified, file="stratified_histograms.pdf", width=12, height=4)

