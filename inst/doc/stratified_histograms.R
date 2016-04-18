## ----warning=FALSE, message=FALSE----------------------------------------
library("ggplot2")
library("cowplot")

## ------------------------------------------------------------------------
set.seed(1)

m <- 10000
binwidth <- 0.025

# generic function to add some properties to gg histograms
gg_hist_aesthetic <-  function(gg_obj, ylim_max=650) {
    gg_obj + 
             scale_x_continuous(expand = c(0.02, 0)) + 
             scale_y_continuous(expand = c(0.02, 0), limits=c(0,ylim_max)) + 
             xlab("p-value")+
             ylab("Counts")+
             theme(axis.title = element_text(face="bold",size=rel(0.1)) )
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
              geom_hline(yintercept=9000*binwidth, size=1.3, col="darkblue")
gg_beta_1 <- gg_hist_aesthetic(gg_beta_1)
gg_beta_1

## ----eval=FALSE----------------------------------------------------------
#  ggsave(plot=gg_beta_1, file="stratified_histograms_beta1.pdf", width=4, height=3)

## ------------------------------------------------------------------------
pv_beta_2 <- data.frame(pvalue=c(runif(5500), rbeta(4500,1,4)))

gg_beta_2 <- ggplot(pv_beta_2, aes(x=pvalue)) + 
              geom_histogram(binwidth = binwidth, boundary = 0, colour="lightgrey", fill="#939598") +
              geom_hline(yintercept=5500*binwidth, size=1.3, col="darkblue")
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

## ------------------------------------------------------------------------
set.seed(1)

m <- 10000
binwidth <- 0.05

# generic function to add some properties to gg histograms
gg_hist_aesthetic <-  function(gg_obj, ylim_max=6) {
    gg_obj + 
             aes(y=..density..)+
             scale_x_continuous(expand = c(0.02, 0), breaks=c(0,0.5,1)) + 
             scale_y_continuous(expand = c(0.02, 0), limits=c(0,ylim_max)) + 
             xlab("p-value")+
             ylab("Frequency")+
             theme(axis.title = element_text(face="bold",size=rel(0.7)))
}

## ------------------------------------------------------------------------
pv_unif <- data.frame(pvalue=runif(m))
gg_unif <- ggplot(pv_unif, aes(x=pvalue)) + 
            geom_histogram(binwidth = binwidth, boundary = 0, colour="lightgrey", fill="#939598") +
            geom_hline(yintercept=1, size=1.3, col="darkblue")
gg_unif <- gg_hist_aesthetic(gg_unif)
 

gg_unif

## ------------------------------------------------------------------------
pv_beta_a <- data.frame(pvalue=c(runif(9000), rbeta(1000,0.5,4)))

gg_beta_a <- ggplot(pv_beta_a, aes(x=pvalue)) + 
              geom_histogram(binwidth = binwidth, boundary = 0, colour="lightgrey", fill="#939598") +
              geom_hline(yintercept=0.9, size=1.3, col="darkblue")
gg_beta_a <- gg_hist_aesthetic(gg_beta_a)
gg_beta_a

## ------------------------------------------------------------------------
pv_beta_b <- data.frame(pvalue=c(runif(8000), rbeta(2000,0.5,11)))

gg_beta_b <- ggplot(pv_beta_b, aes(x=pvalue)) + 
              geom_histogram(binwidth = binwidth, boundary = 0, colour="lightgrey", fill="#939598") +
              geom_hline(yintercept=0.8, size=1.3, col="darkblue")
gg_beta_b <- gg_hist_aesthetic(gg_beta_b, ylim_max=5.5)
gg_beta_b

## ------------------------------------------------------------------------
pv_all <- rbind(pv_unif, pv_beta_b, pv_beta_a) 
gg_all <- ggplot(pv_all, aes(x=pvalue)) + 
              geom_histogram(binwidth = binwidth, boundary = 0, colour="lightgrey", fill="#939598") +
              geom_hline(yintercept=0.9, size=1.3, col="darkblue")
gg_all <- gg_hist_aesthetic(gg_all)
gg_all

## ----fig.width=6, fig.height=2-------------------------------------------
gg_stratified <- plot_grid(gg_all, gg_beta_b, gg_beta_a, gg_unif,
                         nrow=1,
                         labels=c("a)", "b)", "c)","d)"))
gg_stratified

## ----eval=FALSE----------------------------------------------------------
#  ggsave(plot=gg_stratified, file="stratified_histograms.pdf", width=6, height=2)

