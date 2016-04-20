## ----warning=F, message=F------------------------------------------------
library("IHW")
library("IHWpaper")
library("ggplot2")
library("dplyr")
library("wesanderson")

## ---- warning=F----------------------------------------------------------
sim <- du_ttest_sim(80000, 0.9, 2, seed=1)
ddhf_res <- ddhf(sim$pvalue, sim$filterstat, .1)
ws <- ifelse(sim$filterst <= ddhf_res$cutoff_value, 0, 1)
ws <- ws/sum(ws)*length(sim$pvalue)

## ------------------------------------------------------------------------
total_regularization_lambda <- max(ws)-min(ws)
total_regularization_lambda

## ----warning=F-----------------------------------------------------------
ihw_res <- ihw(sim$pvalue, sim$filterstat, .1, lambdas = total_regularization_lambda, nfolds=1L,  nbins=20)

## ----fig.width=5, fig.height=5, warning=F--------------------------------
df <- rbind( data.frame(covariate = sim$filterstat, weight= ws, method="filtering"),
             data.frame(covariate = sim$filterstat, weight= weights(ihw_res, levels_only=FALSE), 
                        method=paste0("IHW; \nlambda=",format(total_regularization_lambda,digits=2))))
weights_filter_plot <- ggplot(df, aes(x=covariate, y=weight, col=method)) + geom_step(size=1.65)+
                      scale_colour_manual(values=wes_palette("Cavalcanti")[c(1,2)]) +
                      theme_classic(16)
weights_filter_plot

## ---- eval=FALSE---------------------------------------------------------
#  pdf("smoothed_threshold.pdf", width=5, height=5)
#  weights_filter_plot
#  dev.off()

## ----message=F, warning=F------------------------------------------------
ihw_res2 <- ihw(sim$pvalue, sim$filterstat, .1, lambdas=seq(0,10,length=20), nfolds=1L, nfolds_internal = 5L, nbins=20, quiet=TRUE)

## ----message=F, warning=F------------------------------------------------
regularization_term(ihw_res2)

## ----fig.width=5, fig.height=5, warning=F--------------------------------
df <- rbind(df,
            data.frame(covariate = sim$filterstat, weight= weights(ihw_res2, levels_only=FALSE),  
                       method=paste0("IHW; \nlambda=",
                                     format(regularization_term(ihw_res2) ,digits=2))))
            
weights_filter_plot <- ggplot(df, aes(x=covariate, y=weight, col=method)) + geom_step(size=1.65)+
                      scale_colour_manual(values=wes_palette("Cavalcanti")[c(1,2,3)]) +
                      theme_classic(16)

weights_filter_plot

## ---- warning=F----------------------------------------------------------
group_by(df, method) %>% summarize(mean_weight = mean(weight))

