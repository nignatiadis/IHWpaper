## ----warning=F, message=F------------------------------------------------
library("ihwPaper")
library("fdrtool")
library("cowplot")

## ------------------------------------------------------------------------
sim <- du_ttest_sim(100, 0.5, 2.5 ,seed=100)

sorted_pvalues <- sort(sim$pvalue)
n  <- length(sorted_pvalues)
unique_pvalues <- unique(sorted_pvalues)
ecdf_values <- cumsum(tabulate(match(sorted_pvalues, unique_pvalues)))/n

df_ecdf <- data.frame(x=unique_pvalues,y=ecdf_values)
gren <- IHW:::presorted_grenander(sorted_pvalues)
df_gren <- data.frame(x=gren$x.knots, y=gren$y.knots)

## ------------------------------------------------------------------------
ggplot(df_ecdf, aes(x=x, y=y)) + geom_step(direction="hv",size=1.3) +
                  scale_x_continuous(expand=c(0,0),lim=c(0,1))+
                  scale_y_continuous(expand=c(0,0))+
                  xlab("p-value") +
                  ylab("Distribution")

## ----eval=FALSE----------------------------------------------------------
#  ggsave(filename="ecdf_plot.pdf", width=7,height=7)

## ------------------------------------------------------------------------
ggplot(df_ecdf, aes(x=x, y=y)) + geom_step(direction="hv",size=1.3) +
  geom_line(data=df_gren, aes(x=x,y=y), color="red",size=1.3) +
                  scale_x_continuous(expand=c(0,0),lim=c(0,1))+
                  scale_y_continuous(expand=c(0,0))+
                  xlab("p-value") +
                  ylab("Distribution")

## ----eval=FALSE----------------------------------------------------------
#  ggsave(filename="ecdf_grenander_plot.pdf", width=7,height=7)

