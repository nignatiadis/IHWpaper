## ----warning=F, message=F------------------------------------------------
library("ggplot2")
library("dplyr")
library("tidyr")
library("gridExtra")
library("cowplot")
library("wesanderson")

## ----warning=F-----------------------------------------------------------
falt1 <- function(t){
    1-pnorm(qnorm(1-t)-mu1)
  }

  lfdr1 <- function(t){
    pi10/(pi10 + (1-pi10)*falt1_density(t))
  }
  falt1_density <- function(t){
      dnorm(qnorm(1-t)-mu1)/dnorm(qnorm(1-t))
  }

#http://stackoverflow.com/questions/10081479/solving-for-the-inverse-of-a-function-in-r
inverse = function (f, lower = -100, upper = 100) {
   function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]$roo
}

sim_pars <- function(eff_size, pi0){
  distrib <- function(t) pi0*t + (1-pi0)*(1-pnorm(qnorm(1-t)-eff_size))
  alt_density <-   function(t) dnorm(qnorm(1-t)-eff_size)/dnorm(qnorm(1-t))
  density <- function(t) pi0 + (1-pi0)*alt_density(t)
  tdr <- function(t)  1- pi0/density(t)
  inverse_f <- function(y) uniroot((function (x) density(x) - y), lower = 0.001, upper = 0.999)[1]$root
  tdr_density <- function(t) pi0*(inverse_f(pi0/(1-t)+0.001) - inverse_f(pi0/(1-t)-0.001))/0.002
  return(list(pi0=pi0,distrib=distrib, alt_density=alt_density,density=density, tdr=tdr, inverse_f=inverse_f,
            tdr_density=tdr_density  ))
}


## ----warning=F-----------------------------------------------------------
plot_tdr <- function(x,pi0,xthreshold=0.05){
  eff_size <- x
  pars <- sim_pars(x,pi0)
  df <- data.frame(t=seq(0.001,1, length=1000))
  df <- mutate(df, f1 = pars$alt_density(t), f=pars$density(t), tdr = pars$tdr(t))
  
  plot_pv <- ggplot(df, aes(x=t, y=f)) + 
                     geom_line()  + 
                     expand_limits(x = 0)  +
                     geom_vline(xintercept=xthreshold, col="blue") + 
                     geom_hline(yintercept=pars$pi0,col="red") +
                     #geom_hline(yintercept=pars$density(xthreshold), linetype=2)+
                     ylab("Density")+
                     xlab("p-value")+
                     scale_y_continuous(limits=c(0,4))+
                     scale_x_continuous(expand=c(0,0),breaks=c(0, 0.25,0.5, 0.75,1),
                                        labels= c("0","0.25","0.5","0.75","1"))   
  
  plot_pv_vs_tdr <- ggplot(df, aes(x=t, y= tdr)) + 
                    geom_line() +
                    ylab("tdr")+
                    xlab("p-value")+
                    expand_limits(x = 0)  +
                    scale_x_continuous(expand=c(0,0),breaks=c(0, 0.25,0.5, 0.75,1),
                                        labels= c("0","0.25","0.5","0.75","1"))   
  # now also generate empirical histogram plots
  a = sim_pars(eff_size,pi0)
  m=100000
  H <- rbinom(m, 1, 1-pi0)
  X <- rnorm(m) + eff_size*H
  PV <- 1 - pnorm(X)
  TDR <- a$tdr(PV)
  
  sim_df <- data.frame(pvalue=PV, tdr=TDR)
  
  hist_pv <- ggplot(sim_df, aes(x=pvalue))
  grid.arrange(plot_pv, plot_pv_vs_tdr, ncol=1)
  
  pv_hist <- ggplot(sim_df, aes(x=pvalue)) + 
        geom_histogram(aes(y = ..density..),binwidth=0.04, 
                       boundary=1.0,colour = "darkgreen", fill = "white") +
        expand_limits(x = 0)  +
        scale_x_continuous(expand=c(0,0),breaks=c(0, 0.25,0.5, 0.75,1),
                                        labels= c("0","0.25","0.5","0.75","1"))   +
        xlab("p-value") + 
        ylab("Density")
  
  tdr_hist <- ggplot(sim_df, aes(x=tdr)) +
                geom_histogram(aes(y = ..density..),
                              binwidth=0.04,boundary=1.0, colour = "darkgreen", fill = "white") +
              expand_limits(x = 0)  +
              scale_x_continuous(expand=c(0,0),breaks=c(0, 0.25,0.5, 0.75,1),
                                        labels= c("0","0.25","0.5","0.75","1"))   +
              xlab("tdr") + 
              ylab("Density")
  list(plot_pv=plot_pv,  plot_pv_vs_tdr=plot_pv_vs_tdr, pv_hist=pv_hist, tdr_hist=tdr_hist)
}

## ----fig.keep="none", warning=F------------------------------------------
myplots <- c(plot_tdr(3.5, 0.6), plot_tdr(3.5,0.85), plot_tdr(1.5,0.6))
#myplots <- myplots[c(1,5,9) + rep(0:3,each=3)]
letter_labels <- letters[c(1,4,7,10) + rep(0:2,each=4)]
myplots <- c(myplots, list(labels = paste0(letter_labels,")"), nrow=3, align="hv"))

## ----warning=F, fig.width=10, fig.height=6-------------------------------
supp_fig <- do.call(plot_grid,myplots)
supp_fig

## ----eval=F--------------------------------------------------------------
#  ggsave(supp_fig, filename="pvalue_vs_tdr_explanation.pdf", width=10, height=6)

## ----warning=F-----------------------------------------------------------

df <- data.frame(t=seq(0.00001,1, length=1000))
x <- 2.5
pi0 <- 0.8
eff_size <- x
pars <- sim_pars(x,pi0)
df1 <- mutate(df, f1 = pars$alt_density(t), CDF=pars$distrib(t), f=pars$density(t), tdr = pars$tdr(t), test = "1")
  
x <- 1.5
pi0 <- 0.9
eff_size <- x
pars <- sim_pars(x,pi0)
df2 <- mutate(df, f1 = pars$alt_density(t), CDF=pars$distrib(t), f=pars$density(t), tdr = pars$tdr(t), test = "2")

df <- rbind(df1,df2)

plot_tdr_inversed1 <- ggplot(df, aes(x=tdr, y= t,color=test)) + 
                    geom_line() +
                    xlab("tdr")+
                    ylab("p-value")+
                    expand_limits(x = c(0,1), y=0)  +
                    scale_x_continuous(expand=c(0,0),breaks=c(0, 0.25,0.5, 0.75,1),
                                        labels= c("0","0.25","0.5","0.75","1"))   +  
                     scale_y_continuous(expand=c(0,0),breaks=c(0, 0.25,0.5, 0.75,1),
                                        labels= c("0","0.25","0.5","0.75","1"))   +  
                    scale_colour_manual(values =wes_palette("Chevalier")[c(1,4)], name='',
                                        labels =  expression(H[1] , H[2]))


plot_density1 <- ggplot(df, aes(x=t, y= f,color=test)) + 
                    geom_line() +
                    xlab("p-value")+
                    ylab("density")+
                    expand_limits(x = 0,y=0)  +
                    scale_x_continuous(expand=c(0,0),breaks=c(0, 0.25,0.5, 0.75,1),
                                        labels= c("0","0.25","0.5","0.75","1"))   +  
                    scale_y_continuous(expand=c(0,0),breaks=c(0, 0.25,0.5, 0.75,1),
                                       labels= c("0","0.25","0.5","0.75","1"), 
                                       limits=c(0,4))   +  
                    scale_colour_manual(values =wes_palette("Chevalier")[c(1,4)], name='',
                                         expression(H[1] , H[2]))


plot_cdf1 <- ggplot(df, aes(x=t, y= CDF,color=test)) + 
                    geom_line() +
                    xlab("p-value")+
                    ylab("CDF")+
                    expand_limits(x = 0,y=0)  +
                    scale_x_continuous(expand=c(0,0),breaks=c(0, 0.25,0.5, 0.75,1),
                                        labels= c("0","0.25","0.5","0.75","1"))   +  
                    scale_y_continuous(expand=c(0,0),breaks=c(0, 0.25,0.5, 0.75,1),
                                       labels= c("0","0.25","0.5","0.75","1"))   +  
                    scale_colour_manual(values =wes_palette("Chevalier")[c(1,4)], name='',
                                        labels =  expression(H[1] , H[2]))


## ----warning=F-----------------------------------------------------------
df <- data.frame(t=seq(0.00001,1, length=1000))
x <- 2.5
pi0 <- 0.8
eff_size <- x
pars <- sim_pars(x,pi0)
df1 <- mutate(df, f1 = pars$alt_density(t), CDF=pars$distrib(t), f=pars$density(t), tdr = pars$tdr(t), test = "1")
  
x <- 2.5
pi0 <- 0.9
eff_size <- x
pars <- sim_pars(x,pi0)
df2 <- mutate(df, f1 = pars$alt_density(t), CDF=pars$distrib(t), f=pars$density(t), tdr = pars$tdr(t), test = "2")
df <- rbind(df1,df2)
  
plot_tdr_inversed <- ggplot(df, aes(x=tdr, y= t,color=test)) + 
                    geom_line() +
                    xlab("tdr")+
                    ylab("p-value")+
                    expand_limits(x = c(0,1),y=0)  +
                    scale_x_continuous(expand=c(0,0),breaks=c(0, 0.25,0.5, 0.75,1),
                                        labels= c("0","0.25","0.5","0.75","1"))   +  
                    scale_y_continuous(expand=c(0,0),breaks=c(0, 0.25,0.5, 0.75,1),
                                       labels= c("0","0.25","0.5","0.75","1"))   +  
                    scale_colour_manual(values =wes_palette("Chevalier")[c(1,4)], name='',
                                          labels = expression(H[1] , H[2]))



plot_density <- ggplot(df, aes(x=t, y= f1,color=test)) + 
                    geom_line() +
                    xlab("p-value")+
                    ylab("alt")+
                    expand_limits(x = 0,y=0)  +
                    scale_x_continuous(expand=c(0,0),breaks=c(0, 0.25,0.5, 0.75,1),
                                        labels= c("0","0.25","0.5","0.75","1"))   +  
                    scale_y_continuous(expand=c(0,0),breaks=c(0, 0.25,0.5, 0.75,1),
                                       labels= c("0","0.25","0.5","0.75","1"),
                                       limits=c(0,0))   +  
                    scale_colour_manual(values =wes_palette("Chevalier")[c(1,4)], name='',
                                        labels =  expression(H[1] , H[2]))


plot_cdf <- ggplot(df, aes(x=t, y= CDF,color=test)) + 
                    geom_line() +
                    xlab("p-value")+
                    ylab("CDF")+
                    expand_limits(x = 0,y=0)  +
                    scale_x_continuous(expand=c(0,0),breaks=c(0, 0.25,0.5, 0.75,1),
                                        labels= c("0","0.25","0.5","0.75","1"))   +  
                    scale_y_continuous(expand=c(0,0),breaks=c(0, 0.25,0.5, 0.75,1),
                                       labels= c("0","0.25","0.5","0.75","1"))   +  
                    scale_colour_manual(values =wes_palette("Chevalier")[c(1,4)], name='',
                                          labels =  expression(H[1] , H[2]))

## ---- fig.width=6, fig.height=4, warning=F-------------------------------
plotgrid <- plot_grid(plot_cdf1, plot_cdf, plot_tdr_inversed1, plot_tdr_inversed, 
                      labels = paste0(letters[1:4],")"), ncol=2, align="hv")
plotgrid

## ---- eval=F-------------------------------------------------------------
#  ggsave(plotgrid,filename = "tdr_size_investing.pdf", width=6, height=4)

