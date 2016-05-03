## ----warning=F, message=F------------------------------------------------
library("ggplot2")
library("grid")
library("dplyr")
library("cowplot")
library("IHWpaper")
library("wesanderson")

## ------------------------------------------------------------------------
colors <- wes_palette("Cavalcanti")[c(3,5)]

## ------------------------------------------------------------------------

null_folder <- system.file("simulation_benchmarks/result_files/ihw_bonf_null",
                        package = "IHWpaper")
null_sim <- rbind_all(lapply(file.path(null_folder,list.files(null_folder)), function(x) readRDS(x))) %>%
            mutate(fdr_method = ifelse(fdr_method == "IHW-Bonferroni E3", "IHW-Bonferroni", fdr_method))

last_vals_a <- group_by(null_sim, fdr_method) %>% summarize(last_vals = max(FDR)) %>%
               mutate(last_vals = last_vals + c(+0.003,-0.003), 
                    label = fdr_method,
                    colour = colors)

panel_a <- ggplot(null_sim, aes(x=alpha, y=FDR, col=fdr_method)) +
                         geom_line(size=1.2) +
                         geom_abline(linetype="dashed") + 
                         xlab(expression(bold(paste("Nominal ",alpha)))) + 
                         ylab("FWER")+
                         scale_x_continuous(limits= c(0.005,0.1), breaks=seq(0.01,0.09,length=3)) +
                         scale_y_continuous(limits= c(0.005,0.1), breaks=seq(0.01,0.09,length=3)) +
                         theme(plot.margin = unit(c(3, 7.5, .2, .2), "lines"))+
                         scale_color_manual(values=colors)+
                         theme(axis.title = element_text(face="bold") )

  
panel_a <- pretty_legend(panel_a, last_vals_a, 0.102)
panel_a

## ------------------------------------------------------------------------
ttest_folder <- system.file("simulation_benchmarks/result_files/ihw_bonf_du_ttest_informative",
                     package = "IHWpaper")
ttest_sim <- rbind_all(lapply(file.path(ttest_folder,list.files(ttest_folder)), function(x) readRDS(x))) %>%
            mutate(fdr_method = ifelse(fdr_method == "IHW-Bonferroni E3", "IHW-Bonferroni", fdr_method))

last_vals_b <- group_by(ttest_sim, fdr_method) %>% summarize(last_vals =  FDR[which.max(eff_size)]) %>%
               mutate(last_vals = last_vals, 
                      label = fdr_method,
                      colour = colors)

panel_b <- ggplot(ttest_sim, aes(x=eff_size, y=FDR, col=fdr_method)) +
                         geom_hline(yintercept=0.1, linetype="dashed") + 
                         geom_line(size=1.2) +
                         xlab("Effect size") + 
                          ylab("FWER")+
                         theme(plot.margin = unit(c(3, 7.5, .2, .2), "lines"))+
                         scale_color_manual(values=colors) + 
                         theme(axis.title = element_text(face="bold") )

panel_b <- pretty_legend(panel_b, last_vals_b, 2.52 )
panel_b

## ------------------------------------------------------------------------
last_vals_c <- group_by(ttest_sim, fdr_method) %>% summarize(last_vals = power[which.max(eff_size)]) %>%
               mutate(last_vals = last_vals, 
                      label = fdr_method,
                      colour = colors)


panel_c <- ggplot(ttest_sim, aes(x=eff_size, y=power, col=fdr_method)) +
                         geom_line(size=1.2) +
                         xlab("Effect size") + 
                         ylab("Power")+
                         theme(plot.margin = unit(c(3, 7.5, .2, .2), "lines"))+
                         scale_color_manual(values=colors)+
                         theme(axis.title = element_text(face="bold") )

panel_c <- pretty_legend(panel_c, last_vals_c, 2.52 )
panel_c

## ------------------------------------------------------------------------
size_investing_folder <- system.file("simulation_benchmarks/result_files/ihw_bonf_wasserman_normal",
                              package = "IHWpaper")
size_investing_sim <- rbind_all(lapply(file.path(size_investing_folder,
                                             list.files(size_investing_folder)), function(x) readRDS(x))) %>%
            mutate(fdr_method = ifelse(fdr_method == "IHW-Bonferroni E3", "IHW-Bonferroni", fdr_method))


last_vals_d <- group_by(size_investing_sim, fdr_method) %>% summarize(last_vals = log10(FDR[which.max(xi_max)])) %>%
               mutate(last_vals = last_vals + c(+2.5,-0.8)*0.08, 
                      label = fdr_method,
                      colour = colors)

panel_d <- ggplot(size_investing_sim, aes(x=xi_max, y=log10(FDR), col=fdr_method)) +
                         geom_line(size=1.2) +
                         xlab(expression(bold(xi[max])))+
                         ylab("log10(FWER)")+
                         geom_hline(yintercept=-1, linetype="dashed") + 
                         theme(plot.margin = unit(c(3, 7.5, .2, .2), "lines"))+
                         scale_color_manual(values=colors) + 
                         theme(axis.title = element_text(face="bold") )

panel_d <- pretty_legend(panel_d, last_vals_d, 6.02 )
panel_d

## ------------------------------------------------------------------------
last_vals_e <- group_by(size_investing_sim, fdr_method) %>% summarize(last_vals = power[which.max(xi_max)]) %>%
               mutate(last_vals = last_vals + c(-0.01,+0.01) , 
                      label = fdr_method,
                      colour = colors)

panel_e <- ggplot(size_investing_sim, aes(x=xi_max, y=power, col=fdr_method)) +
                         geom_line(size=1.2) +
                         ylab("Power") +
                         xlab(expression(bold(xi[max])))+
                         theme(plot.margin = unit(c(3, 7.5, .2, .2), "lines"))+
                         scale_color_manual(values=colors) +
                         theme(axis.title = element_text(face="bold") )


panel_e <- pretty_legend(panel_e, last_vals_e, 6.02 )
panel_e

## ---- fig.width=12, fig.height=8-----------------------------------------
fwer_sim_fig <- plot_grid(panel_a, 
                         panel_b, 
                         panel_d, 
                         ggdraw(),
                         panel_c,
                         panel_e, 
                         nrow=2,
                         labels=c("a)", "b)", "d)","","c)", "e)"))
fwer_sim_fig

## ----eval=FALSE----------------------------------------------------------
#  ggsave(plot=fwer_sim_fig, file="fwer_simulations.pdf", width=12, height=8)

