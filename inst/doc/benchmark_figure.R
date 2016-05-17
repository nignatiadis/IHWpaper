## ----warning=F, message=F------------------------------------------------
library("IHW")
library("dplyr")
library("ggplot2")
library("grid")
library("tidyr")
library("cowplot")
library("RColorBrewer")
library("scales")
library("IHWpaper")

## ------------------------------------------------------------------------
# http://beyoncepalettes.tumblr.com/
beyonce_colors <- c("#b72da0", "#7c5bd2", "#0097ed","#00c6c3",
                   "#9cd78a", "#f7f7a7", "#ebab5f", "#e24344",
                   "#04738d")#,"#d8cdc9")
beyonce_colors[6] <- c("#dbcb09") # thicker yellow
barplot(rep(1,10), col=beyonce_colors)
pretty_colors <- beyonce_colors[c(2,1,3:5)]
pretty_names <- c("IHW", "BH", "Indep. Filt. \n 10 kb", "Indep. Filt. \n 200 kb", "Indep. Filt. \n 1 Mb")
pretty_names <- c("IHW", "BH", "10 kb", "200 kb", "1 Mb")

## ------------------------------------------------------------------------

rnaseq_file <- system.file("real_data_examples/result_files", "RNAseq_benchmark.Rds", package = "IHWpaper")
rnaseq_data <- readRDS(file=rnaseq_file)
panel_a_data <- group_by(rnaseq_data$alpha_df, alpha) %>% summarize(BH = max(bh_rejections), IHW=max(rejections)) %>% 
                     gather(method, rejections, BH, IHW) %>%
                     mutate(method = factor(as.character(method), levels=pretty_names))

## ----fig.width=5, fig.height=5-------------------------------------------
last_vals_a <- group_by(panel_a_data, method) %>% 
               summarize(last_vals = max(rejections)) %>%
               mutate(label = method,
                   colour = pretty_colors[match(label, pretty_names)])


panel_a <- ggplot(panel_a_data, aes(x=alpha,y=rejections,col=method)) +  
                geom_line(size=1.2) +
                xlab(expression(bold(paste("Nominal ",alpha)))) +
                ylab("Discoveries") +
                scale_color_manual(values=pretty_colors)+
                scale_x_continuous(expand=c(0,0), breaks=c(0.06,0.08,0.10)) +
                theme(plot.margin = unit(c(2, 1, 1, 1), "lines"))+
                theme(axis.title = element_text(face="bold" ))


panel_a <- pretty_legend(panel_a, last_vals_a, 0.102)
panel_a

## ------------------------------------------------------------------------
ggsave(panel_a + theme(plot.margin = unit(c(1, 3, 1, 0), "lines"))
      ,filename="DESeq2_benchmark.pdf", width = 4, height=4)

## ------------------------------------------------------------------------
breaks <- rnaseq_data$breaks
break_min <- rnaseq_data$break_min

breaks_left <- c(break_min,breaks[-length(breaks)])

step_df <- mutate(rnaseq_data$alpha_df, baseMean_left = breaks_left[stratum], baseMean_right = breaks[stratum],
                 baseMean_ratio = baseMean_right/baseMean_left , 
                 baseMean_left = baseMean_left * baseMean_ratio^.2,
                 baseMean_right = baseMean_right *baseMean_ratio^(-.2))

stratum_fun <- function(df){
  stratum <- df$stratum
  weight <- df$weight
  stratum_left <- stratum[stratum != length(stratum)]
  weight_left  <- weight[stratum_left]
  baseMean_left <- df$baseMean_right[stratum_left]
  stratum_right <- stratum[stratum != 1]
  weight_right <- weight[stratum_right]
  baseMean_right <- df$baseMean_left[stratum_right]
  data.frame(stratum_left= stratum_left, weight_left= weight_left, 
             stratum_right = stratum_right, weight_right = weight_right,
             baseMean_left = baseMean_left, baseMean_right = baseMean_right)
}

connect_df <- filter(step_df, alpha==0.1) %>% group_by(fold) %>% 
                  do(stratum_fun(.)) %>%
                  mutate(dashed = factor(ifelse(abs(weight_left - weight_right) > 10^(-4) , TRUE, FALSE),
                         levels=c(FALSE,TRUE)))


panel_b <- ggplot(filter(step_df, alpha==0.1), 
                       aes(x=baseMean_left, 
                           xend=baseMean_right, 
                           y=weight, yend=weight, col=fold)) +
                geom_segment(size=0.8)+ 
                geom_segment(data= connect_df, aes(x=baseMean_left, xend=baseMean_right, 
                                                y=weight_left, yend=weight_right, 
                                                linetype=dashed),
                             size=0.8)+
                scale_x_log10(breaks=c(1,10,100,1000,10000))+
                xlab("Mean of normalized counts")+
                ylab("Weight")+
                theme(legend.position=c(0.8,0.4)) +
                theme(plot.margin = unit(c(2, 1, 1, 2), "lines")) +
                guides(linetype=FALSE)+
                scale_color_manual(values=pretty_colors)+
                theme(axis.title = element_text(face="bold" ))


panel_b

## ------------------------------------------------------------------------
proteomics_file <- system.file("real_data_examples/result_files", "proteomics_benchmark.Rds", package = "IHWpaper")
proteomics_data <- readRDS(file=proteomics_file)
panel_c_data <- group_by(proteomics_data$alpha_df, alpha) %>% summarize(BH = max(bh_rejections), IHW=max(rejections)) %>% 
                     gather(method, rejections, BH, IHW) %>%
                     mutate(method = factor(as.character(method), levels=pretty_names))

## ------------------------------------------------------------------------
last_vals_c <- group_by(panel_c_data, method) %>% 
               summarize(last_vals = max(rejections)) %>%
               mutate(label = method,
                   colour = pretty_colors[match(label, pretty_names)])


panel_c <- ggplot(panel_c_data, aes(x=alpha,y=rejections,col=method)) +  
                geom_line(size=1.2) +
                xlab(expression(bold(paste("Nominal ",alpha)))) +
                ylab("Discoveries") +
                scale_color_manual(values=pretty_colors)+
                scale_x_continuous(expand=c(0,0), breaks=c(0.06,0.08,0.10))+
                theme(plot.margin = unit(c(2, 1, 1, 2), "lines"))+
                theme(axis.title = element_text(face="bold" ))


panel_c <- pretty_legend(panel_c, last_vals_c, 0.102)
panel_c

## ------------------------------------------------------------------------
breaks <- proteomics_data$breaks
break_min <- proteomics_data$break_min

breaks_left <- c(break_min,breaks[-length(breaks)])

step_df <- mutate(filter(proteomics_data$alpha_df , alpha==0.1), break_left = breaks_left[stratum],
                 break_right = breaks[stratum],
                 break_ratio = break_right/break_left , 
                 break_left =break_left * break_ratio^.1,
                 break_right = break_right *break_ratio^(-.1))

stratum_fun <- function(df){
  stratum <- df$stratum
  weight <- df$weight
  stratum_left <- stratum[stratum != length(stratum)]
  weight_left  <- weight[stratum_left]
  break_left <- df$break_right[stratum_left]
  stratum_right <- stratum[stratum != 1]
  weight_right <- weight[stratum_right]
  break_right <- df$break_left[stratum_right]
  data.frame(stratum_left= stratum_left, weight_left= weight_left, 
             stratum_right = stratum_right, weight_right = weight_right,
             break_left = break_left, break_right = break_right)
}

connecting_df <- step_df %>% group_by(fold) %>% 
                  do(stratum_fun(.)) %>%
                  mutate(dashed = factor(ifelse(abs(weight_left - weight_right) > 10^(-5) , TRUE, FALSE),
                         levels=c(FALSE,TRUE)))

panel_d <- ggplot(step_df, aes(x=break_left, xend=break_right,y=weight, yend=weight, col=fold)) +
                geom_segment(size=0.8)+                
                geom_segment(data= connecting_df, aes(x=break_left, xend=break_right, 
                                                y=weight_left, yend=weight_right, 
                                                linetype=dashed),
                             size=0.8)+
                scale_x_log10(breaks=c(1, 10,100,400)) +
                xlab("Number of peptides quantified")+
                ylab("Weight")+
                theme(legend.position=c(0.8,0.4)) +
                theme(plot.margin = unit(c(2, 1, 1, 2), "lines"))+
                guides(linetype=FALSE)+
                scale_color_manual(values=pretty_colors)+
                theme(axis.title = element_text(face="bold" ))


panel_d

## ------------------------------------------------------------------------
hqtl_file <- system.file("real_data_examples/result_files", "hQTL_benchmark.Rds", package = "IHWpaper")
hqtl_data <- readRDS(file=hqtl_file)
hqtl_summary <- group_by(hqtl_data$alpha_df, alpha) %>%  gather(method, rejections, 6:10) %>%
                select(alpha, method, rejections)

# will use below later to convert to nicer names for our plot
# see http://stackoverflow.com/questions/7547597/dictionary-style-replace-multiple-items-in-r
map <- setNames(pretty_names,
               c("rejections","bh_rejections","threshold:10000","threshold:2e+05", "threshold:1e+06"))

## ------------------------------------------------------------------------

hqtl_summary <- mutate(hqtl_summary, method = map[method], 
                       method= factor(method, levels = pretty_names))

last_vals_e <- group_by(hqtl_summary, method) %>% 
               summarize(last_vals = max(rejections))  %>%
               mutate(last_vals = last_vals + c(0,0, 250, 0, -350), # offset to make it look nice
                      label = method,
                      colour = pretty_colors[match(label, pretty_names)])


panel_e <- ggplot(hqtl_summary, aes(x=alpha,y=rejections,col=method)) +  
                geom_line(size=1.2) +
                xlab(expression(bold(paste("Nominal ",alpha)))) +
                ylab("Discoveries") +
                scale_x_continuous(expand=c(0,0), breaks=c(0.06,0.08,0.10))+
                scale_color_manual(values=pretty_colors)+
                theme(plot.margin = unit(c(2, 1, 1, 2), "lines"))+
                theme(axis.title = element_text(face="bold" )) 



panel_e <- pretty_legend(panel_e, last_vals_e, 0.102)
panel_e 

## ------------------------------------------------------------------------
ggsave(panel_e + theme(plot.margin = unit(c(1, 3, 1, 0), "lines"))
      ,filename="hqtl_benchmark.pdf", width = 4, height=4)

## ------------------------------------------------------------------------
breaks <-   hqtl_data$breaks
breaks <- breaks[-1]
break_min <- hqtl_data$break_min
breaks_left <- c(break_min,breaks[-length(breaks)])
step_df <- mutate(filter(hqtl_data$alpha_df,alpha==0.1), break_left = breaks_left[stratum],
                 break_right = breaks[stratum],
                 break_ratio = break_right/break_left , 
                 break_left =break_left * break_ratio^.2,
                 break_right = break_right *break_ratio^(-.2))

stratum_fun <- function(df){
  stratum <- df$stratum
  weight <- df$weight
  stratum_left <- stratum[stratum != length(stratum)]
  weight_left  <- weight[stratum_left]
  break_left <- df$break_right[stratum_left]
  stratum_right <- stratum[stratum != 1]
  weight_right <- weight[stratum_right]
  break_right <- df$break_left[stratum_right]
  data.frame(stratum_left= stratum_left, weight_left= weight_left, 
             stratum_right = stratum_right, weight_right = weight_right,
             break_left = break_left, break_right = break_right)
}

connecting_df <- step_df %>% group_by(fold) %>% 
                  do(stratum_fun(.)) %>%
                  mutate(dashed = factor(ifelse(abs(weight_left - weight_right) > 2 , TRUE, FALSE),
                         levels=c(FALSE,TRUE)))

panel_f <- ggplot(step_df, aes(x=break_left, xend=break_right,y=weight, yend=weight, col=fold)) +
                geom_segment(size=0.8)+                
                geom_segment(data= connecting_df, aes(x=break_left, xend=break_right, 
                                                y=weight_left, yend=weight_right, 
                                                linetype=dashed),
                             size=0.8)+
                scale_x_log10(breaks=c(10^4, 10^5,10^6,10^7), 
                              labels = trans_format("log10", math_format(10^.x))) +
                xlab("Genomic distance (bp)")+
                ylab("Weight")+
                theme(legend.position=c(0.8,0.6)) +
                theme(plot.margin = unit(c(2, 1.5, 1, 2.5), "lines"))+
                theme(axis.title = element_text(face="bold" ))+
                scale_color_manual(values=pretty_colors)+
                guides(linetype=FALSE)

panel_f


## ----fig.width=12, fig.height=3------------------------------------------
top_row_main_fig <- plot_grid(panel_a, 
                      panel_c, 
                      panel_e, panel_f,
                      labels= c("a)", "b)", "c)",
                              "d)"),
                      nrow=1,
                      hjust=-3)
top_row_main_fig

## ----eval=FALSE----------------------------------------------------------
#  ggsave(plot=full_fig, file="data_examples.pdf", width=12, height=10)

## ----fig.width=12, fig.height=5------------------------------------------
suppl_fig <- plot_grid(panel_b,
                       panel_d,
                       labels= c("a)", "b)"),
                       nrow=1)
suppl_fig

## ----eval=FALSE----------------------------------------------------------
#  ggsave(plot=suppl_fig, file="data_examples_suppl.pdf", width=12, height=5)

## ----warning=F, message=F------------------------------------------------
library("ggplot2")
library("grid")
library("dplyr")
library("cowplot")
library("IHWpaper")

## ------------------------------------------------------------------------
methods_pretty <- c("BH", "Clfdr", "Greedy Indep. Filt.", "IHW", "IHW naive", "FDRreg", "LSL GBH", "SBH", "TST GBH")
colors <- rep(NA,9 )

conservative_methods <- c("BH", "Clfdr", "IHW", "FDRreg", "LSL GBH")
conservative_idx <- match(conservative_methods, methods_pretty)
colors[conservative_idx] <- beyonce_colors[c(1,3,2,4,5)]

anticonservative_methods <- c("Greedy Indep. Filt.", "IHW naive", "SBH", "TST GBH")
anticonservative_idx     <- match(anticonservative_methods, methods_pretty)
colors[anticonservative_idx] <- beyonce_colors[1:4]


## ------------------------------------------------------------------------
null_grb_file <- system.file("simulation_benchmarks/result_files",
                        "ihw_null_simulation_benchmark_grb.Rds", package = "IHWpaper")
null_e3_file <- system.file("simulation_benchmarks/result_files",
                        "ihw_null_simulation_benchmark_E3.Rds", package = "IHWpaper")
null_file <- system.file("simulation_benchmarks/result_files",
                        "ihw_null_simulation_benchmark.Rds", package = "IHWpaper")
null_df <- rbind(readRDS(null_grb_file),
                 readRDS(null_e3_file),
                 readRDS(null_file)) %>% 
                 filter(fdr_method != "IHW") %>% # just show IHW with E1-E3 
                 mutate(fdr_method = ifelse(fdr_method=="IHW E3", "IHW", fdr_method),  
                           fdr_method = sapply(strsplit(fdr_method," 20"), "[",1),
                           fdr_method = factor(fdr_method, levels= methods_pretty)) 



## ------------------------------------------------------------------------
panel_a_sim_df <-  filter(null_df, fdr_method %in% anticonservative_methods)
#  only for panel a) rename method_pretty to break across lines
levels(panel_a_sim_df$fdr_method)[levels(panel_a_sim_df$fdr_method) ==  "Greedy Indep. Filt."] <- "Greedy"

last_vals_a_sim <- group_by(panel_a_sim_df, fdr_method) %>% summarize(last_vals = max(FDR)) %>%
                    mutate(last_vals = last_vals + c(-0.03,0, 0.1, -0.1),#c(0, 0,0.05, -0.03), 
                    label = fdr_method,
                    colour = colors[anticonservative_idx])

panel_a_sim <- ggplot(panel_a_sim_df, aes(x=alpha, y=FDR, col=fdr_method)) +
                         geom_line(size=1.2) +
                         geom_abline(linetype="dashed") + 
                         xlab(expression(bold(paste("Nominal ",alpha)))) + 
                         scale_x_continuous(limits= c(0.01,0.1), breaks=seq(0.01,0.1,length=4)) +
                         ylim(0,0.9) +
                         theme(plot.margin = unit(c(2, 3, 1, 1), "lines")) +
                         scale_color_manual(values=colors[anticonservative_idx])+
                         theme(axis.title = element_text(face="bold") )

panel_a_sim <- panel_a_sim + annotation_custom(grob=textGrob("Indep. Filt.", hjust=0,
                                               gp=gpar(fontsize=13, 
                                                       col=last_vals_a_sim$colour[1])),
                                               xmin=0.102,
                                               xmax=0.102,
                                               ymin=0.5,
                                               ymax=0.5)

panel_a_sim <- pretty_legend(panel_a_sim, last_vals_a_sim, 0.102)

panel_a_sim


## ------------------------------------------------------------------------
panel_b_sim_df <- filter(null_df, fdr_method %in% conservative_methods)

last_vals_b_sim <- group_by(panel_b_sim_df, fdr_method) %>% summarize(last_vals = max(FDR)) %>% 
                    mutate(last_vals = last_vals + c(0.0085,0.007,-0.0085, 0 ,-0.007 ),  #+  c(0.005,0.005,-0.005, 0 ,-0.005 ), 
                    label = fdr_method,
                    colour = colors[conservative_idx])


panel_b_sim <- ggplot(panel_b_sim_df, aes(x=alpha, y=FDR, col=fdr_method)) +
                         geom_abline(linetype="dashed") + 
                         geom_line(size=1.2) +
                         xlab(expression(bold(paste("Nominal ",alpha)))) + 
                         scale_x_continuous(limits= c(0.01,0.1), breaks=seq(0.01,0.1,length=4)) +
                         theme(plot.margin = unit(c(2, 2.3, 1, 2.5), "lines")) +
                         scale_color_manual(values=colors[conservative_idx])+
                         theme(axis.title = element_text(face="bold") )


panel_b_sim <- pretty_legend(panel_b_sim, last_vals_b_sim, 0.102 )
panel_b_sim

## ------------------------------------------------------------------------
effsize_grb_file <- system.file("simulation_benchmarks/result_files",
                        "ihw_du_ttest_inform_simulation_benchmark_grb.Rds", package = "IHWpaper")
effsize_e3_file <- system.file("simulation_benchmarks/result_files",
                        "ihw_du_ttest_inform_simulation_benchmark_E3.Rds", package = "IHWpaper")
effsize_file <- system.file("simulation_benchmarks/result_files",
                         "ihw_du_ttest_inform_simulation_benchmark.Rds", package = "IHWpaper")
effsize_df <- rbind(readRDS(effsize_grb_file),
                    readRDS(effsize_file),
                    readRDS(effsize_e3_file)) %>%
              filter(fdr_method != "IHW") %>% # just show IHW with E1-E3 
              mutate(fdr_method = ifelse(fdr_method=="IHW E3", "IHW", fdr_method),  
                           fdr_method = sapply(strsplit(fdr_method," 20"), "[",1),
                           fdr_method = factor(fdr_method, levels= methods_pretty)) 

## ------------------------------------------------------------------------
panel_c_sim_df <- filter(effsize_df, fdr_method %in% conservative_methods)

last_vals_c_sim <- group_by(panel_c_sim_df, fdr_method) %>% summarize(last_vals =  FDR[which.max(eff_size)]) %>%
               mutate(last_vals = last_vals + c(0.010,0.005,-0.015, 0 , -0.010),  #+ c(0,0.005,-0.005, 0 ,-0.01 ), 
                      label = fdr_method,
                      colour = colors[conservative_idx])

panel_c_sim <- ggplot(panel_c_sim_df, aes(x=eff_size, y=FDR, col=fdr_method)) +
                         geom_hline(yintercept=0.1, linetype="dashed") + 
                         geom_line(size=1.2) +
                         xlab("Effect size") + 
                         theme(plot.margin = unit(c(2, 3, 1, 2), "lines")) +
                         scale_color_manual(values=colors[conservative_idx])+
                         theme(axis.title = element_text(face="bold") )


panel_c_sim <- pretty_legend(panel_c_sim, last_vals_c_sim, 2.52 )
panel_c_sim

## ------------------------------------------------------------------------
panel_d_sim_df <- filter(effsize_df, fdr_method %in% conservative_methods)

last_vals_d_sim <- group_by(panel_d_sim_df, fdr_method) %>% summarize(last_vals = power[which.max(eff_size)]) %>%
                      mutate(last_vals = last_vals + c(0,0,-0.085, 0.075 , -0.014 ), #+ c(0,-0.015,-0.035, 0.035 ,+0.005 ), 
                      label = fdr_method,
                      colour = colors[conservative_idx])


panel_d_sim <- ggplot(panel_c_sim_df, aes(x=eff_size, y=power, col=fdr_method)) +
                         geom_line(size=1.2) +
                         xlab("Effect size") + 
                         ylab("Power")+
                         theme(plot.margin = unit(c(2, 3.5, 1, .5), "lines")) +
                         scale_color_manual(values=colors[conservative_idx])+
                         theme(axis.title = element_text(face="bold") )

panel_d_sim <- pretty_legend(panel_d_sim, last_vals_d_sim, 2.52 )
panel_d_sim

## ------------------------------------------------------------------------
sizeinvesting_grb_file <- system.file("simulation_benchmarks/result_files",
                        "ihw_wasserman_normal_simulation_benchmark_grb.Rds", package = "IHWpaper")
sizeinvesting_e3_file <-  system.file("simulation_benchmarks/result_files",
                        "ihw_wasserman_normal_simulation_benchmark_E3.Rds", package = "IHWpaper")
sizeinvesting_file <- system.file("simulation_benchmarks/result_files",
                         "ihw_wasserman_normal_simulation_benchmark.Rds", package = "IHWpaper")
sizeinvesting_df <- rbind_all(lapply(c(sizeinvesting_file,
                                       sizeinvesting_e3_file,
                                       sizeinvesting_grb_file), readRDS)) %>%
                    filter(fdr_method != "IHW") %>% # just show IHW with E1-E3 
                    # make names prettier and make sure we use same factor for everything
                    mutate(fdr_method = ifelse(fdr_method=="IHW E3", "IHW", fdr_method),  
                           fdr_method = sapply(strsplit(fdr_method," 20"), "[",1),
                           fdr_method = factor(fdr_method, levels= methods_pretty)) %>% 
                    # add log2 relative to BH
                    group_by(xi_max) %>% 
                    mutate(normalized = log2(power/max(power*(fdr_method=="BH"))))


## ------------------------------------------------------------------------
panel_e_sim_df <- filter(sizeinvesting_df, fdr_method %in% conservative_methods)

last_vals_e_sim <- group_by(panel_e_sim_df, fdr_method) %>% summarize(last_vals = FDR[which.max(xi_max)]) %>%
                      mutate(last_vals = last_vals + c(0.0009,0,-0.0009, 0 ,0 ), 
                      label = fdr_method,
                      colour = colors[conservative_idx])

panel_e_sim <- ggplot(panel_e_sim_df, aes(x=xi_max, y=FDR, col=fdr_method)) +
                         geom_hline(yintercept=0.1, linetype="dashed") + 
                         geom_line(size=1.2) +
                         xlab(expression(bold(xi[max])))+
                         theme(plot.margin = unit(c(3, 7.5, .2, .2), "lines"))+
                         scale_color_manual(values=colors[conservative_idx])+
                         theme(axis.title = element_text(face="bold") )


panel_e_sim <- pretty_legend(panel_e_sim, last_vals_e_sim, 6.02 )
panel_e_sim

## ------------------------------------------------------------------------
panel_f_sim_df <- filter(sizeinvesting_df, fdr_method %in% conservative_methods)

last_vals_f_sim <- group_by(panel_f_sim_df, fdr_method) %>% summarize(last_vals = normalized[which.max(xi_max)]) %>%
                      mutate(last_vals = last_vals +  c(0, 0, 0, +0.008 ,-0.002 ), 
                      label = fdr_method,
                      colour = colors[conservative_idx])

panel_f_sim <- ggplot(panel_f_sim_df, aes(x=xi_max, y=normalized, col=fdr_method)) +
                         geom_line(size=1.2) +
                         ylab(expression(bold(log[2](Power/Power[BH]))))+
                         xlab(expression(bold(xi[max])))+
                         theme(plot.margin = unit(c(3, 7.5, .2, .2), "lines"))+
                         scale_color_manual(values=colors[conservative_idx])+
                         theme(axis.title = element_text(face="bold") )


panel_f_sim <- pretty_legend(panel_f_sim, last_vals_f_sim, 6.02 )
panel_f_sim

## ----fig.width=12, fig.height=3------------------------------------------
bottom_row_main_fig <- plot_grid(panel_a_sim,panel_b_sim,
                          panel_c_sim,panel_d_sim,
                          nrow=1,
                          labels= c("d)","e)","f)","g)"))

bottom_row_main_fig

## ----fig.width=12, fig.height=6------------------------------------------

benchmark_fig <- plot_grid(panel_a,panel_c,
                          panel_e,panel_f,
                          panel_a_sim,panel_b_sim,
                          panel_c_sim,panel_d_sim,
                          nrow=2,
                          labels= c("a)", "b)", "c)",
                              "d)","e)","f)","g)","h)"),
                          hjust=-3)

benchmark_fig

## ----eval=FALSE----------------------------------------------------------
#  ggsave(plot=benchmark_fig, file="main_simulations.pdf", width=12, height=6)

## ------------------------------------------------------------------------
# Supplementary Effect Size Simulation

## Supplementary panel a)


## ------------------------------------------------------------------------
sup_panel_a_df <- filter(effsize_df, fdr_method %in% anticonservative_methods)

sup_last_vals_a <- group_by(sup_panel_a_df, fdr_method) %>%
                   summarize(last_vals =  FDR[which.max(eff_size)]) %>%
                   mutate(last_vals = last_vals + c(-0.005 ,+0.016,0, -0.002 ), 
                      label = fdr_method,
                      colour = colors[anticonservative_idx])

sup_panel_a <- ggplot(sup_panel_a_df, aes(x=eff_size, y=FDR, col=fdr_method)) +
                         geom_hline(yintercept=0.1, linetype="dashed") + 
                         geom_line(size=1.2) +
                         xlab("Effect size") + 
                         theme(plot.margin = unit(c(3, 7.5, .2, .2), "lines"))+
                         scale_color_manual(values=colors[anticonservative_idx])+
                         theme(axis.title = element_text(face="bold") )


sup_panel_a <- pretty_legend(sup_panel_a, sup_last_vals_a, 2.52 )
sup_panel_a

## ------------------------------------------------------------------------
sup_panel_b_df <- filter(effsize_df, fdr_method %in% anticonservative_methods)

sup_last_vals_b <- group_by(sup_panel_b_df, fdr_method) %>%
                   summarize(last_vals = power[which.max(eff_size)]) %>%
                   mutate(last_vals = last_vals + c(+0.015,+0.03,-0.06, -0.025 ), 
                      label = fdr_method,
                      colour = colors[anticonservative_idx])


sup_panel_b <- ggplot(sup_panel_b_df, aes(x=eff_size, y=power, col=fdr_method)) +
                         geom_line(size=1.2) +
                         xlab("Effect size") + 
                         ylab("Power")+
                         theme(plot.margin = unit(c(3, 7.5, .2, .2), "lines"))+
                         scale_color_manual(values=colors[anticonservative_idx])+
                         theme(axis.title = element_text(face="bold") )


sup_panel_b <- pretty_legend(sup_panel_b, sup_last_vals_b, 2.52 )
sup_panel_b

## ------------------------------------------------------------------------
sup_panel_c_df <- filter(sizeinvesting_df, fdr_method %in% anticonservative_methods)

sup_last_vals_c <- group_by(sup_panel_c_df, fdr_method) %>% 
                   summarize(last_vals = FDR[which.max(xi_max)]) %>%
                   mutate(last_vals = last_vals +  c(0.00005,0, 0.001,-0.0015), 
                      label = fdr_method,
                      colour = colors[anticonservative_idx])

sup_panel_c <- ggplot(sup_panel_c_df, aes(x=xi_max, y=FDR, col=fdr_method)) +
                         geom_hline(yintercept=0.1, linetype="dashed") + 
                         geom_line(size=1.2) +
                         xlim(3,6)+
                         xlab(expression(bold(xi[max])))+
                         theme(plot.margin = unit(c(3, 7.5, .2, .2), "lines"))+
                         scale_color_manual(values=colors[anticonservative_idx])+
                         theme(axis.title = element_text(face="bold") )


sup_panel_c <- pretty_legend(sup_panel_c, sup_last_vals_c, 6.02 )
sup_panel_c

## ------------------------------------------------------------------------
sup_panel_d_df <- filter(sizeinvesting_df, fdr_method %in% anticonservative_methods)

sup_last_vals_d <- group_by(sup_panel_d_df, fdr_method) %>% 
                  summarize(last_vals = normalized[which.max(xi_max)]) %>%
                  mutate(last_vals = last_vals + c(0.00,0.0, +0.005, -0.005), 
                      label = fdr_method,
                      colour = colors[anticonservative_idx])

sup_panel_d <- ggplot(sup_panel_d_df, aes(x=xi_max, y=normalized, col=fdr_method)) +
                         geom_line(size=1.2) +
                         ylab(expression(bold(log[2](Power/Power[BH]))))+
                         xlab(expression(bold(xi[max])))+
                         xlim(3,6)+
                         ylim(-0.1,+0.3)+
                         theme(plot.margin = unit(c(3, 7.5, .2, .2), "lines"))+
                         scale_color_manual(values=colors[anticonservative_idx])+
                         theme(axis.title = element_text(face="bold") )


sup_panel_d <- pretty_legend(sup_panel_d, sup_last_vals_d, 6.02 )
sup_panel_d

## ----fig.width=14, fig.height=10-----------------------------------------
sup_sim_fig <- plot_grid(sup_panel_a,
                          sup_panel_c,
                          panel_e_sim,
                          sup_panel_b,
                          sup_panel_d,
                          panel_f_sim,
                          ncol=3,
                          labels= c("a)", "c)", "e)",
                              "b)","d)","f)"))
  sup_sim_fig

## ----eval=FALSE----------------------------------------------------------
#  cowplot::ggsave(plot=sup_sim_fig, file="suppl_simulations.pdf", width=14, height=10)

