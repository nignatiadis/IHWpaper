[![Travis-CI Build Status](https://travis-ci.org/nignatiadis/IHWpaper.svg?branch=master)](https://travis-ci.org/nignatiadis/IHWpaper)


# IHWpaper

This is a package which reproduces all IHW paper figures (as well as some additional figures for presentations etc.).


# Vignettes in the IHWpaper package

All vignettes can be found in the /vignettes directory. For the vignettes which reproduce figures in our manuscript, we also provide a link to the rendered html page below:

* Figure 1: [Stratified Histograms](https://rawgit.com/nignatiadis/IHWpaper/master/inst/doc/stratified_histograms.html)
* Figure 2: [Real data examples](https://rawgit.com/nignatiadis/IHWpaper/master/inst/doc/real_data_examples.html)
* Figures 3 and S1: [Simulation Results](https://rawgit.com/nignatiadis/IHWpaper/master/inst/doc/simulations_vignette.html)
* Figure 4: Schematic figure, no code used for its generation
* Figures S2 and S3: [Explaining tdr](https://rawgit.com/nignatiadis/IHWpaper/master/inst/doc/explaining_tdr.html)
* Figure S4: [IHW-Bonferroni simulations](https://rawgit.com/nignatiadis/IHWpaper/master/inst/doc/IHW_bonferroni_simulations.html)

# Installing the package and reproducing all simulations

You can install the package as follows:

```{r}
library("devtools")
# install IHW
install_github("vladchimescu/lpsymphony", subdir="lpsymphony")
install_github("nignatiadis/IHW")
# install FDRreg, version on CRAN is outdated
install_github(repo= "jgscott/FDRreg", subdir="R_pkg/", ref = "a63cebae6faecb1fb0ebee634195296f39faa11b")
# Bioconductor prerequisites
source("http://bioconductor.org/biocLite.R")
biocLite(c("genefilter","DESeq2","qvalue","Biobase",
            "BiocParallel","airway","pasilla", "BiocStyle"))
# finally install this package
install_github("nignatiadis/ihwPaper")
```

Afterwards you can view the pre-built vignettes or rebuild them yourself.

Notice that the vignettes load binary .Rds files which are stored in inst/real_data_examples/result_files
and inst/simulations_benchmarks/result_files. These have been generated by longer lasting jobs, the code for which can be found in inst/real_data_examples and inst/simulation_benchmarks.

