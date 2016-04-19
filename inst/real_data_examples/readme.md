*IHWpaper/inst/extdata/real_data/bottomly_eset.RData
From Recount project (http://bowtie-bio.sourceforge.net/recount/)

*IHWpaper/inst/extdata/real_data/hammer_eset.RData
From Recount project (http://bowtie-bio.sourceforge.net/recount/)

*IHWpaper/inst/extdata/real_data/science_signaling.csv
Proteomics example, Csv extracted from "2002548TableS1.xlsx" table in supplementary materials for:
Hyperplexing: A Method for Higher-Order Multiplexed Quantitative Proteomics Provides 
a Map of the Dynamic Response to Rapamycin in Yeast
Noah Dephoure and Steven P. Gygi*
(http://stke.sciencemag.org/content/5/217)

*IHWpaper/inst/extdata/real_data/hqtl_pvalue_filtered.Rds
For the hQTL example, here only p-values <= 0.005 are retained and the columns
pvalue, dist and group  (pvalue for the snp-peak pair, distance of that pair
and stratum into which it was categorized based on distance). The original number
of hypotheses in each stratum are stored in the attribute "m_groups".


```r
#full data-frame with all hypotheses (a few GB), not made available here

hqtl <- readRDS("qtls_chrom_21.Rds")

my_breaks <- c(-1,
               seq(from=10000,to=290000, by=10000) , 
               seq(from=300000, to=0.9*10^6, by=100000),
               seq(from=10^6, to=50*10^6, by=10^7))

hqtl <- mutate(hqtl, group=cut(dist, my_breaks))

m_groups <- table(hqtl$group)

hqtl_filt <- filter(hqtl, pvalue <= 5*10^(-3)) %>%
             select(pvalue, dist, group)

attr(hqtl_filt, "m_groups") <- m_groups
attr(hqtl_filt, "breaks") <- my_breaks


saveRDS(hqtl_filt, file="hqtl_pvalue_filtered.Rds")
```