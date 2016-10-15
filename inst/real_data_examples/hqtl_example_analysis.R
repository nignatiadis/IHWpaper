################################################################################
# Note: This is the only example in the paper which is not reproducible yet,
# since many of the files are too large to include. We will provide a solution
# at a later time point, but include the code so all the steps can at least
# be followed.
################################################################################

##########################################################
# Step 1: Get the hQTL p-values using matrix eQTL
##########################################################

# redo judith's analysis but with only one chromosome
library("DESeq2")
library("preprocessCore")
library("MatrixEQTL")
library("peer")
library("genefilter")
library("dplyr")

## set the project folder - it has to contain the object fileAnnoDF.rda 
projectFolder='/g/huber/users/ignatiadis/judith_hqtl/'

## load the fileAnnotationObject
load(paste(projectFolder,'fileAnnoDF.rda',sep=''))

## define the folder for the results
resultDir=paste(projectFolder,'/results/',sep='')
if(!file.exists(resultDir))system(paste('mkdir',resultDir))
plotFolder = paste(resultDir, 'plots/',sep='')
if(!file.exists(plotFolder))system(paste('mkdir',plotFolder))

# load the covariates
covariates_file_used = fileAnnoDF$FileLocation[fileAnnoDF$Filename=='covariates_file_used']
load(covariates_file_used)


# pick smallest chromosome...
snpchr <- 21
mod    <- 'H3K27AC'

cat('processing',mod,'...\n')

## load the peak annotations for the respective mark
peak_location_file_name = fileAnnoDF$FileLocation[fileAnnoDF$Filename=='peak_location_file_name'&fileAnnoDF$Mark==mod]
cat('reading peaks position file',mod,'...\n')
peakpos_all = read.table(peak_location_file_name, header = TRUE, stringsAsFactors = FALSE);
peakpos_all$chr=gsub('chr','',peakpos_all$chr)

## load the hpeaks signals for the respective mark
hpeaksfileName = fileAnnoDF$FileLocation[fileAnnoDF$Filename=='hpeaksfileName'&fileAnnoDF$Mark==mod]
cat('reading peaks signal file',mod,'...\n')
load(hpeaksfileName)
hpeaksMat = as.matrix(hpeaks)

cat('processing chromosome',snpchr,'...\n')

## define your output file
QTLresultFile = paste(resultDir, 'cisQTLs_',mod,'_chr',snpchr,'.txt',sep='')
logFile = paste(QTLresultFile,'.log',sep='')


snps_location_file_name = fileAnnoDF$FileLocation[fileAnnoDF$Filename=='snps_location_file_name'&fileAnnoDF$Chr==snpchr]
SNP_file_name_used = fileAnnoDF$FileLocation[fileAnnoDF$Filename=='SNP_file_name_used'&fileAnnoDF$Chr==snpchr]

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
snpspos=snpspos[,1:3]

## use only the peaks of the chromosome in question
peakpos=peakpos_all[peakpos_all$chr%in%snpchr,]

## load the SNP genotype information file
load(SNP_file_name_used)

## match columns of SNPs, covariates and peaks (there is one individual too many in the hpeaks)
snps$ColumnSubsample(snps$columnNames%in%colnames(hpeaksMat))
cvrt2<-cvrt
cvrt2$ColumnSubsample(cvrt2$columnNames%in%snps$columnNames)
cvrt2$ColumnSubsample(match(snps$columnNames,cvrt2$columnNames))

hpeaksMatT = hpeaksMat[rownames(hpeaksMat)%in%peakpos$id,]

## redefine the hpeaks, for a given chromosome, this speeds up the QTL calling
hpeaks2 = SlicedData$new();
hpeaks2$fileDelimiter = " "; # the TAB character
hpeaks2$fileOmitCharacters = "NA"; # denote missing values;
hpeaks2$fileSkipRows = 1; # one row of column labels
hpeaks2$fileSkipColumns = 1; # one column of row labels
hpeaks2$fileSliceSize = 2000; # read file in pieces of 2,000 rows
hpeaks2$CreateFromMatrix(hpeaksMatT[match(peakpos$id,rownames(hpeaksMatT)),match(snps$columnNames,colnames(hpeaksMatT))]);
hpeaks2$ResliceCombined()


## Main function (matrix eQTL)
me = Matrix_eQTL_main(
  gene = hpeaks2, 
  snps = snps, 
  cvrt = cvrt2,
  output_file_name = QTLresultFile,
  pvOutputThreshold  = 1,
  useModel = modelLINEAR,  
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = NULL,
  pvOutputThreshold.cis = 0,#pvOutputThreshold,
  snpspos = snpspos, 
  genepos = peakpos[,1:4],
  cisDist = 0,
  pvalue.hist = 402)

saveRDS(me, file="chrom21H3K27AC") #save intermediate result


tmp <- readRDS("chrom21H3K27AC")
qtls <- tmp$all$eqtls
qtls$beta <- NULL
qtls$statistic <- NULL
qtls$FDR <- NULL
snpspos <- select(snpspos, snp, pos) %>% rename(snps = snp, snp_pos=pos)
qtls <- left_join(qtls, snpspos)

peakpos <- select(peakpos, id, start, end) %>% rename(gene=id, peak_start=start, peak_end=end)
peakpos <- mutate(peakpos, gene=as.factor(gene))

qtls <- left_join(qtls, peakpos)
# calculate distances
qtls <- mutate(qtls, dist = pmin(abs(peak_start-snp_pos), abs(peak_end-snp_pos)))
saveRDS(qtls, file="qtls_chrom_21.Rds")

##########################################################
# Step 2: Apply IHW and compare to BH/Indep. Filtering
##########################################################

library("IHW")


qtls <- readRDS(file = "qtls_chrom_21.Rds")
print("qtls loaded into memory")
  
## up to 300k in 10 bins
my_breaks <- c(-1, 
                 seq(from=10000,to=290000, by=10000) , 
                 seq(from=300000, to=0.9*10^6, by=100000),
                 seq(from=10^6, to=50*10^6, by=10^7))
myf <- cut(qtls$dist, my_breaks)
  
alphas <- seq(0.05,0.1,length=5)
filter_thresholds <- c(10^4, 2*10^5, 10^6)
dfs <- list()


for (i in seq_along(alphas)){
  alpha <- alphas[i]
  print(paste0("alpha:",alpha))
  x <- ddhw(qtls$pvalue, myf, alpha, lambdas=Inf,quiet=F)
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
  df$bh_rejections <- sum(p.adjust(qtls$pvalue, method="BH") <= alpha, na.rm=T)
  for (filter_t in filter_thresholds){
    filter_f <- tail(which( my_breaks <= filter_t),1)-1
    filt_pvalue <- qtls$pvalue[as.numeric(myf) <= filter_f]
    df[paste0("threshold:", filter_t)] <- sum(p.adjust(filt_pvalue, method="BH") <= alpha, na.rm=T)
  }
  dfs[[i]] <- df
}

res <- list(alpha_df = bind_rows(dfs),
            breaks   = my_breaks,
            break_min = 5000)

saveRDS(res, file = paste0(projectFolder, "hQTL_benchmark.Rds"))

