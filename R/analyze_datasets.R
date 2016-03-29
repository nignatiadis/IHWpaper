#' analyze_dataset: Basically performs preprocessing and then returns analyzed RNASeq dataset (diff. expression)
#'           , i.e. the DESeq2 result whose p-values and baseMean statistics can then be used with DDHW
#'
#'
#' @param dataset  Character, name of dataset to be preprocessed, only 4 choices currently available
#' @param res      (default TRUE): return result table, rather than DESeq2 object
#'
#' @return Preprocessed dataset
#' @examples pasilla <- analyze_dataset("pasilla")
#' @export
analyze_dataset <- function(dataset=c("pasilla","airway","bottomly","pasilla"), res= T) {

    if (dataset=="airway") {
      requireNamespace("airway", quietly = TRUE)
      data("airway", package="airway")
      dds <- DESeqDataSet(se = airway, design = ~ cell + dex)
      dds <- DESeq(dds)

    } else if (dataset == "pasilla") {
      requireNamespace("pasilla", quietly=TRUE)
      data("pasillaGenes", package="pasilla")
      countData <- counts(pasillaGenes)
      colData <- pData(pasillaGenes)[,c("condition","type")]
      dds <- DESeqDataSetFromMatrix(countData = countData,
                                    colData = colData,
                                    design = ~ condition)
      dds$condition <- factor(dds$condition,
                              levels=c("untreated","treated"))
      dds <- DESeq(dds)
    } else if (dataset == "bottomly") {
      load(system.file("real_data_examples/raw_data", "bottomly_eset.RData", package = "IHWpaper"))
      countData <- exprs(bottomly.eset)
      colData <- pData(bottomly.eset)
      dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData , design = ~strain)
      dds <- DESeq(dds)

    } else if (dataset == "hammer") {
      load(system.file("real_data_examples/raw_data", "hammer_eset.RData", package = "IHWpaper"))
      countData <- exprs(hammer.eset)
      colData <- pData(hammer.eset)
      dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData , design = ~protocol)
      dds <- dds[, dds$Time == "2 months"]
      dds <- DESeq(dds)
    } else {
      stop("No such dataset currently available!")
    }

  if (res){
    vars <- rowVars(counts(dds,normalized=T))
    resDf <- as.data.frame(results(dds))
    resDf$var <- vars
    resDf$max <- rowMax(counts(dds,normalized=T))
    resDf$min <- rowMin(counts(dds,normalized=T))
    resDf <-subset(resDf, !is.na(pvalue)) #maybe remove that since IHW should be able to handle NAs anyway
    return (resDf)
  } else {
    return (dds)
  }

}