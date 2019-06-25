



#' calculate TMEscore
#'
#' @param eset normalized expression data use function: scale(); with human gene symbol in rows and patients ID in columns;
#' @param pdata phenotype data with patients' ID or other clinical data in columns;the number of rows in pdata should be equal to the number of columns in eset
#' @param method  use PCA(default) or z-score function to calculate TMEscore
#'
#' @return pdata with TMEscoreA, TMEscoreB and TMEscore
#' @export
#' @author Dongqiang Zeng
#' @examples
#' library('tidyverse')
#' # system.file("extdata","example.RData",package = "TMEscore")
#' tmescore<-tmescore(eset = eset,pdata = pdata,method = "pca")
#'

tmescore<-function(eset = eset,pdata = pdata, method = method){

  #check genes in expressionset
  feature_ids <- rownames(eset)
  signature_genes<-c(unlist(signature$TMEscoreA),unlist(signature$TMEscoreB))
  if (! all(feature_ids %in% signature_genes)){
    tbl <- table(feature_ids %in% signature_genes)
    msg1 <- sprintf("%i gene is shared, %i gene is specified", tbl[[2]],tbl[[1]])
    warning(msg1)
  }

  ###############################
  if(!dim(eset)[2]==dim(pdata)[1])
    msg2<-print("expression set data(eset) and phenotype data(pdata) did not have the same length, please modify")

  ###############################
  source("gsScore.R")
  goi <- names(signature);goi
  for (sig in goi) {
    pdata[, sig] <- NA
    genes <- signature[[sig]]
    genes <- genes[genes %in% rownames(eset)]
    tmp <- eset[genes, , drop=FALSE]
    pdata[, sig] <- gsScore(tmp,summarizationFunction = method)
  }
  ####################
  #' TMEscore = TMEscoreA - TMEscoreB
  pdata$TMEscore<-pdata$TMEscoreA-pdata$TMEscoreB
  #####################

}







