






##' Calculate score across genes and samples
##'
##' This wrapper function combines filtering out genes with low reads in a number of samples (recommended for limma:voom) with normalization
##' @param gm normalized count matrix; rows are all genes in the signature that shall be summarized into one score; columns are samples
##' @param summarizationFunction character vector defining whether signature scores shall be based on principal component 1 ("PC", default) or z-scores (other value)
##' @return numeric vector of length ncol(gm); a score summarizing the rows of gm
##' @author Dorothee Nickles
##' @export
gsScore <- function(gm, summarizationFunction="pca") {
  if (summarizationFunction == "pca") {
    pc <- prcomp(t(gm),
                 retx=TRUE)
    gss <- pc$x[,1] * sign(cor(pc$x[,1], colMeans(gm)))
  } else {
    gss <- colMeans(gm)
  }
  return(gss)
}



#' Calculate TMEscore
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
#' tmescore<-tmescore(eset = eset_acrg,pdata = pdata,method = "pca")
#'

tmescore<-function(eset, pdata, method = "pca"){

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
  return(pdata)
}

