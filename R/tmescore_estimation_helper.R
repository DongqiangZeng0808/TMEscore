








#' Calculating TMEscore using PCR data
#'
#' @param pdata phenotype data of input sample;
#' if phenotype data is NULL, create a data frame with `Index` and `ID` contain column names of eset
#' @param eset normalizaed  transcriptomic data: normalized (CPM, TPM, RPKM, FPKM, etc.)
#' @param signature List of gene signatures
#' @param mini_gene_count filter out signatures with genes less than minimal gene in expression set;
#' @param column_of_sample  Defines in which column of pdata the sample identifier can be found
#' @param adjust_eset remove variables with missing value, sd =0, and Inf value
#' @param scale default is FALSE
#' @param method default is TRUE
#' @param log2trans
#'
#' @author Dongqiang Zeng
#' @return data frame with pdata and signature scores for gene sets; signatures in columns, samples in rows
#' @export
#'
#' @examples
#'
tmescore_estimation_helper<-function(pdata = NULL,
                                     eset,
                                     signature,
                                     mini_gene_count,
                                     column_of_sample,
                                     scale = FALSE,
                                     method = "mean",
                                     adjust_eset = TRUE,
                                     log2trans = FALSE){

  message(paste0("\n", ">>> Calculating signature score using mean value of signature genes"))

  #creat pdata if NULL
  if(is.null(pdata)){
    pdata<-data.frame("Index" = 1:length(colnames(eset)),"ID" = colnames(eset))
  }else{
    pdata<-as.data.frame(pdata)

    if("ID"%in%colnames(pdata) & !column_of_sample=="ID"){
      colnames(pdata)[which(colnames(pdata)== "ID")]<-"ID2"
      message("In order to prevent duplicate names, the 'ID' column of original pdata was rename into 'ID2' ")
    }

    if(column_of_sample%in%colnames(pdata)){
      colnames(pdata)[which(colnames(pdata)==column_of_sample)]<-"ID"
    }

  }
  #match phenotype data and gene expression set
  ###########################
  pdata<-pdata[pdata$ID%in%colnames(eset),]
  eset<-eset[,colnames(eset)%in%pdata$ID]
  eset<-eset[,match(pdata$ID,colnames(eset))]
  ###########################
  #normalization
  if(log2trans){
    eset<- IOBR::log2eset(eset = eset)
    if(ncol(eset) <5000) IOBR::check_eset(eset)
  }


  if(adjust_eset){
    feas<-IOBR::feature_manipulation(data=eset,is_matrix = T)
    eset<-eset[rownames(eset)%in%feas,]
  }

  if(scale) eset<-scale(eset,center = T,scale = T)
  ###########################
  ###########################
  if(mini_gene_count<=2) mini_gene_count <- 2
  signature<-signature[lapply(signature,function(x) sum(x%in%rownames(eset)==TRUE))>= mini_gene_count]
  ###########################
  #calculating signature score
  goi <- names(signature)
  for (sig in goi) {
    pdata[, sig] <- NA
    genes <- signature[[sig]]
    genes <- genes[genes %in% rownames(eset)]
    tmp <- eset[genes, , drop=FALSE]

    if(method=="mean"){
      pdata[, sig] <- colMeans(tmp)
    }else{
      pdata[, sig] <- colSums(tmp)
    }

  }

  if ("TMEscoreA"%in%goi &"TMEscoreB" %in% goi) {
    pdata[,"TMEscore"]<-pdata[,"TMEscoreA"]-pdata[,"TMEscoreB"]
  }
  pdata<-tibble::as_tibble(pdata)
  return(pdata)
}
###################################################








##' Calculate siganture score using PCA function
##'
##'
##' @param eset normalized count matrix; rows are all genes in the signature that shall be summarized into one score; columns are samples
##' @param methods character vector defining whether signature scores shall be based on principal component 1 ("PC", default) or z-scores (other value)
##'
##' @return numeric vector of length ncol(eset); a score summarizing the rows of eset
##' @author Dorothee Nickles, Dongqiang Zeng
##' @export
sigScore <- function(eset, methods = "PCA") {

  if (methods == "PCA") {
    # message(paste0("Calculating siganture score using PCA function"))
    pc <- prcomp(t(eset))
    sigs <- pc$x[,1] * sign(cor(pc$x[,1], colMeans(eset)))
  } else {
    # message(paste0("Calculating siganture score using mean of signature genes"))
    sigs <- colSums(eset)
  }
  return(sigs)
}
#####################################
