



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
#' @param log2trans default is FALSE
#' @param replace_na default is FALSE
#' @param save_data default is FALSE
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
                                     mini_gene_count = 2,
                                     column_of_sample = "ID",
                                     replace_na = FALSE,
                                     scale = FALSE,
                                     method = "mean",
                                     adjust_eset = FALSE,
                                     log2trans = FALSE,
                                     save_data = FALSE){

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

  if(sum(is.na(eset)>0)) message(">>> Parameter `adjust_eset` must be FALSE, if variables with NA wanted to be preserved")

  if(adjust_eset){
    feas<-IOBR::feature_manipulation(data=eset,is_matrix = T)
    eset<-eset[rownames(eset)%in%feas,]
  }


  if(log2trans&replace_na==FALSE) stop("If paramater `log2trans` is TRUE, `replace_na` must be TRUE to replace NA")

  if(replace_na){

    if(sum(is.na(eset)>0)){

      message(">>> Retain NA variables, replaced by mean value of all observations")
      teset<-as.data.frame(t(eset))
      for(i in c(1:ncol(teset))){
        teset[is.na(teset[,i]),i]<-mean(teset[!is.na(teset[,i]),i])
      }
      eset<-t(teset)

      if(log2trans){

        feas<-IOBR::feature_manipulation(data=eset,is_matrix = T)
        eset<-eset[rownames(eset)%in%feas,]

        eset<- IOBR::log2eset(eset = eset)
        if(ncol(eset) <5000) IOBR::check_eset(eset)
      }

    }else{
      message(">>> There are no missing values")
    }
  }

  if(scale){
    if(sum(is.na(eset))>0){
      message("Before scaling, NA will be repleased by mean value")

      teset<-as.data.frame(t(eset))
      for(i in c(1:ncol(teset))){
        teset[is.na(teset[,i]),i]<-mean(teset[!is.na(teset[,i]),i])
      }
      eset<-t(teset)
      message(paste0("Counts of NA after replacement = "), sum(is.na(eset)))

      feas<-IOBR::feature_manipulation(data=eset,is_matrix = T)
      eset<-eset[rownames(eset)%in%feas,]

    }
    eset<-scale(eset, center = T,scale = T)
  }
  ###########################
  if(save_data) write.csv(eset, "eset.csv")
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
      pdata[, sig] <- colMeans(tmp,na.rm = T)
    }else{
      pdata[, sig] <- colSums(tmp, na.rm = T)
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
    sigs <- pc$x[,1] * sign(cor(pc$x[,1], colMeans(eset, na.rm = T)))
  } else {
    # message(paste0("Calculating siganture score using mean of signature genes"))
    sigs <- colSums(eset, na.rm = T)
  }
  return(sigs)
}
#####################################
