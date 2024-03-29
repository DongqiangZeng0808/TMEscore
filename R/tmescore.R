




#' Extract best cut-off and add new binary object to data frame
#'
#' @param pdata phenotype data with surcival information and features
#' @param variable selected featrues that will be classified into two groups
#' @param time column name of survival time
#' @param status column name of event of follow up
#' @param PrintResult logit, print results of survival analysis before and after classifying.
#'
#' @return pdata with binary variables
#' @export
#'
#' @examples
#'
best_cutoff<-function(pdata,variable,time = "time",status = "status",PrintResult = T){

  colnames(pdata)[which(colnames(pdata)==time)]<-"time"
  colnames(pdata)[which(colnames(pdata)==status)]<-"status"
  y<-Surv(pdata$time,pdata$status)
  iscutoff<-surv_cutpoint(pdata,time = "time",event = "status",variables =variable)
  aa<-paste("best cutoff = ",iscutoff$cutpoint$cutpoint)
  plot(iscutoff,variable,palette="npg")
  bb<-summary(coxph(y~pdata[,which(colnames(pdata)==variable)],data =pdata))
  variable2<-paste(variable,"_binary",sep = "")
  pdata[,variable2]<-ifelse(pdata[,variable]<iscutoff$cutpoint$cutpoint,"Low","High")
  pdata[,variable2]<-as.factor(pdata[,variable2])
  cc<-summary(pdata[,variable2])
  dd<-summary(coxph(y~pdata[,which(colnames(pdata)==variable2)],data =pdata))
  if(PrintResult) {
    print(list(best_cutoff=aa,cox_continuous_object=bb,summary_binary_variable = cc,cox_binary_object=dd))
  }
  return(pdata)
}





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
    sigs <- colMeans(eset)
  }
  return(sigs)
}




#' Calculate TMEscore
#'
#' @param eset normalized expression data use function: scale(); with human gene symbol in rows and patient identifier in columns;
#' @param pdata phenotype data with patient identifier or other clinical data in columns;
#' the number of rows in pdata should be equal to the number of columns in eset
#' @param method  use PCA (default) or z-score function to calculate TMEscore
#' @param print_gene_pro print signature gene propotion, default is FALSE
#' @param column_of_sample Defines in which column of pdata the sample identifier can be found.
#' @param classify classifying TMEscore into two groups using survival data and function of suvmmer.
#'
#' @return pdata with TMEscoreA, TMEscoreB and TMEscore
#' @export
#' @author Dongqiang Zeng
#' @examples
#'
#' tmescore<-tmescore(eset = eset_stad,pdata = pdata_stad,column_of_sample = "ID")
#'

tmescore<-function(eset, pdata = NULL, column_of_sample = "ID",method = "PCA",print_gene_pro = FALSE, classify = FALSE){


  #to prevent crashing on duplicated gene symbols, add unique numbers to identical names
  dups <- dim(eset)[1] - length(unique(rownames(eset)))
  if(dups > 0) {
    warning(paste(dups," duplicated gene symbol(s) found in mixture file! \n
                  Expression set will be order by the mean value, and genes with the highest mean value will be preserved",sep=""))
    eset<-eset[order(apply(eset,1,mean), decreasing = T),]
    eset<-eset[!duplicated(rownames(eset)),]
  }
  ###################################
  data(signature)
  #check genes in expression set
  ###################################
  freq1<-length(intersect(signature$TMEscoreA,rownames(eset)))/length(signature$TMEscoreA)
  if(freq1<0.5){
    msg1<- paste0(paste0(sprintf(" Only %1.2f%%", 100*freq1)," of TMEscoreA signature genes appear on gene matrix,\n interpret results with caution"))
    warning(msg1)
  }else if(freq1>=0.5 & print_gene_pro){
    message(paste0(paste0(sprintf("  %1.2f%%", 100*freq1)," of TMEscoreA signature genes appear on gene matrix")))
  }
  ##################################
  freq2<-length(intersect(signature$TMEscoreB,rownames(eset)))/length(signature$TMEscoreB)
  if(freq2<0.5){
    msg1<- paste0(paste0(sprintf(" Only %1.2f%%", 100*freq2)," of TMEscoreB signature genes appear on gene matrix,\n interpret results with caution"))
    warning(msg1)
  }else if(freq2>=0.5 & print_gene_pro){
    message(paste0(paste0(sprintf("  %1.2f%%", 100*freq2)," of TMEscoreB signature genes appear on gene matrix")))
  }
  ##################################

  #filter signatures
  mini_gene_count <- 2 # count of overlapping genes should more than 2
  tme_signature<-signature[c("TMEscoreA","TMEscoreB")]
  tme_signature<-tme_signature[lapply(tme_signature,function(x) sum(x%in%rownames(eset)==TRUE))>= mini_gene_count]
  #################################


  if(!"TMEscoreA"%in% names(tme_signature) & !"TMEscoreB"%in% names(tme_signature)){
    stop("TMEscore signature genes did not appear on expression matirx ")
  } else if("TMEscoreA"%in% names(tme_signature) & !"TMEscoreB"%in% names(tme_signature)){
    warning("Only TMEscoreA will be estimated for insufficient TMEscoreB signature genes appear on gene matrix")
  } else if(!"TMEscoreA"%in% names(tme_signature) & "TMEscoreB"%in% names(tme_signature)){
    warning("Only TMEscoreB will be estimated for insufficient TMEscoreA signature genes appear on gene matrix")
  }
  #################################
  #creat pdata if not provided
  if(is.null(pdata)){
    pdata<-data.frame("Index" = 1:length(colnames(eset)),"ID" = colnames(eset))
    ###############################
    if(!dim(eset)[2]==dim(pdata)[1])
      message("Expression set data (eset) and phenotype data (pdata) have different number of samples,
            \n only samples with expression set provided will be estimated")
  }

  #match phenotype data and gene expression set
  ###########################
  colnames(pdata)[which(colnames(pdata)==column_of_sample)]<-"ID"
  pdata<-pdata[pdata$ID%in%colnames(eset),]
  eset<-eset[,colnames(eset)%in%pdata$ID]
  eset<-eset[,match(pdata$ID,colnames(eset))]
  ###############################

  # if(max(eset)<50 & array) eset<- preprocessCore::normalize.quantiles(eset)

  if(max(eset)>100) eset<-log2(eset+1)
  ###############################
  goi <- names(tme_signature)
  for (sig in goi) {
    pdata[, sig] <- NA
    genes <- tme_signature[[sig]]
    genes <- genes[genes %in% rownames(eset)]
    tmp <- eset[genes, , drop=FALSE]
    pdata[, sig] <- sigScore(tmp,method = method)
  }
  ####################
  if ("TMEscoreA"%in%goi &"TMEscoreB" %in% goi) {
    pdata[,"TMEscore"]<-pdata[,"TMEscoreA"]-pdata[,"TMEscoreB"]
  }
  #####################
  if(classify & !c("time"%in%colnames(pdata) & "status"%in%colnames(pdata)) ) {
   message("If survival data are provided, please transform column names of follow up data into 'time' and 'status' ")
  }

  if(classify & "time"%in%colnames(pdata) &"status"%in%colnames(pdata)) {
    pdata<-best_cutoff(pdata = pdata, variable = "TMEscore",time = "time",status = "status",PrintResult = F)
  }

  return(pdata)
}

