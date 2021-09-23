



#' TMEscore calculation
#'
#' @param eset
#' @param scale
#' @param method
#' @param mini_gene_count
#' @param print_gene_pro
#' @param coef
#'
#' @return
#' @export
#' @author Dongqiang Zeng
#' @examples
tmescore_pcr4<-function(eset, scale = F, method = "mean", mini_gene_count = 2, print_gene_pro = T, coef = 1){


  eset<-as.data.frame(eset)
  # eset<-matrix(as.numeric(eset), dim(eset), dimnames = dimnames(eset))
  # eset<-eset[complete.cases(eset),]

  eset10<-eset*coef
  #################################
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

  score<-tmescore_estimation_helper(eset = eset10,
                                    signature = tme_signature,
                                    mini_gene_count = mini_gene_count,
                                    scale = scale,
                                    method = method)


  # score$TMEscoreA<-score$TMEscoreA*10
  # score$TMEscoreB<-score$TMEscoreB*10
  #
  score$TMEscore<- score$TMEscore+20

  return(score)

}
