



#' TMEscore calculation
#'
#' @param eset
#' @param scale
#' @param method
#' @param mini_gene_count
#' @param print_gene_pro
#' @param save_eset
#' @param file_name
#' @param sig
#' @param max
#'
#' @return
#' @export
#' @author Dongqiang Zeng
#' @examples
tmescore_pcr<-function(eset, scale = T, sig = "pcr", method = "mean", max = 20, mini_gene_count = 2, print_gene_pro = T, save_eset = F, file_name = NULL){


  eset<-as.data.frame(eset)
  # eset<-matrix(as.numeric(eset), dim(eset), dimnames = dimnames(eset))
  # eset<-eset[complete.cases(eset),]



  if(max(eset, na.rm = T)>50){
    message(paste0(">>> The maximum value of Fc is = ",max(eset, na.rm = T), ". Data correction is preferred!!" ))
    message(">>> The maximum value can be limited by the parameter `max`.")
  }else{
    message(paste0(">>> The maximum value of Fc is = ",max(eset, na.rm = T)))
    message(paste0(">>> The minimum value of Fc is = ",min(eset, na.rm = T)))
  }

  if(!is.null(max)){
    eset[eset>max]<-max
    message(paste0(">>> The maximum value was limited to:: ",max, " by the parameter `max`."))
  }

  eset10<-eset
  #################################

  if(sig=="pcr"){
    tme_signature<-signature[c("TMEscoreA_pcr","TMEscoreB_pcr")]
    names(tme_signature)<-c("TMEscoreA","TMEscoreB")

    message("The modified TME_signature was chosen to estimate TMEscore: for qPCR data.")

  }else if(sig=="jitc"){
    tme_signature<-signature[c("TMEscoreA","TMEscoreB")]

    message("The JITC TME_signature was chosen to estimate TMEscore: for RNAseq or Array data.")

  }else if(sig=="cir"){
    tme_signature<-signature[c("TMEscoreA_CIR","TMEscoreB_CIR")]
    names(tme_signature)<-c("TMEscoreA","TMEscoreB")
    message("The CIR TME_signature was chosen to estimate TMEscore: for RNAseq or Array data.")
  }
  ################################

  freq1<-length(intersect(tme_signature$TMEscoreA,rownames(eset)))/length(tme_signature$TMEscoreA)
  if(freq1<0.5){
    msg1<- paste0(paste0(sprintf(" Only %1.2f%%", 100*freq1)," of TMEscoreA signature genes appear on gene matrix,\n interpret results with caution"))
    warning(msg1)
  }else if(freq1>=0.5 & print_gene_pro){
    message(paste0(paste0(sprintf("  %1.2f%%", 100*freq1)," of TMEscoreA signature genes appear on gene matrix")))
  }
  ##################################
  freq2<-length(intersect(tme_signature$TMEscoreB,rownames(eset)))/length(tme_signature$TMEscoreB)
  if(freq2<0.5){
    msg1<- paste0(paste0(sprintf(" Only %1.2f%%", 100*freq2)," of TMEscoreB signature genes appear on gene matrix,\n interpret results with caution"))
    warning(msg1)
  }else if(freq2>=0.5 & print_gene_pro){
    message(paste0(paste0(sprintf("  %1.2f%%", 100*freq2)," of TMEscoreB signature genes appear on gene matrix")))
  }
  ##################################

  #filter signatures
  mini_gene_count <- 2 # count of overlapping genes should more than 2

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
                                    method = method,
                                    save_data = save_eset,
                                    file_name = file_name)

  if(method =="mean"){
    score$TMEscoreA<-score$TMEscoreA*10
    score$TMEscoreB<-score$TMEscoreB*10
    score$TMEscore<-score$TMEscore*10 + 20
  }else if(method=="PCA"){
    score$TMEscoreA<-score$TMEscoreA
    score$TMEscoreB<-score$TMEscoreB
    score$TMEscore<-score$TMEscore+10
  }


  return(score)

}
