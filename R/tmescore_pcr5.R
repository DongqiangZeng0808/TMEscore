





#' Title
#'
#' @param eset
#' @param scale
#'
#' @return
#' @export
#'
#' @examples
tmescore_pcr5<-function(eset, scale = FALSE){


  if(scale) message("All score estimated by method 1-4 will be scaled......")

  message("<<< Calculating TMEscore using method-1 >>>")
  score1<-tmescore_pcr(eset = dat,save_eset = F)
  colnames(score1)[2:4]<-paste0(colnames(score1)[2:4], "_m1")


  message("                                           ")
  message("<<< Calculating TMEscore using method-2 >>>")
  score2<-tmescore_pcr2(eset = dat,save_eset = F)
  colnames(score2)[2:4]<-paste0(colnames(score2)[2:4], "_m2")
  if(scale) score2[,2:4]<-scale(score2[,2:4], center = T, scale =  T)

  message("                                           ")
  message("<<< Calculating TMEscore using method-3 >>>")
  score3<-tmescore_pcr3(eset = dat,save_eset = F)
  colnames(score3)[2:4]<-paste0(colnames(score3)[2:4], "_m3")
  if(scale) score3[,2:4]<-scale(score3[,2:4], center = T, scale =  T)

  message("                                           ")
  message("<<< Calculating TMEscore using method-4 >>>")
  score4<-tmescore_pcr4(eset = dat,save_eset = F)
  colnames(score4)[2:4]<- paste0(colnames(score4)[2:4], "_m4")
  # print(head(score4))
  if(scale) score4[,2:4]<-scale(score4[,2:4], center = T, scale =  T)


  score <- inner_join(score1, score2, by = "ID") %>%
    inner_join(., score3, by = "ID") %>%
    inner_join(., score4, by = "ID")

  score$TMEscoreA<-base:: rowSums(score[,c("TMEscoreA_m1","TMEscoreA_m2","TMEscoreA_m3","TMEscoreA_m4")])
  score$TMEscoreB<-base::rowSums(score[,c("TMEscoreB_m1","TMEscoreB_m2","TMEscoreB_m3","TMEscoreB_m4")])
  score$TMEscore<-base::rowSums(score[,c("TMEscore_m1","TMEscore_m2","TMEscore_m3","TMEscore_m4")])

  if(scale) score[,c("TMEscore","TMEscoreA","TMEscoreB")]<-scale(score[,c("TMEscore","TMEscoreA","TMEscoreB")], center = T, scale =  T)

  return(score)
}
