




#' Title
#'
#' @param score
#' @param method
#'
#' @return
#' @export
#'
#' @examples
tmescore_report<-function(score, method){

  data(cutoff)

  # res<-data.frame(NULL)
  #
  # if(nrow(score)>1){
  #
  #   for (i in 1:nrow(res)) {
  #
  #   }
  #
  # }


  score$method<-method
  score$cutoff_ORR<- cutoff[cutoff$method==method&cutoff$regimen=="All"&cutoff$type=="ORR"&cutoff$name=="TMEscore",]$value
  score$cutoff_PFS<- cutoff[cutoff$method==method&cutoff$regimen=="All"&cutoff$type=="PFS"&cutoff$name=="TMEscore",]$value

  score$cutoff_TMEscoreA<- cutoff[cutoff$method==method&cutoff$regimen=="All"&cutoff$type=="ORR"&cutoff$name=="TMEscoreA",]$value
  score$cutoff_TMEscoreB<- cutoff[cutoff$method==method&cutoff$regimen=="All"&cutoff$type=="ORR"&cutoff$name=="TMEscoreB",]$value

  ##############################
  if(score$TMEscore < score$cutoff_ORR){
    score$benefit_ORR<-"No"
  }else{
    score$benefit_ORR<-"Yes"
  }
  if(score$TMEscore < score$cutoff_PFS){
    score$benefit_PFS<-"No"
  }else{
    score$benefit_PFS<-"Yes"
  }


  if(score$TMEscoreA < score$cutoff_TMEscoreA){
    score$TMEscoreA_location<-"low_TMEscoreA"
  }else{
    score$TMEscoreA_location<-"high_TMEscoreA"
  }

  if(score$TMEscoreB < score$cutoff_TMEscoreB){
    score$TMEscoreB_location<-"low_TMEscoreB"
  }else{
    score$TMEscoreB_location<-"high_TMEscoreB"
  }


  if(score$TMEscoreA < score$cutoff_TMEscoreA & score$TMEscoreB < score$cutoff_TMEscoreB ){
    score$prediction<-"Unreliable_prediction"
  }else{
    score$prediction<-"Reliable_prediction"
  }
  ##############################

  return(score)
}
