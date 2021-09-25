




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
  score$method<-method
  score$cutoff_ORR<- cutoff[cutoff$method==method&cutoff$regimen=="All"&cutoff$type=="ORR"&cutoff$name=="TMEscore",]$value
  score$cutoff_PFS<- cutoff[cutoff$method==method&cutoff$regimen=="All"&cutoff$type=="PFS"&cutoff$name=="TMEscore",]$value
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
  return(score)
}
