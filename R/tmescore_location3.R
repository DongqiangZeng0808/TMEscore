




#' Location of TMEscore of each patient
#'
#' @param score  score
#' @param vars default is TMEscore, TMEscroeA and TMEscoreB
#' @param palette palette of response
#' @param showplot default is FALSE
#' @param path path to save result
#' @param palette_line palette of line
#' @param panel "PFS" or "ORR"
#'
#' @return
#' @export
#'
#' @examples
tmescore_location<-function(score, vars = c("TMEscore", "TMEscoreA", "TMEscoreB"), panel = "PFS", palette = "nrc", palette_line = "jama",
                            showplot = TRUE, path = NULL ){


  score<-as.data.frame(score)
  if(!is.null(path)){
    if(!file.exists(path)) dir.create(path)
  }else{
    path<- "TMEscore-location"
    if(!file.exists(path)) dir.create(path)
  }

  if(is.null(cols)){
    cols<-IOBR::palettes(category = "box", palette = palette,
                         show_message = FALSE, show_col = FALSE)
  }else{
    cols<-cols
  }
  ###############################
  data(ref_score)
  ###############################
  pats<-as.character(score$ID)

  for (i in 1:length(pats)) {

    pat<-pats[i]
    print(paste0(">>> Processing patient: ", pat))


    for (j in 1:length(vars)) {

      var<-vars[j]

      if(panel == "ORR"){
        if(var=="TMEscore"){
          cutoff_mono<-6.62
          cutoff_com<-6.15
          cutoff_all<-6.15
          pat_score<-score[score$ID==pat,]$TMEscore
        }else if(var == "TMEscoreA"){
          cutoff_mono<-10.94
          cutoff_com<-11.00
          cutoff_all<-10.94
          pat_score<-score[score$ID==pat,]$TMEscoreA
        }else if(var == "TMEscoreB"){
          cutoff_mono<-10.94
          cutoff_com<- 10.25
          cutoff_all<-10.25
          pat_score<-score[score$ID==pat, ]$TMEscoreB
        }
      }else if(panel == "PFS"){

        if(var=="TMEscore"){
          cutoff_mono<-6.5
          cutoff_com<-6.84
          cutoff_all<-6.92
          pat_score<-score[score$ID==pat,]$TMEscore
        }else if(var == "TMEscoreA"){
          cutoff_mono<-11.64
          cutoff_com<-9.56
          cutoff_all<-10.9
          pat_score<-score[score$ID==pat,]$TMEscoreA
        }else if(var == "TMEscoreB"){
          cutoff_mono<-10.5
          cutoff_com<- 8.3
          cutoff_all<-8.3
          pat_score<-score[score$ID==pat, ]$TMEscoreB
        }
      }



      pat_score<-round(pat_score, 2)
      message(paste0(">>> ", var, " of ", pat, " is ", pat_score))
      target<-sym(var)
      cols<-IOBR::palettes(category = "box", palette = palette,
                           show_message = FALSE, show_col = FALSE, alpha = 0.75)

      cols2<-IOBR::palettes(category = "box", palette = palette_line,
                           show_message = FALSE, show_col = FALSE, alpha = 1)

      p<-ggplot(ref_score, aes(x= !!target,fill= BOR)) +
        geom_histogram(aes(y=..density..), binwidth=.5, colour = "black")+
        scale_fill_manual(values= cols)+
        geom_density(alpha=.2, fill="grey", weight = 2)+
        labs(title=  paste0(target, " of ", pat),
             subtitle= paste0(var, " of patient = ",pat_score ),
             caption = paste0(" Data of qPCR ",panel, ";  ","BC: best cutoff;   ", "mono: monotherapy;   ", "com: combination;   ", date()))+

        xlab(paste0(target))+
        theme_light()+
        IOBR:: design_mytheme(legend.position = "bottom",axis_angle = 0)


        if(var == "TMEscore"){
          p<-p+geom_vline(aes(xintercept = cutoff_all),
                          linetype="dashed",color = cols2[1], size = 0.70)+
            annotate(geom = "text", fontface = "plain", color= cols2[1],
                     x = cutoff_all+0.4, y=0.36,
                     label = paste0('BC of all = ', cutoff_all), size=3.5,angle = 90)+

            geom_vline(aes(xintercept = cutoff_mono),
                       linetype="dashed",color =  cols2[2], size = 0.70)+
            annotate(geom = "text", fontface = "plain", color= cols2[2],
                     x = cutoff_mono-0.4, y=0.35,
                     label = paste0('BC of mono = ', cutoff_mono), size=3.5,angle = 90)+

            geom_vline(aes(xintercept = cutoff_com),
                       linetype="dashed",color = cols2[3], size = 0.70)+
            annotate(geom = "text", fontface = "plain", color= cols2[3],
                     x = cutoff_com-0.4, y=0.35,
                     label = paste0('BC of com = ', cutoff_com), size=3.5,angle = 90)+
            geom_vline(aes(xintercept = pat_score),
                       linetype="dashed",color = "black", size = 0.70)+
            annotate(geom = "text", fontface = "plain", color= "black",
                     x = pat_score-0.4, y=0.33,
                     label = paste0( var, ' of patient = ', pat_score), size=3.5,angle = 90)
        }else{

          p<-p+
          geom_vline(aes(xintercept = cutoff_all),
                        linetype="dashed",color = cols2[1], size = 0.70)+
          annotate(geom = "text", fontface = "plain", color= cols2[1],
                   x = 7, y = 1, hjust = 0,
                   label = paste0('Best cutoff of all = ', cutoff_all), size=3.5)+

          geom_vline(aes(xintercept = cutoff_mono),
                     linetype="dashed",color =  cols2[2], size = 0.70)+
          annotate(geom = "text", fontface = "plain", color= cols2[2],
                   x = 7, y=0.90, hjust = 0,
                   label = paste0('Best cutoff of mono = ', cutoff_mono), size=3.5)+

          geom_vline(aes(xintercept = cutoff_com),
                     linetype="dashed",color = cols2[3], size = 0.70)+
          annotate(geom = "text", fontface = "plain", color= cols2[3],
                   x = 7, y=0.76,hjust = 0,
                   label = paste0('Best cutoff of com = ', cutoff_com), size=3.5)+
          geom_vline(aes(xintercept = pat_score),
                     linetype="dashed",color = "black", size = 0.70)+
          annotate(geom = "text", fontface = "plain", color= "black",
                   x = 7, y=0.62,hjust = 0,
                   label = paste0( var, ' of patient = ', pat_score))
        }


      if(showplot) print(p)
      ggsave(p,filename =paste0(i,"-",j,"-",pat,"-",var,".pdf"),
             width = 6.64,height = 5.76, path = path)

    }


  }

}
