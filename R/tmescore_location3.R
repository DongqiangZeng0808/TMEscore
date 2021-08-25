




#' Location of TMEscore of each patient
#'
#' @param score  score
#' @param vars default is TMEscore, TMEscroeA and TMEscoreB
#' @param palette palette of response
#' @param showplot default is FALSE
#' @param path path to save result
#' @param palette_line palette of line
#' @param panel "PFS" or "ORR"
#' @param x_label default is 2.5
#'
#' @return
#' @export
#'
#' @examples
tmescore_location3<-function(score, vars = c("TMEscore", "TMEscoreA", "TMEscoreB"), panel = "PFS", palette = "nrc", palette_line = "jama", showplot = TRUE, path = NULL, tmescore_x = 2.5, tmescore_ab_x = 10){


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
  data(ref_score3)
  ###############################
  pats<-as.character(score$ID)

  for (i in 1:length(pats)) {

    pat<-pats[i]
    print(paste0(">>> Processing patient: ", pat))


    for (j in 1:length(vars)) {

      var<-vars[j]

      if(panel == "ORR"){
        if(var=="TMEscore"){
          cutoff_mono<-9.99
          cutoff_com<-10.19
          cutoff_all<-9.99
          pat_score<-score[score$ID==pat,]$TMEscore
        }else if(var == "TMEscoreA"){
          cutoff_mono<-3.02
          cutoff_com<-3.61
          cutoff_all<-3.21
          pat_score<-score[score$ID==pat,]$TMEscoreA
        }else if(var == "TMEscoreB"){
          cutoff_mono<-3.14
          cutoff_com<- 1.46
          cutoff_all<-3.01
          pat_score<-score[score$ID==pat, ]$TMEscoreB
        }
      }else if(panel == "PFS"){

        if(var=="TMEscore"){
          cutoff_mono<-10.20
          cutoff_com<-10.29
          cutoff_all<-10.29
          pat_score<-score[score$ID==pat,]$TMEscore
        }else if(var == "TMEscoreA"){
          cutoff_mono<-4.49
          cutoff_com<-3.24
          cutoff_all<-3.24
          pat_score<-score[score$ID==pat,]$TMEscoreA
        }else if(var == "TMEscoreB"){
          cutoff_mono<-3.62
          cutoff_com<- 1.03
          cutoff_all<-1.03
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

      pat_split<-unlist(stringr::str_split(pat, pattern = "_"))
      subt<-paste0("Name: ", pat_split[1], "; H: ", pat_split[2], "; Stage: ", pat_split[3],
                   "; Tissue: ", pat_split[4], "; Batch: ", sub(pat_split[5], pattern = "batch", replacement = ""),
                   "; Regimen: ", pat_split[6], "; BOR: ", pat_split[7])

      p<-ggplot(ref_score3, aes(x= !!target,fill= BOR)) +
        geom_histogram(aes(y=..density..), binwidth=.5, colour = "black")+
        scale_fill_manual(values= cols)+
        geom_density(alpha=.2, fill="grey", weight = 2)+
        labs(title=  paste0(target, " = ", pat_score),
             subtitle= paste0(subt),
             caption = paste0(" Data of qPCR: ",panel, ";  ","BC: best cutoff;   ", "mono: monotherapy;   ", "com: combination;   ", date()))+

        xlab(paste0(target))+
        theme_light()+
        IOBR:: design_mytheme(legend.position = "bottom",axis_angle = 0,plot_title_size = 1.5)


        if(var == "TMEscore"){
          p<-p+geom_vline(aes(xintercept = cutoff_all),
                          linetype="dashed",color = cols2[1], size = 0.70)+
            annotate(geom = "text", fontface = "plain", color= cols2[1],
                     x = tmescore_x, y=1.05,hjust = 0,
                     label = paste0('BC of all = ', cutoff_all), size=4.5)+

            geom_vline(aes(xintercept = cutoff_mono),
                       linetype="dashed",color =  cols2[2], size = 0.70)+
            annotate(geom = "text", fontface = "plain", color= cols2[2],
                     x = tmescore_x, y= 0.9,hjust = 0,
                     label = paste0('BC of mono = ', cutoff_mono), size=4.5)+

            geom_vline(aes(xintercept = cutoff_com),
                       linetype="dashed",color = cols2[3], size = 0.70)+
            annotate(geom = "text", fontface = "plain", color= cols2[3],
                     x = tmescore_x, y=0.76,hjust = 0,
                     label = paste0('BC of com = ', cutoff_com), size=4.5)+
            geom_vline(aes(xintercept = pat_score),
                       linetype="dashed",color = "black", size = 0.70)+
            annotate(geom = "text", fontface = "plain", color= "black",
                     x = tmescore_x, y = 1.35,hjust = 0,
                     label = paste0( var, ' of patient = ', pat_score), size=4.5)
        }else{

          p<-p+
          geom_vline(aes(xintercept = cutoff_all),
                        linetype="dashed",color = cols2[1], size = 0.70)+
          annotate(geom = "text", fontface = "plain", color= cols2[1],
                   x = tmescore_ab_x, y = 1, hjust = 0,
                   label = paste0('Best cutoff of all = ', cutoff_all), size=4.5)+

          geom_vline(aes(xintercept = cutoff_mono),
                     linetype="dashed",color =  cols2[2], size = 0.70)+
          annotate(geom = "text", fontface = "plain", color= cols2[2],
                   x = tmescore_ab_x, y=0.80, hjust = 0,
                   label = paste0('Best cutoff of mono = ', cutoff_mono), size=4.5)+

          geom_vline(aes(xintercept = cutoff_com),
                     linetype="dashed",color = cols2[3], size = 0.70)+
          annotate(geom = "text", fontface = "plain", color= cols2[3],
                   x = tmescore_ab_x, y=0.6,hjust = 0,
                   label = paste0('Best cutoff of com = ', cutoff_com), size=4.5)+
          geom_vline(aes(xintercept = pat_score),
                     linetype="dashed",color = "black", size = 0.70)+
          annotate(geom = "text", fontface = "plain", color= "black",
                   x = tmescore_ab_x, y=1.2,hjust = 0,
                   label = paste0( var, ' of patient = ', pat_score),size = 4.5)
        }


      if(showplot) print(p)
      ggsave(p,filename =paste0(i,"-",j,"-",pat,"-",var,".png"),
             width = 7.64,height = 5.76, path = path, dpi = 300)

    }


  }

}
