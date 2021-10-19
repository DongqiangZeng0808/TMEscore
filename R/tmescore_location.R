









#' Location of TMEscore of each patient
#'
#' @param score  score
#' @param vars default is TMEscore, TMEscroeA and TMEscoreB
#' @param palette palette of response
#' @param showplot default is FALSE
#' @param path path to save result
#' @param palette_line palette of line
#' @param panel "PFS" or "ORR"
#' @param tmescore_x 7
#' @param tmescore_ab_x default is 3.8
#' @param tmescoreab_y_index
#' @param method default is method2
#' @param tmescore_y_index
#' @param fig.type default is png
#'
#' @return
#' @export
#'
#' @examples
tmescore_location<-function(score,
                            vars = c("TMEscore", "TMEscoreA", "TMEscoreB"),
                            panel = "ORR",
                            method = "method2",
                            palette = "nrc",
                            palette_line = "jama",
                            showplot = TRUE,
                            path = NULL,
                            tmescore_x = 7,
                            tmescore_ab_x = 2,
                            tmescore_y_index = -0.5,
                            tmescoreab_y_index = 0.5,
                            fig.type = "png"){


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
  if(method=="method1"){
    data(ref_score1)
    ref_score<-ref_score1
    binwidth = 0.45
  }else if(method == "method2"){
    data(ref_score2)
    ref_score<-ref_score2
    binwidth = 0.15
  }else if(method == "method3"){
    data(ref_score3)
    ref_score<-ref_score3
    binwidth = 0.85
  }else if(method == "method4"){
    data(ref_score4)
    ref_score<-ref_score4
    binwidth = 0.55
  }else if(method == "method5"){
    data(ref_score5)
    ref_score<-ref_score5
    binwidth = 0.5

  }else if(method == "method6"){
    data(ref_score6)
    ref_score<-ref_score6
    binwidth = 5.5

  }else if(method == "method7"){
    data(ref_score7)
    ref_score<-ref_score7
    binwidth = 0.85

  }

  ###############################
  pats<-as.character(score$ID)

  for (i in 1:length(pats)) {

    pat<-pats[i]
    print(paste0(">>> Processing patient: ", pat))


    for (j in 1:length(vars)) {

      var<-vars[j]

      data("cutoff")
      cutoff_mono<-cutoff[cutoff$method==method&cutoff$regimen=="monotherapy"&cutoff$type==panel&cutoff$name==var,]$value
      cutoff_com<-cutoff[cutoff$method==method&cutoff$regimen=="combination"&cutoff$type==panel&cutoff$name==var,]$value
      cutoff_all<-cutoff[cutoff$method==method&cutoff$regimen=="All"&cutoff$type==panel&cutoff$name==var,]$value
      pat_score<-as.numeric(score[score$ID==pat,colnames(score)==var])

      pat_score<-round(pat_score, 2)
      message(paste0(">>> ", var, " of ", pat, " is ", pat_score))
      target<-sym(var)
      cols<-IOBR::palettes(category = "box", palette = palette,
                           show_message = FALSE, show_col = FALSE, alpha = 0.75)

      cols2<-IOBR::palettes(category = "box", palette = palette_line,
                            show_message = FALSE, show_col = FALSE, alpha = 1)

      pat_split<-unlist(stringr::str_split(pat, pattern = "_"))
      subt<-paste0("Name: ", pat_split[1], " | H: ", pat_split[2], " | C: ", pat_split[3],
                   " | T: ", pat_split[4], " | B: ", sub(pat_split[5], pattern = "batch", replacement = ""),
                   " | R: ", pat_split[6], " | BOR: ", pat_split[7])




      if(var == "TMEscore"){
        p<-ggplot(ref_score, aes(x= !!target,fill= BOR)) +
          geom_histogram(aes(y=..density..), binwidth = binwidth, colour = "black")+
          scale_fill_manual(values= cols)+
          geom_density(alpha=.2, fill="grey", weight = 2)+
          labs(title=  paste0(target, " = ", pat_score),
               subtitle= paste0(subt),
               caption = paste0(" Cutoff of: ",panel, ";  ","BC: best cutoff;   ", "mono: monotherapy;   ",
                                "com: combination;   ", date()))+


          theme_light()+
          IOBR:: design_mytheme(legend.position = "bottom",axis_angle = 0,plot_title_size = 2)+
          xlab(NULL)

        p<-p+geom_vline(aes(xintercept = cutoff_all),
                        linetype="dashed",color = cols2[1], size = 0.70)+
          annotate(geom = "text", fontface = "plain", color= cols2[1],
                   x = tmescore_x, y=0.9+tmescore_y_index, hjust = 0,
                   label = paste0('BC of all = ', cutoff_all), size=4.5)+

          geom_vline(aes(xintercept = cutoff_mono),
                     linetype="dashed",color =  cols2[2], size = 0.70)+
          annotate(geom = "text", fontface = "plain", color= cols2[2],
                   x = tmescore_x, y= 0.8+tmescore_y_index,hjust = 0,
                   label = paste0('BC of mono = ', cutoff_mono), size=4.5)+

          geom_vline(aes(xintercept = cutoff_com),
                     linetype="dashed",color = cols2[3], size = 0.70)+
          annotate(geom = "text", fontface = "plain", color= cols2[3],
                   x = tmescore_x, y=0.7+tmescore_y_index,hjust = 0,
                   label = paste0('BC of com = ', cutoff_com), size=4.5)+
          geom_vline(aes(xintercept = pat_score),
                     linetype="dashed",color = "black", size = 0.70)+
          annotate(geom = "text", fontface = "plain", color= "black",
                   x = tmescore_x, y = 1.05+tmescore_y_index,hjust = 0,
                   label = paste0( var, ' of patient = ', pat_score), size=4.5)
      }else{

        p<-ggplot(ref_score, aes(x= !!target,fill= BOR)) +
          geom_histogram(aes(y=..density..), binwidth = 0.15, colour = "black")+
          scale_fill_manual(values= cols)+
          geom_density(alpha=.2, fill="grey", weight = 2)+
          labs(title=  paste0(target, " = ", pat_score),
               subtitle= paste0(subt),
               caption = paste0(" Cutoff of: ",panel, ";  ","BC: best cutoff;   ", "mono: monotherapy;   ",
                                "com: combination;   ", date()))+


          theme_light()+
          IOBR:: design_mytheme(legend.position = "bottom",axis_angle = 0,plot_title_size = 2)+
          xlab(NULL)

        p<-p+
          geom_vline(aes(xintercept = cutoff_all),
                     linetype="dashed",color = cols2[1], size = 0.70)+
          annotate(geom = "text", fontface = "plain", color= cols2[1],
                   x = tmescore_ab_x, y = 1.1+tmescoreab_y_index, hjust = 0,
                   label = paste0('Best cutoff of all = ', cutoff_all), size=3)+

          geom_vline(aes(xintercept = cutoff_mono),
                     linetype="dashed",color =  cols2[2], size = 0.70)+
          annotate(geom = "text", fontface = "plain", color= cols2[2],
                   x = tmescore_ab_x, y=1+tmescoreab_y_index, hjust = 0,
                   label = paste0('Best cutoff of mono = ', cutoff_mono), size=3)+

          geom_vline(aes(xintercept = cutoff_com),
                     linetype="dashed",color = cols2[3], size = 0.70)+
          annotate(geom = "text", fontface = "plain", color= cols2[3],
                   x = tmescore_ab_x, y=0.9+tmescoreab_y_index,hjust = 0,
                   label = paste0('Best cutoff of com = ', cutoff_com), size=3)+
          geom_vline(aes(xintercept = pat_score),
                     linetype="dashed",color = "black", size = 0.70)+
          annotate(geom = "text", fontface = "plain", color= "black",
                   x = tmescore_ab_x, y=1.25+tmescoreab_y_index,hjust = 0,
                   label = paste0( var, ' of patient = ', pat_score),size = 3)
      }


      if(showplot) print(p)
      ggplot2:: ggsave(plot = p, filename = paste0(i,"-",j,"-",pat,"-",var,".",fig.type),
             width = 7.64, height = 5.76, path = path)

    }

  }

}
