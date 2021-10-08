









#' Title
#'
#' @param score
#' @param pcr_data
#' @param BOR
#' @param index
#' @param jitter
#' @param fig.type
#' @param fig.path
#' @param palette_box
#' @param palette_roc
#' @param palette_his
#' @param binwidth
#' @param addtable
#'
#' @return
#' @export
#'
#' @examples
tmescore_value<-function(score, pcr_data, BOR = "BOR", index = "1-method1", jitter = TRUE,
                      binwidth = 3, addtable = T,
                      palette_box = "jama", palette_roc = "jama",
                      palette_his = 'jco',
                      fig.type = "png", fig.path = "1-testing-score-value" ){


  #' @转置数据
  ################################
  eset<-pcr_data
  eset<-rownames_to_column(as.data.frame(t(eset)),var = "ID")
  summary(eset)
  ##################################

  file_name<- mydb::creat_folder(fig.path)

  score<-as.data.frame(score)

  colnames(score)[which(colnames(score)==BOR)]<-"BOR"

  message(">>>>>> Data frame `score` should contaied variables: ID,TMEscoreA,TMEscoreB,TMEscore,BOR,name,Tissue_type,tissue, batch......")

  input<-score[,c("ID","TMEscoreA","TMEscoreB","TMEscore","BOR","name","Tissue_type","tissue","batch","cohort")]


  input$BOR3<-plyr:: revalue(input$BOR, c("CR" = "CRPR", "PR"= "CRPR","SD"= "SD","PD"="PD"))
  input$BOR2<-plyr:: revalue(input$BOR, c("CR" = "CRPR", "PR"= "CRPR","SD"= "SDPD","PD"="SDPD"))
  input$BOR01<-plyr:: revalue(input$BOR, c("CR" = 1, "PR"= 1,"SD"= 0,"PD"= 0))
  #################################


  pp<-sig_box(data = input, signature = "TMEscore", variable = "BOR3", jitter = jitter, palette = palette_box,hjust = 0.5)
  ggsave(pp, filename = paste0(index, "-1-1-TMEscore-BOR3-boxplot.",fig.type), path = file_name$folder_name,
         width = 5, height = 7)

  pp<-sig_box(data = input, signature = "TMEscore", variable = "BOR2", jitter = jitter, palette = palette_box,hjust = 0.5)
  ggsave(pp, filename = paste0(index, "-1-2-TMEscore-BOR2-boxplot.",fig.type), path = file_name$folder_name,
         width = 5, height = 7)


  pp<-sig_box(data = input, signature = "TMEscoreA", variable = "BOR3", jitter = jitter, palette = palette_box,hjust = 0.5)
  ggsave(pp, filename = paste0(index, "-1-3-TMEscoreA-BOR3-boxplot.",fig.type), path = file_name$folder_name,
         width = 5, height = 7)

  pp<-sig_box(data = input, signature = "TMEscoreA", variable = "BOR2", jitter = jitter, palette = palette_box,hjust = 0.5)
  ggsave(pp, filename = paste0(index, "-1-4-TMEscoreA-BOR2-boxplot.",fig.type), path = file_name$folder_name,
         width = 5, height = 7)

  pp<-sig_box(data = input, signature = "TMEscoreB", variable = "BOR3", jitter = jitter, palette = palette_box,hjust = 0.5)
  ggsave(pp, filename = paste0(index, "-1-5-TMEscoreB-BOR3-boxplot.",fig.type), path = file_name$folder_name,
         width = 5, height = 7)

  pp<-sig_box(data = input, signature = "TMEscoreB", variable = "BOR2", jitter = jitter, palette = palette_box,hjust = 0.5)
  ggsave(pp, filename = paste0(index, "-1-6-TMEscoreB-BOR2-boxplot.",fig.type), path = file_name$folder_name,
         width = 5, height = 7)


  pp<-sig_box(data = input, signature = "TMEscore", variable = "tissue", jitter = jitter, palette = palette_box,hjust = 0.5)

  if(addtable){
    ps<-as.matrix(addmargins(table(input$BOR3, input$tissue)))
    addTab <- as.data.frame(as.matrix(ifelse(round(ps, 3) < 0.001, "0",
                                             round(ps, 3))))
    df <- tibble(x = 0, y = 0, tb = list(addTab))
    pp <- pp +ggpp::geom_table_npc(data = df, aes(npcx = x, npcy = y, label = tb), table.rownames = TRUE)
    #############################################
  }

  ggsave(pp, filename = paste0(index, "-2-1-TMEscore-tissue-boxplot.",fig.type), path = file_name$folder_name,
         width = 5, height = 7)



  pp<-sig_box(data = input, signature = "TMEscore", variable = "Tissue_type", jitter = jitter, palette = palette_box,hjust = 0.5)

  if(addtable){
    ps<-as.matrix(addmargins(table(input$BOR3, input$Tissue_type)))
    addTab <- as.data.frame(as.matrix(ifelse(round(ps, 3) < 0.001, "0",
                                             round(ps, 3))))
    df <- tibble(x = 0, y = 0, tb = list(addTab))
    pp <- pp +ggpp::geom_table_npc(data = df, aes(npcx = x, npcy = y, label = tb), table.rownames = TRUE)
    #############################################
  }

  ggsave(pp, filename = paste0(index, "-2-2-TMEscore-Tissue_type-boxplot.",fig.type), path = file_name$folder_name,
         width = 5, height = 7)


  pp<-sig_box(data = input, signature = "TMEscore", variable = "batch", jitter = jitter, palette = palette_box,hjust = 0.5)

  if(addtable){
    ps<-as.matrix(addmargins(table(input$BOR3, input$batch)))
    addTab <- as.data.frame(as.matrix(ifelse(round(ps, 3) < 0.001, "0",
                                             round(ps, 3))))
    df <- tibble(x = 0, y = 0, tb = list(addTab))
    pp <- pp +ggpp::geom_table_npc(data = df, aes(npcx = x, npcy = y, label = tb), table.rownames = TRUE)
    #############################################
  }
  ggsave(pp, filename = paste0(index, "-2-3-TMEscore-batch-boxplot.",fig.type), path = file_name$folder_name,
         width = 5, height = 7)


  if(nlevels(as.factor(input$cohort))>2){

    if(addtable){
      pp<-sig_box(data = input, signature = "TMEscore", variable = "cohort", jitter = jitter, palette = palette_box,hjust = 0.5)
      ps<-as.matrix(addmargins(table(input$BOR3, input$cohort)))
      addTab <- as.data.frame(as.matrix(ifelse(round(ps, 3) < 0.001, "0",
                                               round(ps, 3))))
      df <- tibble(x = 0, y = 0, tb = list(addTab))
      pp <- pp +ggpp::geom_table_npc(data = df, aes(npcx = x, npcy = y, label = tb), table.rownames = TRUE)
      #############################################
    }
    ggsave(pp, filename = paste0(index, "-2-4-TMEscore-cohort-boxplot.",fig.type), path = file_name$folder_name,
           width = 5, height = 7)

  }


  sig_roc(data = input,
          response = "BOR01",
          variables =c("TMEscore","TMEscoreA","TMEscoreB"),
          fig.path = file_name$abspath,
          palette = palette_roc,
          file.name = paste0(index, "-3-TMEscore-ROC-plot"))


  #' @画AUC曲线图---------------------------

  multi_roc(data = input,
            response_name = "BOR01",
            var_names = c("TMEscore","TMEscoreA","TMEscoreB"),
            index = paste0(index, "-4"),
            ProjectID = "qPCR",
            target_signature = "TMEscore",
            palette = palette_roc,
            show_col = FALSE,
            compare = T,
            path = file_name$folder_name)



  #####################################
  target<-sym("TMEscore")
  cols<-IOBR::palettes(category = "box", palette = palette_his,
                       show_message = FALSE, show_col = FALSE, alpha = 0.75)
  # input$bor<-revalue(input$BOR, c("CR"= "CRPR", "SD" = "SD", "PD" = "PD", "PR"= "CRPR"))

  p<-ggplot(input, aes(x= !!target,fill= BOR3)) +
    geom_histogram(aes(y=..density..), binwidth = binwidth, colour = "black")+
    scale_fill_manual(values= cols)+
    geom_density(alpha=.2, fill="grey", weight = 0.7)+
    labs(title=  paste0("Distribution of TMEscore"),
         # subtitle= paste0(var, " of patient " ),
         caption = paste0(" Data of qPCR: ", ";  ","BC: best cutoff;   ",
                          "mono: monotherapy;   ", "com: combination;   ", date()))+

    xlab(paste0(target))+
    theme_light()+
    IOBR:: design_mytheme(legend.position = "bottom",axis_angle = 0)

  ggsave(p, filename = paste0(index,"-5-Distribution-of-TMEscore.",fig.type),
         width = 6.5, height = 5.5, path = file_name$folder_name)
  ################################################



  #################################
  #' @加载缩减后的基因
  (load("H:/01-1-GC-Project/3-4-GC-TMEscore-PanCan-TCGA/12-TaqMan-QPCR/13-Gene_reduction/0-Paten-genes.RData"))
  paten_gene
  genes_heat<-as.character(unlist(paten_gene[2:3]))

  message(">>> Signature genes after reduction will be present, except reference genes and ICB genes......")
  print(paten_gene)
  #################################

  input<-merge(input[!is.na(input$BOR),], eset, by="ID")
  # write_xlsx(input, paste0(file_name$abspath,"0-TMEscore-and-genes-method1.xlsx"))
  #################################
  summary(genes_heat%in%colnames(input))
  #################################
  idd<-c(1:ncol(input))[colnames(input)%in% c("TMEscoreA","TMEscoreB","TMEscore", genes_heat)]
  idd
  #################################
  head(input)
  for(i in idd){
    input[is.na(input[,i]),i]<-mean(input[!is.na(input[,i]),i])
  }
  message(paste0("Counts of NA after replacement = "), sum(is.na(eset)))
  ################################

  ################################
  input[,c("TMEscoreA","TMEscoreB","TMEscore", genes_heat)] <-scale(input[,c("TMEscoreA","TMEscoreB","TMEscore", genes_heat)])
  # help(sig_heatmap)
  ###############################


  #' @标准化后比较与其他基因的优劣-----------------
  #############################
  sig_roc(data = input,
          response = "BOR01",
          variables = c("TMEscore","CD274","PDCD1","CTLA4"),
          palette =  palette_roc,
          fig.path = file_name$abspath,
          file.name = paste0(index, "-6-TMEscore-and-ICBs-ROC-plot") )


  multi_roc(data = input,
            response_name = "BOR01",
            var_names = c("TMEscore","CD274","PDCD1","CTLA4"),
            index = paste0(index, "-7"),
            ProjectID = "qPCR",
            target_signature = "TMEscore",
            palette = palette_roc,
            show_col = FALSE,
            compare = T,
            path = file_name$folder_name)
  ###################################



  #################################
  input$lable<-paste0(input$ID,"_",input$name,"_",input$Tissue_type,"_",input$tissue,"_",input$batch)

  p<-sig_heatmap(input = input,
                 ID = "lable",
                 features = c("TMEscoreA","TMEscoreB","TMEscore", genes_heat),
                 group = "BOR3",
                 index = paste0(index,"-8"),
                 palette = 2,
                 show_heatmap_col_name = T,
                 path = file_name$folder_name,
                 width = 14,
                 height = 9,
                 angle_col = 60)

  tiff(filename = paste0(file_name$abspath,index,"-8-BOR3-tidyheatmap.tiff"), units = "in",
      width = 14, height = 9, res = 300)
  print(p)
  dev.off()
  #################################



  input$label2<-paste0(input$ID,"_",input$name,"_",input$Tissue_type,"_",input$tissue,"_",input$batch,"_",input$BOR)
  input$label2<-gsub(input$label2, pattern = "batch", replacement = "B")

  input$BOR4<-ifelse(input$BOR2=="CRPR","1","2")

  library(ggrepel)
  c<-ggplot(input,aes(x=BOR2,y= TMEscore,fill=BOR2,label=label2))+
    geom_boxplot(notch = F,alpha=0.95,outlier.colour = "black",
                 outlier.shape = NA,outlier.size = 5)+
    geom_jitter(width = 0,size=5)+
    scale_fill_manual(values= c("#EFC000FF","#0073C2FF" ,"#2E9FDF")) +
    ggtitle("TME signature score","stratified by BOR2")+
    ylab(paste0("TMEscore"))+
    geom_text_repel(
      data          = subset(input, input$BOR4== "2"),
      nudge_x       = 1.5,
      segment.size  = 0.1,
      segment.color = "grey50",
      direction     = "y",
      hjust         = 0,
      size          =3
    )+
    geom_text_repel(
      data          = subset(input, input$BOR4== "1"),
      nudge_x       = -1.7,
      segment.size  = 0.1,
      segment.color = "grey50",
      direction     = "y",
      hjust         = 0,
      size          =3
    )
  p<-c+design_mytheme()
  ggsave(p, filename = paste0(index,"-9-Distribution-of-TMEscore.",fig.type),
         width = 8.5, height = 12, path = file_name$folder_name)
  ################################################
  #################################


  #' @划分feature-genes-并加上标签
  score<-score[order(score$TMEscore,decreasing = F),]
  score$order<-1:nrow(score)
  score$cluster<- score[,BOR]
  summary(as.factor(score$cluster))
  ####################################
  score$label<-paste0(score$ID,"_",score$name,"_",score$batch)
  # score$label[score$cluster=="Notspecific"]<-NA
  ####################################
  head(score)
  ####################################
  score$BOR3<-plyr:: revalue(score$BOR, c("CR" = "CRPR", "PR"= "CRPR","SD"= "SD","PD"="PD"))
  ####################################
  pp<-ggdotchart(score, x = "label", y = "TMEscore",
                 color = "BOR3",                                # Color by groups
                 palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
                 sorting = "descending",                       # Sort value in descending order
                 add = "segments",                             # Add segments from y = 0 to dots
                 add.params = list(color = "lightgray", size = 2), # Change segment color and size
                 group = "BOR3",                                # Order by groups
                 dot.size = 6,                                 # Large dot size
                 label = round(score$TMEscore,1),                        # Add mpg values as dot labels
                 font.label = list(color = "white", size = 9,
                                   vjust = 0.5),               # Adjust label parameters
                 ggtheme = theme_pubr()                        # ggplot2 theme
  ) + geom_hline(yintercept = 0, linetype = 2, color = "lightgray")
  pp<-pp+design_mytheme(axis_angle = 60, hjust = 1)

  ggsave(pp,file = paste0(index,"-10-TMEscore-distribution.png"),
         width = 17,height = 8, path = file_name$folder_name)
  ####################################

  return(input)


}










