





#' Title
#'
#' @param score 
#' @param pcr_data 
#' @param patient_ids 
#' @param batch 
#' @param width 
#' @param index 
#' @param fig.path 
#' @param BOR 
#'
#' @return
#' @export
#'
#' @examples
tmescore_heatmap<-function(score, pcr_data, BOR = "BOR", patient_ids = NULL, batch = NULL, width,  index = 1, fig.path = "1-TMEscore-and-genes-heatmap"){
  
  
  
  #' @转置数据
  ################################
  eset<-pcr_data
  eset<-rownames_to_column(as.data.frame(t(eset)),var = "ID")
  summary(eset)
  ##################################
  
  file_name<- mydb::creat_folder(fig.path)
  
  score<-as.data.frame(score)
  
  colnames(score)[which(colnames(score)==BOR)]<-"BOR"
  
  message(">>> Data frame `score` should contaied variables: ID,TMEscoreA,TMEscoreB,TMEscore,BOR,name,Tissue_type,tissue, batch......")
  
  input<-score[,c("ID","TMEscoreA","TMEscoreB","TMEscore","BOR","name","Tissue_type","tissue","batch")]
  
  
  
  input<-score[,c("ID","TMEscoreA","TMEscoreB","TMEscore","BOR","name","Tissue_type","tissue","batch")]
  input$BOR3<-plyr:: revalue(input$BOR, c("CR" = "CRPR", "PR"= "CRPR","SD"= "SD","PD"="PD"))
  input$BOR2<-plyr:: revalue(input$BOR, c("CR" = "CRPR", "PR"= "CRPR","SD"= "SDPD","PD"="SDPD"))
  input$BOR01<-plyr:: revalue(input$BOR, c("CR" = 1, "PR"= 1,"SD"= 0,"PD"= 0))
  #################################
  #################################
  #' @加载缩减后的基因
  (load("H:/01-1-GC-Project/3-4-GC-TMEscore-PanCan-TCGA/12-TaqMan-QPCR/13-Gene_reduction/0-Paten-genes.RData"))
  paten_gene
  genes_heat<-as.character(unlist(paten_gene[2:3]))
  
  message(">>> Signature genes after reduction will be present, except reference genes and ICB genes......")
  # print(paten_gene)
  #################################
  
  input<-merge(input[!is.na(input$BOR),], eset, by="ID")
  # write_xlsx(input, paste0(file_name$abspath,"0-TMEscore-and-genes-method1.xlsx"))
  #################################
  summary(genes_heat%in%colnames(input))
  #################################
  idd<-c(1:ncol(input))[colnames(input)%in% c("TMEscoreA","TMEscoreB","TMEscore", genes_heat)]
  idd
  #################################
  
  message(">>> Genes with NA was replaced by mean value of none.na values.....")
  head(input)
  for(i in idd){
    input[is.na(input[,i]),i]<-mean(input[!is.na(input[,i]),i])
  }
  # message(paste0("Counts of NA after replacement = "), sum(is.na(eset)))
  ################################
  
  message(paste0(">>> Each Gene of ", dim(input)[1], " patients were scaled and normalizated....."))
  
  ################################
  input[,c("TMEscoreA","TMEscoreB","TMEscore", genes_heat)] <-scale(input[,c("TMEscoreA","TMEscoreB","TMEscore", genes_heat)])
  # help(sig_heatmap)
  ###############################
  
  #################################
  input$lable<-paste0(input$ID,"_",input$name,"_",input$Tissue_type,"_",input$tissue,"_",input$batch)
  
  
  if(!is.null(patient_ids))  input2<-input[input$ID%in%patient_ids,]
  
  if(!is.null(batch)) input2<-input[input$batch==batch,]

  if(is.null(patient_ids)&is.null(batch)){
    
    message(">>>> All patient will be display, user can adjust `patient_ids` or `batch` to choose patients")
    input2<- input
  }  
  
  if(nlevels(as.factor(input2$BOR3))<2) stop("Selected patients with no BOR difference.... user should add more patients with different BOR...")
  ##################################
  
  p<-sig_heatmap(input = input2,
                 ID = "lable",
                 features = c("TMEscoreA","TMEscoreB","TMEscore", genes_heat),
                 group = "BOR3", 
                 index = index,
                 palette = 2,
                 show_heatmap_col_name = T,
                 path = file_name$folder_name,
                 width = width, 
                 height = 9, 
                 angle_col = 60)
  
  tiff(filename = paste0(file_name$abspath,index,"-BOR3-tidyheatmap.tiff"), units = "in",
       width = width, height = 9, res = 300)
  print(p)
  dev.off()
  #################################
  
  return(input2)
  
}
