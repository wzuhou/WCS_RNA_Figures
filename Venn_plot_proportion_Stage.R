######################
#Venn plot
library(eulerr)
TISSUE=c("ADR","PIT","FAT",'GON','HEA' ,'HYP','LIV', 'PEC') #j
LHS=c('EGG_LAY','PBM','WINTER')
full_TISSUE=c("Adrenal","Pituitary","Fat",'Gonad','Heart' ,'Hypothalamus' ,'Liver', 'Pectoralis Muscle')
TS='PIT'
#EFFECT='GvsN' # c('GvsN','Sex','LHS')
{dat <- tidyr::crossing(x=1:3,y=1:3) #Case=y
  dat2 <- dat %>%
    filter(!row_number()==2)%>%
    group_by(grp = paste0(pmin(x, y), pmax(x, y))) %>%
    filter(!x==y) %>%
    filter(row_number() == max(row_number())) %>%
    ungroup() %>%
    # filter(!grp%in%c("13","24"))%>%
    select(-grp)}
dat2

full_name<-data.frame(
  TS=TISSUE,
  full_TS=full_TISSUE
)


for (TS in TISSUE) {
  
  setwd(paste0('C:/Users/zwu33/Downloads/ROSLIN/WCS/WCS_RNA/',TS,'/effectStage/Interaction/'))
  for (d in (1:3)){
    #d=2
    Case= LHS[dat2$y[d]] ; Case #EGG_LAY
    Base= LHS[dat2$x[d]] ; Base # WINTER
    EFFECT=paste(Case,Base,sep="vs")
    compare=paste0(full_name[full_name[,1] == TS,2],'\n',EFFECT)
    
    
    set1 <- read.table(paste0('./GENE_interact_',TS,'_GWCS-',EFFECT,'_Volcano_ALL.txt'),header = F)$V1 #GWCS
    set2 <- read.table(paste0('./GENE_interact_',TS,'_NWCS-',EFFECT,'_Volcano_ALL.txt'),header = F)$V1 #NWCS
    #set3 <- read.table(paste0('./GENE_interact_',TS,'_WINTER-GWCSvsNWCS_Volcano_ALL.txt'),header = F)$V1 #Fat up
    # set4 <- read.table(paste0("M:/ROSLIN/RNA_2LALO/FAT/GENE_FAT_Conditions_4_VS_3_Volcano_DOWN.txt"),header=F) #Fat down
    
    x4 <- list(
      GWCS = set1, 
      NWCS = set2
      #Winter= set3
    )
    
    ###nVennR-proportional-atypical
    ### I did not use this one
    library(nVennR)
    obj<-plotVenn(x4, setColors=c('#3B9AB2','#F21A00'), borderWidth = 0) #, "mediumorchid",'gray50'
    df2<-listVennRegions(obj)
    df3<-lapply(df2,length)
    
    df4<-stringr::str_split(names(df3),'[[(]]',simplify = T) [,2]
    df4<-stringr::str_replace(df4,'[[)]]' , "" )
    df4<-stringr::str_replace_all(df4,',' , "&" )
    df4<-stringr::str_replace_all(df4,' ' , "" )
    
    df3<-as.data.frame(df3)
    names(df3)<-df4
    df3
    
    ######venneuler Proportional 
    #Sys.setenv(JAVA_HOME="C:\\Program Files (x86)\\Java\\jre1.8.0_311") #change R to 32 bit #install.packages("rJava")  
    # library(rJava)
    # library("venneuler")
    # plot(venneuler(c("Snowstorm"=1539,"Extreme_spring"=70,"Extreme_spring&Snowstorm"=92
    #                  )))
    
    a <-paste0(capture.output(for (i in 1:length(df3)){
      GRP=colnames(df3)[i]
      VALUE=df3[i]
      #print(GRP)
      cat(paste0('"',GRP,'"',"=",VALUE,','))
    } ), collapse="\n") 
    
    a<-stringr::str_sub(a, end=-2)
    #or str_sub(iris$Species, 1, str_length(iris$Species)-3)
    #b <- paste0(capture.output(cat(a)))
    
    
    #####manually copy-paste#####
    # fit1<-euler(c(
    #   "Winter"=111,"PreBasic_Molt"=219,"PreBasic_Molt&Winter"=14,"Breeding"=35,"Breeding&Winter"=2,"Breeding&PreBasic_Molt"=10,"Breeding&PreBasic_Molt&Winter"=1
    # ))
    
    fit1<-eval(parse(text = paste0('euler(c(', a ,'))')
    ))
    
    p1 <- plot(fit1,
               #quantities = TRUE,
               quantities=list(cex=1.5),
               # fill = "transparent",
               fills = list(fill = c( '#3B9AB2','#F21A00'), alpha = 0.6),
               lty = 1, #1:4,
               legend = T,
               #labels = list(font =2,cex=0.8), #labels = list(font =2,cex=0.8,pos = 4,y =c(3,-3 ),just = c( "bottom")),
               main = list(label=compare,y=-0.2,cex=1.2))
    #p1
    #dev.off()
    
    ggplot2::ggsave(file=paste0("Venn_",TS,'_',EFFECT,'.png'), p1,
                    width=5, height=4,dpi= 350) #, res=600
  }}
p1
dev.off()


#############################
#####Focus on Gonad##########
#############################
#############################
# Focus on Sex###############
#############################
library(eulerr)
TISSUE=c("ADR","PIT","FAT",'GON','HEA' ,'HYP','LIV', 'PEC') #j
full_TISSUE=c("Adrenal","Pituitary","Fat",'Gonad','Heart' ,'Hypothalamus' ,'Liver', 'Pectoralis Muscle')
TS='GON'

full_name<-data.frame(
  TS=TISSUE,
  full_TS=full_TISSUE
)
for (TS in TISSUE) {
  setwd(paste0('C:/Users/zwu33/Downloads/ROSLIN/WCS/WCS_RNA/',TS,'/effectStage/Indep'))
  for (d in (1:3)){
   # d=2
    Case= LHS[dat2$y[d]] ; Case #EGG_LAY
    Base= LHS[dat2$x[d]] ; Base # WINTER
    EFFECT=paste(Case,Base,sep="vs")
    compare=paste0(full_name[full_name[,1] == TS,2],'\n',EFFECT) ;compare
    
    set1 <- read.table(paste0('./GENE_independent_',TS,'_GWCS_F','-',EFFECT,'_Volcano_ALL.txt'),header = F)$V1 #GWCS
    set2 <- read.table(paste0('./GENE_independent_',TS,'_GWCS_M','-',EFFECT,'_Volcano_ALL.txt'),header = F)$V1 #GWCS
    set3 <- read.table(paste0('./GENE_independent_',TS,'_NWCS_F','-',EFFECT,'_Volcano_ALL.txt'),header = F)$V1 #NWCS
    set4 <- read.table(paste0('./GENE_independent_',TS,'_NWCS_M','-',EFFECT,'_Volcano_ALL.txt'),header = F)$V1 #NWCS
    
    x4 <- list(
      GWCS_F = set1, 
      GWCS_M = set2,
      NWCS_F = set3, 
      NWCS_M = set4
      #Winter= set3
    )
    
    
    library(nVennR)
    obj<-plotVenn(x4, setColors=c("#ECCBAE", "#046C9A" ,"#D69C4E", "#ABDDDE"), borderWidth = 0) #, "mediumorchid",'gray50'
    df2<-listVennRegions(obj)
    df3<-lapply(df2,length)
    
    df4<-stringr::str_split(names(df3),'[[(]]',simplify = T) [,2]
    df4<-stringr::str_replace(df4,'[[)]]' , "" )
    df4<-stringr::str_replace_all(df4,',' , "&" )
    df4<-stringr::str_replace_all(df4,' ' , "" )
    
    df3<-as.data.frame(df3)
    names(df3)<-df4
    df3
######
    
    a <-paste0(capture.output(for (i in 1:length(df3)){
      GRP=colnames(df3)[i]
      VALUE=df3[i]
      #print(GRP)
      cat(paste0('"',GRP,'"',"=",VALUE,','))
    } ), collapse="\n") 
    
    a<-stringr::str_sub(a, end=-2)
    #or str_sub(iris$Species, 1, str_length(iris$Species)-3)
    #b <- paste0(capture.output(cat(a)))
    
    
    #####manually copy-paste#####
    # fit1<-euler(c(
    #   "Winter"=111,"PreBasic_Molt"=219,"PreBasic_Molt&Winter"=14,"Breeding"=35,"Breeding&Winter"=2,"Breeding&PreBasic_Molt"=10,"Breeding&PreBasic_Molt&Winter"=1
    # ))
    
    fit1<-eval(parse(text = paste0('euler(c(', a ,'))')
    ))
    
    p1 <- plot(fit1,
               #quantities = TRUE,
               quantities=list(cex=1.5),
               # fill = "transparent",
               fills = list(fill = c("#ECCBAE", "#046C9A" ,"#D69C4E", "#ABDDDE"), alpha = 0.6),
               lty = 1, #1:4,
               legend = T,
               #labels = list(font =2,cex=0.8), #labels = list(font =2,cex=0.8,pos = 4,y =c(3,-3 ),just = c( "bottom")),
               main = list(label=compare,y=-0.2,cex=1.2))
    #p1
    #dev.off()
    
    ggplot2::ggsave(file=paste0("Venn_",TS,'_',EFFECT,'.png'), p1,
                    width=8, height=6,dpi= 350) #, res=600
  }
  }
p1
dev.off()
