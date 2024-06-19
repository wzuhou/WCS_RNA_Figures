######################
#Venn plot
library(eulerr)

TISSUE=c("ADR","PIT","FAT",'GON','HEA' ,'HYP','LIV', 'PEC') 
full_TISSUE=c("Adrenal","Pituitary","Fat",'Gonad','Heart' ,'Hypothalamus' ,'Liver', 'Pectoralis Muscle')
full_name<-data.frame(
  TS=TISSUE,
  full_TS=full_TISSUE
)


EFFECT='GvsN'  # c('GvsN','Sex','LHS')
TS='PIT'

##For loop, iterate among tissues
for (TS in TISSUE) {
  compare=full_name[full_name[,1] == TS,2]
{  setwd(paste0('./WCS/WCS_RNA/',TS,'/effectBreed/Interaction/'))
  set1 <- read.table(paste0('./GENE_interact_',TS,'_EGG_LAY-GWCSvsNWCS_Volcano_ALL.txt'),header = F)$V1 
  set2 <- read.table(paste0('./GENE_interact_',TS,'_PBM-GWCSvsNWCS_Volcano_ALL.txt'),header = F)$V1
  set3 <- read.table(paste0('./GENE_interact_',TS,'_WINTER-GWCSvsNWCS_Volcano_ALL.txt'),header = F)$V1 
}

x4 <- list(
  Breeding = set1, 
  PreBasic_Molt = set2,
  Winter= set3
)

###nVennR-proportional-atypical
### I did not use this one
library(nVennR)
obj<-plotVenn(x4, setColors=c("#FF0000", "#00A08A","#F2AD00"), borderWidth = 0) #, "mediumorchid",'gray50'
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

#b <- paste0(capture.output(cat(a)))

fit1<-eval(parse(text =  paste0('euler(c(', a ,'))') 
))

p1 <- plot(fit1,
         #quantities = TRUE,
         quantities=list(cex=1.5),
         # fill = "transparent",
         fills = list(fill = c("#FC4E07",'darkblue',"#00AFBB","#E7B800"), alpha = 0.5),
         lty = 1, #1:4,
         legend = T,
         #labels = list(font =2,cex=0.8), #labels = list(font =2,cex=0.8,pos = 4,y =c(3,-3 ),just = c( "bottom")),
         main = list(label=compare,y=-0.2,cex=1.2))
#p1
ggplot2::ggsave(file=paste0("Venn_",TS,'_',EFFECT,'.png'), p1,
    width=5, height=4,dpi= 350) #, res=600
}
p1
dev.off()
