## Rscript
{library(ggplot2)
  library(tidyverse)
  library(ggrepel)
  library(wesanderson)}

TS='PEC'
EFFECT='GvsN' # c('GvsN','Sex','LHS')

{
setwd(paste0('./WCS/WCS_RNA/',TS,'/effectBreed/Interaction/'))
#1 input 
dfBR <- read.table(paste0("DESeq2_out_interact_",TS,"_EGG_LAY-GWCSvsNWCS.txt"),header = T)[,c('padj','log2FoldChange','Gene')]
head(dfBR)
dfBR$cluster = 'BR'

dfBR$label <- ifelse(dfBR$padj<0.01 & abs(dfBR$log2FoldChange) >= 0.5849625,"Adj P<0.01","NS")
head(dfBR)

#2 input 
dfPBM <- read.table(paste0("DESeq2_out_interact_",TS,"_PBM-GWCSvsNWCS.txt"),header = T)[,c('padj','log2FoldChange','Gene')]
head(dfPBM)
dfPBM$cluster = 'PBM'

dfPBM$label <- ifelse(dfPBM$padj < 0.01 & abs(dfPBM$log2FoldChange) >= 0.5849625,"Adj P<0.01","NS")
head(dfPBM)

#3 input 
dfWI <- read.table(paste0("DESeq2_out_interact_",TS,"_WINTER-GWCSvsNWCS.txt"),header = T)[,c('padj','log2FoldChange','Gene')]
head(dfWI)
dfWI$cluster = 'WI'

dfWI$label <- ifelse(dfWI$padj<0.01 & abs(dfWI$log2FoldChange)>= 0.5849625,"Adj P<0.01","NS")
head(dfWI)

#combine
df <- rbind(dfBR,dfPBM,dfWI)
top10sig1 <- filter(df,cluster=="BR") %>% distinct(Gene,.keep_all = T) %>% top_n(10,-log(abs(padj)))
head(top10sig1)

top10sig2 <- filter(df,cluster=="PBM") %>% distinct(Gene,.keep_all = T) %>% top_n(10,-log(abs(padj)))
head(top10sig2)

top10sig3 <- filter(df,cluster=="WI") %>% distinct(Gene,.keep_all = T) %>% top_n(10,-log(abs(padj)))
head(top10sig3)

#extreact cluster's Top10,combine
top10sig <- rbind(top10sig1,top10sig2,top10sig3)
#top10sig<-top10sig0
#New col,Top10==2，Otherwise1；
df$size <- case_when(!(df$Gene %in% top10sig$Gene)~ 1,
                     df$Gene %in% top10sig$Gene ~ 2)

#NON-Top10gene list
dt <- filter(df,size==1)
head(dt)

#Using log2FoldChange  to determine the lim ##max and min of y axis
dfbar <- data.frame(x=c("BR" , "PBM",'WI'),
                  y=c(max(df[df$'cluster'== 'BR' , c('log2FoldChange')]),
                      max(df[df$'cluster'== 'PBM' , c('log2FoldChange')]),
                      max(df[df$'cluster'== 'WI' , c('log2FoldChange')])))
dfbar1 <- data.frame(x=c("BR" , "PBM",'WI'),
                   y=c(min(df[df$'cluster'== 'BR' , c('log2FoldChange')]) ,
                       min(df[df$'cluster'== 'PBM' , c('log2FoldChange')]),
                       min(df[df$'cluster'== 'WI' , c('log2FoldChange')])))

#X-axis cluster color
dfcol<-data.frame(x=c("BR" , "PBM",'WI'),
                  y=0) #,label=c(0:8))
mycol <-c(wes_palette("Darjeeling1",n=5,type ="discrete"),wes_palette("Darjeeling2",n=3,type ="discrete"))
mycol <- c("#FF0000", "#00A08A","#F2AD00")


ggplot()+
  geom_col(data = dfbar,mapping = aes(x = x,y = y),fill = "gray90",alpha = 0.5)+
  geom_col(data = dfbar1,mapping = aes(x = x,y = y),fill = "gray90",alpha = 0.5)+
  geom_jitter(data = dt, position = position_jitter(seed = 1),
              aes(x = cluster, y = log2FoldChange, color = label),
              #width =0.4,
              size = 0.8)+
  geom_jitter(data = top10sig,position = position_jitter(seed = 1),
              aes(x = cluster, y = log2FoldChange, color = label),
              #width =0.4,
              size = 1)+
  geom_tile(data = dfcol,
            aes(x=x,y=y),
            height=0.4,color = "black",fill = mycol,alpha = 0.6,show.legend = F)+
  #Add label to Top10 genes
  geom_text_repel(
    data=top10sig,position = position_jitter(seed = 1),
    aes(x=cluster,y=log2FoldChange,label=Gene),
    point.padding = 0, # additional padding around each point
    box.padding = 0.3,
    #size=2.5,
    arrow =arrow(length = unit(0.008, "npc"),type = "open", ends = "last"),
    max.overlaps = 16)+
  #color of dots
  scale_color_manual(values = c("#C93312","#899DA4"))+
  #X/Y axis, Cluster label
  labs(x="",y="log2FC",title = paste0(TS,'_',EFFECT))+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=x),
            size =6,color ="white",fontface='bold')+
  #Theme
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "gray10",
                               size = 1),
    axis.line.x = element_blank(),
    axis.title.x =element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = 'none' , #"top",
    #legend.direction = "vertical",    legend.justification = c(1,0),    legend.text = element_text(size = 8),
    plot.title = element_text(size=13.5,face = "bold") #,hjust = 0.5
    ) -> p_HYP
#p_HYP
ggsave(paste0('Joint_volcano_interact_',TS,'-GWCSvsNWCS.png'),p_HYP,width = 9,height = 5.5,dpi=350)

}
