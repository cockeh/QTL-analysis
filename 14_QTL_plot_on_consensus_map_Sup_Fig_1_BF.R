# Lplot

library(ggplot2)

setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/")
map<- read.csv("vesca2consensus_map_integerposns_2017-07-09.csv",header=FALSE)
map$V4<-map$V3/1000000
map$V1<- gsub('-', '.', map$V1)

Sef<-read.csv("All_Significant_QTL_hk_EF_BF.csv")
Sef$Pop <- rep("ExF",nrow(Sef))
colnames(Sef)[21] <-"Year"
Sef$Year<- gsub('Sig_valuesmil', '20', Sef$Year)
Sef$Year<- gsub('.csv', '', Sef$Year)
Sef<-Sef[,c(5:7,9:12,21,22)]
Srh<-read.csv("All_Significant_QTLhkBF.csv")
Srh$Pop <- rep("RxH",nrow(Srh))
colnames(Srh)[19] <-"Year"
Srh$Year<- gsub('Sig_valuesRGH_mil_', '20', Srh$Year)
Srh$Year<- gsub('.csv', '', Srh$Year)
Srh<-Srh[,c(5:10,19,20)]
Srh<- merge(Srh,map,by.x="Rname",by.y="V1", all.x=TRUE)
Srh<-Srh[c(2,1,10,3:8)]
colnames(Srh)[3]<-"pos"
Sall <- rbind(Sef,Srh)
mapS<- map[which(map$V1 %in% Sall$Rname),]
mapSrh<- map[which(map$V1 %in% Srh$Rname),]
mapSef<- map[which(map$V1 %in% Sef$Rname),]

mapSC<-merge(mapS,Sall,by.x = "V1", by.y = "Rname",sort = FALSE)
mapSC$linewidth<- mapSC$sig
mapSC$linewidth<-as.character(mapSC$linewidth)
for(i in 1:length(mapSC$linewidth)){
  mapSC$linewidth[i]<-nchar(mapSC$linewidth[i])
} 
mapSC$linewidth<-as.numeric(mapSC$linewidth)
mapSC$linewidth2<-mapSC$linewidth /2

map$V4<-map$V3/1000000
mapSC$V4<-mapSC$V3/1000000

position_dodge(width = 1)
plot<- ggplot()+
  #ggtitle("Mildew QTL")+ 
  geom_point(data=map, aes(x=V2, y=V4), color='grey',pch="_", size = 3)+
  geom_line(data=map,aes(x=V2, y=V4),size=0.1, colour ="gray35")+
  geom_point(data=mapSC, aes(x=V2, y=V4, colour=Pop, group=Pop, shape=Year), size=4) + 
  #, stroke=linewidth2
  scale_shape_manual(values=c(1,2,5,6,7,9,10))+
  scale_color_manual(values=c('red3','slateblue'))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),legend.key = element_rect(fill = "white",size = 0.01), panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_y_reverse()+
  labs(x="Chromosome",y="Location (Mb)")
#guides(colour=FLASE)

setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/figs/")
pdf("Sup Fig 1_BF.pdf",width=12,height=8,paper='special') 
plot
dev.off()
