# Lplot

library(ggplot2)

setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/")
map<- read.csv("vesca2consensus_map_integerposns_2017-07-09.csv",header=FALSE)
map$V4<-map$V3/1000000
map$V1<- gsub('-', '.', map$V1)

Srh<-read.csv("Significant_QTL_hk_BFRHcombmil.csv")
Srh<-Srh[c(4:9)]
Sef<-read.csv("Significant_QTL_hk_BFEFcombmil.csv")
Sef<-Sef[c(4:6,8:11)]
Sef$pop<-rep("ExF",nrow(Sef))

Srh<- merge(Srh,map,by.x="Rname",by.y="V1", all.x=TRUE)
Srh<-Srh[c(2,1,8,3:6)]
colnames(Srh)[3]<-"pos"
Srh$pop<-rep("RxH",nrow(Srh))

Sgwas<- read.csv("No_Mil2.txt_istraw90.out_fix_min10_pheno_0.2_strat.assoc.linear.csv")
Sgwas$SNP<- gsub('-', '.', Sgwas$SNP)
Sgwas<- Sgwas[Sgwas$UNADJ == 0.0000120,]
position_of_markers<-read.csv("vesca2consensus_map_integerposns_2017-07-09.csv", header=FALSE)
GWAS_qtl<-merge(Sgwas,position_of_markers, by.x="SNP",by.y="V1", all.x=TRUE) 
GWAS_qtl<-GWAS_qtl[,c(11,1,12,3)]
GWAS_qtl$sig<-rep("****",nrow(GWAS_qtl))
GWAS_qtl$pop<-rep("GWAS",nrow(GWAS_qtl))
colnames(GWAS_qtl)<-c("lg","Rname","pos","pval","sig","pop") 
Sall<-rbind(Sef,Srh)
Sall<-Sall[c(1:4,6,8)]
Sall<-rbind(Sall,GWAS_qtl)

mapS<- map[which(map$V1 %in% Sall$Rname),]
mapSC<-merge(mapS,Sall,by.x = "V1", by.y = "Rname",sort = FALSE)
mapSC$linewidth<- mapSC$sig
mapSC$linewidth<-as.character(mapSC$linewidth)
for(i in 1:length(mapSC$linewidth)){
  mapSC$linewidth[i]<-nchar(mapSC$linewidth[i])
} 
mapSC$linewidth<-as.numeric(mapSC$linewidth)
mapSC$linewidth2<-mapSC$linewidth /2
mapSC$Mb<-mapSC$V3/1000000

position_dodge(width = 1)

plot<- ggplot()+
  geom_point(data=map, aes(x=V2, y=V4), color='grey',pch="_", size = 3)+
  geom_line(data=map,aes(x=V2, y=V4),size=0.1, colour ="gray35")+
  geom_point(data=mapSC, aes(x=V2, y=Mb, colour=pop, group=pop,stroke=linewidth) , pch=4, size=4) + 
  #geom_point(data=Correct_pos, aes(x=V2, y=Mb,stroke=linewidth), pch=4, size=4) + 
  #scale_shape_manual(values=c(1,2,5,6,7,9,10))+
  scale_color_manual(values=c('red3',"green4",'slateblue4'))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),legend.key = element_rect(fill = "white",size = 0.01), panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_y_reverse()+
  labs(x="Chromosome",y="Location (Mb)")+
  theme(legend.position='none')

setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/figs/")
pdf("Fig_5_BF.pdf",width=12,height=8,paper='special') 
plot
dev.off()










