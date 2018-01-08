library(ggplot2)

setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/")
map<- read.csv("vesca2consensus_map_integerposns_2017-07-09.csv",header=FALSE)
map$V4<-map$V3/1000000
map$V1<- gsub('-', '.', map$V1)

setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/exf_Hap/")
myFiles <- list.files(pattern="MS35*")
y=data.frame()
for (i in 1:length(myFiles)){
  x<-read.csv(myFiles[i])
  lg<-gsub("_genotypes_cult.csv","",myFiles[i])
  lg<-gsub("MS35","",lg)
  x$LG<-rep(lg,nrow(x))
  y<-rbind(y,x)
}

ExFmark<-y
ExFmark$V4<-ExFmark$cM/1000000

setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/rxh_Hap/")
myFiles <- list.files(pattern="MS35*")
y=data.frame()
for (i in 1:length(myFiles)){
  x<-read.csv(myFiles[i])
  lg<-gsub("_genotypes_cult.csv","",myFiles[i])
  lg<-gsub("MS35","",lg)
  x$LG<-rep(lg,nrow(x))
  y<-rbind(y,x)
}

RxHmark<-y
RxHmark$V1<-gsub("-",".",RxHmark$V1)
RxHmark$LG<-gsub("3B.2","3B",RxHmark$LG)
RxHmark$LG<-gsub("4D.2","4D",RxHmark$LG)
RxHmark2<-merge(RxHmark,map,by.x="V1",by.y="V1", all.x= T)

setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/")
Srh<-read.csv("Significant_QTL_hk_BFRHcombmil.csv")
Srh$type<-rep("i90k",nrow(Srh))
SrhE<-read.csv("PSelectMSRHcombmilMS_Eqiv.csv")
SrhE$type<-rep("i35k",nrow(SrhE))
SR<-rbind(Srh,SrhE)
SR$Rname<- gsub('-', '.', SR$Rname)
Srh<-SR[c(4:9,18,19)]
Sef<-read.csv("Significant_QTL_hk_BFEFcombmil.csv")
Sef$type<-rep("i90k",nrow(Sef))
SefE<-read.csv("PSelectMSEFcombmilMS_Eqiv.csv")
SefE$type<-rep("i35k",nrow(SefE))
SE<-rbind(Sef,SefE)
SE$Rname<- gsub('-', '.', SE$Rname)
Sef<-SE[c(4:6,8:11,20,21)]
colnames(Sef)[8]<-"year"
Sef$pop<-rep("ExF",nrow(Sef))
Srh<- merge(Srh,map,by.x="Rname",by.y="V1", all.x=TRUE)
Srh<-Srh[c(2,1,10,3:8)]
colnames(Srh)[3]<-"pos"
colnames(Srh)[8]<-"year"
Srh$pop<-rep("RxH",nrow(Srh))
Sall<-rbind(Sef,Srh)
Sall$V4<-Sall$pos/1000000

# mapS<- map[which(map$V1 %in% Sall$Rname),]
# mapSC<-merge(mapS,Sall,by.x = "V1", by.y = "Rname",sort = FALSE)
# mapSC<-mapSC[-10,]
position_dodge(width = 1)
plot<- ggplot()+
  geom_point(data=map, aes(x=V2, y=V4), color='grey',pch="_", size = 3)+
  geom_point(data=ExFmark, aes(x=LG, y=V4), color='red3',pch="_", size = 2)+
  geom_point(data=RxHmark2, aes(x=V2.y, y=V4.y), color='slateblue4',pch="_", size = 2)+
  geom_line(data=map,aes(x=V2, y=V4),size=0.1, colour ="darkgrey")+
  #geom_point(data=mapSC, aes(x=V2, y=V3, colour=exp), size=4) +
  geom_point(data=Sall, aes(x=lg, y=V4, colour=pop, shape=type,lwd= 2), size=5) + 
  #geom_point(data=Correct_pos, aes(x=V2, y=Mb), size=4) + 
  scale_shape_manual(values=c(1,2))+
  scale_color_manual(values=c('red','slateblue4'))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),legend.key = element_rect(fill = "white",size = 0.01), panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_y_reverse()+
  labs(x="Chromosome",y="Location (Mb)")
#theme(legend.position='none')

setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/figs/")
pdf("Sig Fig 2.pdf",width=12,height=8,paper='special') 
plot
dev.off()
