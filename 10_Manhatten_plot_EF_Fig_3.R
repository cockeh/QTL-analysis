install.packages("glmnet")
install.packages("ggplot2")
install.packages("plyr")

library(ggplot2)
library(scales)
library(gridExtra)
library(gtable)
library(grid)
library(dplyr)
library(data.table)
library(plyr)
library(glmnet)
library(ggplot2)

setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/")
X13<-read.csv("Sig_valuesmil13.csv")
X13$Year <- rep("13a",nrow(X13)) 
X13D<-read.csv("Sig_valuesmil13D.csv")
X13D$Year <- rep("13b",nrow(X13D)) 
X14<-read.csv("Sig_valuesmil14.csv")
X14$Year <- rep(14,nrow(X14)) 
X11<-read.csv("Sig_valuesmil11.csv")
X11$Year <- rep(11,nrow(X11)) 
X12<-read.csv("Sig_valuesmil12.csv")
X12$Year <- rep("12a",nrow(X12)) 
X12D<-read.csv("Sig_valuesmil12D.csv")
X12D$Year <- rep("12b",nrow(X12D)) 
BLUP<-read.csv("Sig_valuesEFcombmil.csv")
BLUP$Year <- rep("Combined",nrow(BLUP)) 
allmil<-rbind(X11,X12,X12D,X13,X13D,X14,BLUP)
allmil$extra <-allmil$pos
# allmil <- allmil[order(allmil$pos),]
# allmil <- allmil[order(allmil$lg),]

PosVal<- data.frame(Pv = numeric())
for (i in unique(allmil$lg)) {
  x <- tail(allmil[allmil$lg == i,], n=1)
  z<- x$pos 
  PosVal[i,] = z   
}
setDT(PosVal, keep.rownames = TRUE)[]
newrow<- NA
PosVal <- PosVal[order(PosVal$rn),]
PosVal = rbind(newrow,PosVal,fill=TRUE)
nc  <- nrow(PosVal)
PosVal[1, 3] <- 0

shift <- function(x, n){
  c(x[-(seq(n))])
}
PosVal<- rbind(PosVal, "29"=NA, fill=TRUE)
PosVal$rn <- shift(PosVal$rn, 1)
#error message is not an issue
PosVal$x <- NULL
PosVal <- PosVal[-30, ]
PosVal$lg <- PosVal$rn
PosVal$Pv2<-cumsum(PosVal$Pv)

PosVal$lty<- rep(3,nrow(PosVal))
thick<- c(1,5,9,13,17,21,25,29)
PosVal$lty[thick]<-gsub("3","1",PosVal$lty[thick])

plyr1 <- join(allmil, PosVal, by = 'lg')
plyr1$position <- plyr1$extra + plyr1$Pv2
EMallmil<- plyr1[plyr1$mtype == "lmxll",]
threshold<-0.05
threshold2<-0.01

for (i in 1:nrow(PosVal)){
  i2<-i+1
  textpos<- PosVal[i,4]+((PosVal[i2,4]-PosVal[i,4])/2)
  print(textpos)
  PosVal$V5[i]<-textpos
}
print(unlist(PosVal$V5))

EMallmil$xforline<-rep(6,nrow(EMallmil))

EMp <- ggplot(EMallmil , aes(x=position, y=-log10(pval), colour=Year, group=Year, size=Year)) +
  #geom_point(pch=1) +
  geom_line()+
  ggtitle("Emily")+ 
  scale_color_manual(values=c('#7CAE00','#00BFC4',"#00BA38","#F8766D","#619CFF","#FF61CC","black","black"))+
  geom_vline(data=PosVal,aes(xintercept=Pv2,linetype=lty, colour ="black"))+
  #scale_color_manual(values=c('grey8','grey20','grey32','grey44','grey56','grey68',"black"))+
  #scale_y_log10(limits = c(0.0000001,1))+
  scale_size_manual(values=c(0.3,0.3,0.3,0.3,0.3,0.3,0.4))+
  geom_hline(aes(yintercept=1.3))+
  geom_hline(aes(yintercept=2), lty=2)+
  annotate("text", x = unlist(PosVal$V5[1:28]), y = -0.2, label = c("1A","1B","1C","1D","2A","2B","2C","2D","3A","3B","3C","3D","4A","4B","4C","4D","5A","5B","5C","5D","6A","6B","6C","6D","7A","7B","7C","7D"))+
  scale_x_continuous(expand = c(0, 0))+
  coord_cartesian(xlim=c(0,B29),ylim=c(-0.5, 6))+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'), axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x="Position (Mb)",y=expression(-log[10](italic(p))))+
  theme(legend.position='none')
EMp


FEallmil<-plyr1[plyr1$mtype == "nnxnp",]

FEp <- ggplot(FEallmil , aes(x=position, y=-log10(pval), colour=Year, group=Year, size=Year)) +
  #geom_point(pch=1) +
  geom_line()+
  ggtitle("Fenella")+ 
  scale_color_manual(values=c('#7CAE00','#00BFC4',"#00BA38","#F8766D","#619CFF","#FF61CC","black","black"))+
  geom_vline(data=PosVal,aes(xintercept=Pv2,linetype=lty, colour ="black"))+
  #scale_color_manual(values=c('grey8','grey20','grey32','grey44','grey56','grey68',"black"))+
  #scale_y_log10(limits = c(0.0000001,1))+
  scale_size_manual(values=c(0.3,0.3,0.3,0.3,0.3,0.3,0.4))+
  geom_hline(aes(yintercept=1.3))+
  geom_hline(aes(yintercept=2), lty=2)+
  geom_vline(aes(xintercept= 0))+
  annotate("text", x = unlist(PosVal$V5[1:28]), y = -0.2, label = c("1A","1B","1C","1D","2A","2B","2C","2D","3A","3B","3C","3D","4A","4B","4C","4D","5A","5B","5C","5D","6A","6B","6C","6D","7A","7B","7C","7D"))+
scale_x_continuous(expand = c(0, 0))+
  coord_cartesian(xlim=c(0,B29),ylim=c(-0.5, 6))+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'), axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x="Position (Mb)",y=expression(-log[10](italic(p))))+
  theme(legend.position='none')
FEp

Ballmil<-plyr1[plyr1$mtype == "hkxhk",]

Bp <- ggplot(Ballmil , aes(x=position, y=-log10(pval), colour=Year, group=Year, size=Year)) +
  #geom_point(pch=1) +
  geom_line()+
  ggtitle("Shared")+ 
  scale_color_manual(values=c('#7CAE00','#00BFC4',"#00BA38","#F8766D","#619CFF","#FF61CC","black","black"))+
  #scale_color_manual(values=c('grey8','grey20','grey32','grey44','grey56','grey68',"black"))+
  #scale_y_log10(limits = c(0.0000001,1))+
  scale_size_manual(values=c(0.3,0.3,0.3,0.3,0.3,0.3,0.4))+
  geom_hline(aes(yintercept=1.3))+
  geom_hline(aes(yintercept=2), lty=2)+
  geom_vline(aes(xintercept= 0))+
  geom_vline(data=PosVal,aes(xintercept=Pv2,linetype=lty, colour ="black"))+
  annotate("text", x = unlist(PosVal$V5[1:28]), y = -0.2, label = c("1A","1B","1C","1D","2A","2B","2C","2D","3A","3B","3C","3D","4A","4B","4C","4D","5A","5B","5C","5D","6A","6B","6C","6D","7A","7B","7C","7D"))+
scale_x_continuous(expand = c(0, 0))+
  coord_cartesian(xlim=c(0,B29),ylim=c(-0.5, 6))+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'), axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x="Position (Mb)",y=expression(-log[10](italic(p))))+
  theme(legend.position='none')
Bp

setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/figs")
pdf("Fig3.pdf",width=12,height=9.6,paper='special') 
# run multiplot from here http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot(EMp,FEp,Bp, cols=1)
dev.off()

