
library(ggplot2)
library(graphics)
library(gplots)

setwd("/Users/helencockerton/Desktop/QTL/Files for git hub")

score<-read.csv("./data/RH_EF_relative_AUDPC2.csv")

Mean<-with(na.omit(score),tapply(audpc_r,year,mean)) 
Std<-with(na.omit(score),tapply(audpc_r,year,sd))
SE<- Std/(sqrt(10))
CI<- 1.96*SE
par(oma=c(3,3,3,3))
par(las=2)
par(mar=c(7, 6, 4, 2)+ 0.1)
barplot2(Mean,plot.ci=TRUE,ci.l=Mean-SE,ci.u=Mean+SE, axis.lty=1,ylim=c(0,450), border = "black", col= c("#373737","#696969","#696969","#9B9B9B","#9B9B9B","#CDCDCD","#696969","#9B9B9B","#CDCDCD","white"), plot.grid= TRUE, ylab="Relative AUDPC", names.arg=c("2011","2012a", "2012b","2013a","2013b","2014","2012","2013","2014","2016")) 
box()

