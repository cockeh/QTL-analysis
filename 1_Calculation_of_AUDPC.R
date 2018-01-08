# Calculation of AUDPC
# Script written before the discorey of for loops.... and yet it works so why fix it

library(gplots)
library(lattice)
library(ggplot2)
library(agricolae)
library(plyr)
library(colorspace)
library(car)
library(lme4)


EF_11  <-	9
EF_12	<-	9
EF_12D	<-	9
EF_13  <-  14
EF_13D	<-	26
EF_14	<-	57
RH_12  <-  21
RH_13	<-	14
RH_14	<-	57
RH_16	<-	17


setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data")
pheno<-read.csv("exf2011.csv")
head(pheno)

days<-c(327,336)
pheno[,1]=as.factor(pheno[,1])
pheno[,5]=audpc(pheno[,3:4],days)
pheno[,6]= pheno[,5]/ EF_11 *100
colnames(pheno)[5]=c('audpc')
colnames(pheno)[6]=c('audpc_r')

write.csv(pheno, file = "emxf2011A.csv", row.names=FALSE)

pheno<-read.csv("emxf2012.csv")
head(pheno)

days<-c(327,336)
pheno[,1]=as.factor(pheno[,1])
pheno[,2]=as.factor(pheno[,2])

pheno[,5]=audpc(pheno[,3:4],days)
pheno[,6]=pheno[,5]/ EF_12 *100
colnames(pheno)[5]=c('audpc')
colnames(pheno)[6]=c('audpc_r')
write.csv(pheno, file="emxf2012A.csv", row.names=FALSE)


pheno<-read.csv("emxf2012D.csv")
head(pheno)

days<-c(327,336)
pheno[,1]=as.factor(pheno[,1])
pheno[,2]=as.factor(pheno[,2])

pheno[,5]=audpc(pheno[,3:4],days)
pheno[,6]=pheno[,5]/ EF_12D *100
colnames(pheno)[5]=c('audpc')
colnames(pheno)[6]=c('audpc_r')
write.csv(pheno, file="emxf2012DA.csv", row.names=FALSE)

pheno<-read.csv("emxf2013.csv")
head(pheno)

days<-c(325,339)
pheno[,1]=as.factor(pheno[,1])
pheno[,2]=as.factor(pheno[,2])

pheno[,5]=audpc(pheno[,3:4],days)
pheno[,6]=pheno[,5]/ EF_13 *100
colnames(pheno)[5]=c('audpc')
colnames(pheno)[6]=c('audpc_r')
write.csv(pheno, file="emxf2013A.csv", row.names=FALSE)

pheno<-read.csv("emxf2013D.csv")
head(pheno)

days<-c(218,242)
pheno[,1]=as.factor(pheno[,1])
pheno[,2]=as.factor(pheno[,2])

pheno[,5]=audpc(pheno[,3:4],days)
pheno[,6]=pheno[,5]/ EF_13D *100
colnames(pheno)[5]=c('audpc')
colnames(pheno)[6]=c('audpc_r')
write.csv(pheno, file="emxf2013DA.csv", row.names=FALSE)

pheno<-read.csv("emxf2014.csv")
head(pheno)

days<-c(281,338)
pheno[,1]=as.factor(pheno[,1])
pheno[,2]=as.factor(pheno[,2])

pheno[,5]=audpc(pheno[,3:4],days)
pheno[,6]=pheno[,5]/ EF_14 *100
colnames(pheno)[5]=c('audpc')
colnames(pheno)[6]=c('audpc_r')
write.csv(pheno, file="emxf2014A.csv", row.names=FALSE)

pheno<-read.csv("rgxh2012.csv")
head(pheno)

days<-c(312,333)
pheno[,1]=as.factor(pheno[,1])
pheno[,2]=as.factor(pheno[,2])

pheno[,5]=audpc(pheno[,3:4],days)
pheno[,6]=pheno[,5]/ RH_12 *100
colnames(pheno)[5]=c('audpc')
colnames(pheno)[6]=c('audpc_r')
write.csv(pheno, file="rgxh2012A.csv", row.names=FALSE)

pheno<-read.csv("rgxh2013.csv")
head(pheno)

days<-c(325,339)
pheno[,1]=as.factor(pheno[,1])
pheno[,2]=as.factor(pheno[,2])

pheno[,5]=audpc(pheno[,3:4],days)
pheno[,6]=pheno[,5]/ RH_13 *100
colnames(pheno)[5]=c('audpc')
colnames(pheno)[6]=c('audpc_r')
write.csv(pheno, file="rgxh2013A.csv", row.names=FALSE)

pheno<-read.csv("rgxh2014.csv")
head(pheno)

days<-c(281,338)
pheno[,1]=as.factor(pheno[,1])
pheno[,2]=as.factor(pheno[,2])

pheno[,5]=audpc(pheno[,3:4],days)
pheno[,6]=pheno[,5]/ RH_14 *100
colnames(pheno)[5]=c('audpc')
colnames(pheno)[6]=c('audpc_r')
write.csv(pheno, file="rgxh2014A.csv", row.names=FALSE)

pheno<-read.csv("rgxh2016.csv")
head(pheno)

days<-c(358,370,375)
pheno[,1]=as.factor(pheno[,1])
pheno[,2]=as.factor(pheno[,2])

pheno[,6]=audpc(pheno[,3:5],days)
pheno[,7]=pheno[,6]/ RH_16 *100
colnames(pheno)[6]=c('audpc')
colnames(pheno)[7]=c('audpc_r')
write.csv(pheno, file="rgxh2016A.csv", row.names=FALSE)

my_files<-list.files(pattern="*A.csv")
x=NULL
for (i in 1:length(my_files)){
  print(my_files[i])
  data<-my_files[i]
  dataR<-na.omit(read.csv(data))
  print(colnames(dataR))
  drops <- "Mil3"
  dataR<-dataR[ , !(names(dataR) %in% drops)]
  data<-gsub("A.csv","",data)
  popn<-gsub("20*","",data)
  popn<-gsub("A.csv","",popn)
  dataR$pop<- rep(popn,nrow(dataR))
  year<- gsub("emxf","",data)
  year<- gsub("rgxh","",data)
  dataR$year<- rep(year,nrow(dataR))
  x<-rbind(x,dataR)
}

name_data<-"RH_EF_relative_AUDPC.csv"
write.csv(x, name_data, row.names=TRUE)
