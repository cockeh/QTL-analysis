# Script Number 2
# Calc BLUP for mildew data

library(nlme)

setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/")
All <- read.csv("RH_EF_relative_AUDPC.csv")
exf<-All[grep("emxf", All$pop), ]
rxh<-All[grep("rgxh", All$pop), ]
# Standardise genotype names 

rxh$Genotype<-gsub("RG","RH",rxh$Genotype)
rxh$Genotype<-gsub("Redgntlt","Redgauntlet",rxh$Genotype)

# Below taken from http://rstudio-pubs-static.s3.amazonaws.com/1856_bf7c5943ff4f4d97983e823c7a21b4a9.html

BLUP <- function(data = All, model = audpc_r ~ Genotype, random = ~1 | year, trait = "audpc_r") {
  lmeout1 <- lme(model, data = data, random = random)
  ped.hat1 <- lmeout1$coef$fixed
  ped.hat1[-1] <- ped.hat1[-1] + ped.hat1[1]
  #names(ped.hat1)[1] = "EF001"
  names(ped.hat1) <- gsub("Genotype2", "", names(ped.hat1))
  tped <- data.frame(Genotype = names(ped.hat1), trait = ped.hat1)
  names(tped)[2] <- trait
  return(tped)
}

exf<-na.omit(exf)
rxh<-na.omit(rxh)

eB <- BLUP(data = exf, model = audpc_r ~ Genotype, random = ~1 | year, trait = "mildew")
rB <- BLUP(data = rxh, model = audpc_r ~ Genotype, random = ~1 | year, trait = "mildew")

hist(eB$mildew, breaks = 30, col = "red")
hist(rB$mildew, breaks = 30, col = "red")
eB$Genotype<-gsub("Genotype","",eB$Genotype)
rB$Genotype<-gsub("Genotype","",rB$Genotype)
 
 #eB<-with(exf,tapply(audpc_r,Genotype,mean)) 
 eB<-as.data.frame(eB)
 eB<-na.omit(eB)
 colnames(eB)<-c("X","Mean")
 eB[1,1]<-"EF001"
# eE<-cbind(eM,eB)
# summary(lm(eE$eM~eE$mildew))
# plot(lm(eE$eM~eE$mildew))
# plot(eE$mildew,eE$eM)

#rB<-with(rxh,tapply(audpc_r,Genotype,mean)) 
rB<-as.data.frame(rB)
rB<-na.omit(rB)
colnames(rB)<-c("X","Mean")

write.csv(eB,"emxfcombA._M.csv", row.names=FALSE)
write.csv(rB,"rgxhcombA._M.csv", row.names=FALSE)
 