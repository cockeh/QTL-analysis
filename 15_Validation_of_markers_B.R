# Validation of markers (REML) code


library(data.table)
library(plyr)
library(ggplot2)
library(gplots)
setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/")
mil<-read.csv("mildew_AUDPC.csv")
geno<-read.csv("found_geno_name.csv")
geno$probeset_id<-gsub("PHR-","AX",geno$probeset_id)
geno$probeset_id<-gsub("NMH-","AX",geno$probeset_id)
geno$probeset_id<-gsub("MHR-","AX",geno$probeset_id)
decode<-read.csv("id_decode.csv")
Sef<-read.csv("PSelectMSEFcombmilMS_Eqiv.csv")
Sef$Pop <- rep("EXF",nrow(Sef))
Sef<-Sef[,-c(6,7)]
Srh<-read.csv("PSelectMSRHcombmilMS_Eqiv.csv")
Srh$Pop <- rep("RXH",nrow(Srh))
markers <- rbind(Sef,Srh)
decode$snp_id_90<-gsub("-",".",decode$snp_id_90)
markers$Rname<-gsub("-",".",markers$Rname)
colnames(decode)[4] <-"Rname"
marker_new<- merge(decode, markers, by.x = "snp_id_90", by.y ="Rname", all.y=T)
marker_new$probe_id_35<-gsub("-","",marker_new$probe_id_35)
marker<-marker_new$probe_id_35
ty<-t(geno)
colnames(ty)<-ty[1,]
ty<-ty[-1,]
ty<-as.data.frame(ty)
#row.names(ty)<-pheno$V1

locM<-merge(mil,ty,by.x="Genotype",by.y=0, all.x=T)
locM<-locM[!is.na(locM$audpc),]
locM[locM== -1] <- "NA"

write.csv(locM,"locM_val.csv", row.names=FALSE)

marker_pos<-data.frame()
for (i in c(1:nrow(marker_new))){
y<-as.data.frame(grep(marker_new$probe_id_35[i], colnames(locM)))
print(y)
marker_pos<-rbind(marker_pos,y)
}

library(lme4)
library(lmertest)
for(i in c(1,2,4:nrow(marker_pos))){
mp<-marker_pos[i,1]
locM2<-locM[!is.na(locM[mp]),]
locM2[mp]<-as.factor(unlist(locM2[mp]))
l1 <- lmer(locM2$audpc ~ unlist(locM2[mp]) + (1|locM2$Genotype) + (1|locM2$Block))
l2 <- lmer(locM2$audpc ~ 1 + (1|locM2$Genotype) + (1|locM2$Block),REML=FALSE)
print(i)
print(colnames(locM2[mp]))
print(anova(l1,l2)[8])
}
## 6,9
mp<-marker_pos[c(6),1]

M1<-locM[,c(1,7,13117)]
M1[c(2:4)] <- if (M1[c(2:4)] == 1) M1[c(2:4)] <- "AB" else if (M1[j,i] == 0) M1[j,i] <- "AA" else if (M1[j,i] == 2) M1[j,i] <- "BB" else "NA"
par(mfrow=c(2,3))
for (i in 3: ncol(M1)){
Mean<-with(M1,tapply(audpc,M1[i],mean)) 
Std<-as.data.frame(with(M1,tapply(audpc,M1[i],sd)))
Std<-na.omit(Std)
Num<- as.data.frame(na.omit(count(M1[i]), vars = 1))
All<-cbind(Std,Num)
All$SE<-All[1] / (sqrt(All$freq))
# g <- ggplot(M1, aes(group=M1[i],x=M1$audpc)) + 
#   geom_bar(stat="identity")
SE<- as.numeric(unlist(All[4]))
barplot2(Mean,plot.ci=TRUE,ci.l=Mean-SE,ci.u=Mean+SE, axis.lty=1, main = colnames(M1)[i], border = "black", cex.names = par("cex.axis"), col= c("azure3"), plot.grid= TRUE, ylab="AUDPC") 
box()
}

M1<-locM[,c(1,7,13117)]
Mean<-with(M1,tapply(audpc,marker.1,mean)) 
Mean<-as.data.frame(Mean)
rownames(Mean)<-c("AB","BB")
Std<-with(M1,tapply(audpc,marker.1,sd))
Std<-as.data.frame(Std)
Num<- count(M1$marker.1, vars = 1)
Num<-Num[-3,]
SE<-cbind(Std,Num)
SE$SE<-as.numeric(SE$Std) / (sqrt(SE$freq))
SE<-SE[-4,]
barplot2(Mean[,1],plot.ci=TRUE,ylim=c(0,50),ci.l=Mean[,1]-SE$SE,ci.u=Mean[,1]+SE$SE, axis.lty=1, main = "marker.1",axisnames=FALSE, border = "black", col= c("azure3"), plot.grid= TRUE, ylab="AUDPC") 
box()
text(0.5, -10, "AB", cex = 1.0)

# Just for significant marker

mod<-lm(audpc~marker.1,data=M1)
summary(mod)

#  I don’t think it makes sense to include other QTL as cofactors in the validation REML analysis because no interaction was observed at the population level. 


