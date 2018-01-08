library(MASS)
library(data.table)
setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/")
myFiles2 <- list.files(pattern="hk.csv")
myFiles <- myFiles2[2:8] 

locM<-read.csv("emxflocM.csv")

V1 <- gsub("PSelect","",myFiles)
V1 <- gsub("hk.csv","",V1)
V1 <-unique(V1)

#Recode hk markers based on phenotype - asses using aov
# Mk<- c("marker.1","marker.2","marker.3")
# for (marker in 1:length(Mk)){
#   Mki<-Mk[marker]
#   locM[Mki] <- data.frame(lapply(locM[Mki], as.character), stringsAsFactors=FALSE)
# }
# 
# for (i in 1:length(locM$marker.1)) {locM$marker.1[i] <- if (locM$marker.1[i] == "hh") locM$ marker.1[i] <- "A" else "B"}
# for (i in 1:length(locM$marker.2)) {locM$marker.2[i] <- if (locM$marker.2[i] == "hh") locM$ marker.2[i] <- "A" else if (locM$marker.2[i] == "kk") locM$ marker.2[i] <- "A" else "B"}

# for (marker in 1:length(Mk)){
#   Mki<-Mk[marker]
#   locM[Mki] <- data.frame(lapply(locM[Mki], as.factor), stringsAsFactors=FALSE)
# }


for(i in 1:length(myFiles)){
  V1a<-V1[i]
  QTL<-read.csv(myFiles[i], header =T) 
  QTL$Rname<- gsub('-', '.', QTL$Rname)
  sig_markers<-QTL$Rname
  marker_names <- paste0(sig_markers,collapse='+')
  fullformula <- as.formula(paste0(V1,"~",marker_names))
  qtl_model <- stepAIC(lm(formula=fullformula,data = locM))
  qtl_model <- lm(formula=fullformula,data = locM)
  step_model <- stepAIC(qtl_model,
                        scope=list(upper=fullformula,lower=as.formula(score~1)),
                        direction="both",
                        steps=9999)
  summary(qtl_model)
  
  
  while(TRUE)
  {
    #extract the chosen markers and their pvalues
    df <- as.data.frame(coef(summary(step_model))[-1,])
    df$marker <- rownames(df)
    df$pvalue <- df[,4]
    df <- df[,c("marker","pvalue")]
    df <- df[order(-df$pvalue),]
    rownames(df) <- 1:nrow(df)
    sig_markers <- df$marker
    
    print(paste0("largest p value=",df$pvalue[1]))
    print(paste0("number of markers=",length(df$marker)))
    
    #stop if least significant marker has pvalue <= 0.05
    BF<-0.05/nrow(QTL)
    if(df$pvalue[1] <= BF) break
    
    #exclude least significant marker
    markerpref<-c("kk","hk","kh","hh","np","lm","B","C")
    for (i in markerpref){df$marker <- lapply(df$marker, function(x) {
      gsub(i, "", x) 
    })}
    for (i in markerpref){sig_markers <- lapply(sig_markers, function(x) {
      gsub(i, "", x) 
    })}
    
    sig_markers <- as.data.frame(sig_markers)
    sig_markers <- sig_markers[sig_markers!=df$marker[1]]
    
    #run step again
    marker_names <- paste0(sig_markers,collapse='+')
    fullformula <- as.formula(paste0(V1a,"~",marker_names))
    qtl_model <- lm(formula=fullformula,data=locM)
    step_model <- stepAIC(qtl_model,
                          scope=list(upper=fullformula,lower=as.formula(score~1)),
                          direction="both",
                          steps=9999)
  }
  
  sumA<-summary(step_model)
  filename<-paste("Model_Summary_BF",V1a,".txt", sep="")
  capture.output(sumA, file = filename)
  
  summary(step_model)
  markerpref<-c("lm","np","hk","kh","kk","B","C")
  y <- variable.names(step_model)
  y<- as.data.frame(y)
  for (i in markerpref){y <- lapply(y, function(x) {
    gsub(i, "", x) 
  })
  }
  y<- as.data.frame(y)
  y<- unique(y)
  y <- y[-1,]
  y<- as.data.frame(y)
  y
  
  QTLs<- QTL[which(QTL$Rname %in% y[,1]),]
  QTLs
  
  filename<-paste("Significant_QTL_hk_BF",V1a,".csv", sep="")
  write.csv(QTLs, file=filename,row.names= F)
  
}

#Combine significant QTL 


mil11<-read.csv("Significant_QTL_hk_BFmil11.csv", header=TRUE)
mil12<-read.csv("Significant_QTL_hk_BFmil12.csv", header=TRUE)
mil12D<-read.csv("Significant_QTL_hk_BFmil12D.csv", header=TRUE)
mil13<-read.csv("Significant_QTL_hk_BFmil13.csv", header=TRUE)
mil13D<-read.csv("Significant_QTL_hk_BFmil13D.csv", header=TRUE)
mil14<-read.csv("Significant_QTL_hk_BFmil14.csv", header=TRUE)
allmil<-rbind(mil11,mil12,mil12D,mil13,mil13D,mil14)
write.csv(allmil, file="All_Significant_QTL_hk_EF_BF.csv")

## CHange direction to calc effect size
# 
# 
# locM$marker.1 <- as.character(locM$marker.1)
# locM$marker.2 <- as.character(locM$marker.2)
# for (i in 1:length(locM$marker.21)) {
# locM$marker.1[i] <- if (locM$marker.1[i] == "lm") locM$ marker.1[i] <- "A" else "B"}
# for (i in 1:length(locM$marker.2)) {
# locM$marker.2[i] <- if (locM$marker.2[i] == "np") locM$ marker.2[i] <- "A" else "B"}
# locM$marker.21 <-as.factor(locM$marker.1)
# locM$marker.2 <-as.factor(locM$marker.2)
# effect_model<-lm(formula = EFcombmil ~ "marker.1","marker.2","marker.3", 
#                  data = locM)
# summary(effect_model)

