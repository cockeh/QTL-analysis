library(MASS)
library(data.table)
setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/")
myFiles2 <- list.files(pattern="hk.csv")
myFiles <- myFiles2[2:8] 

locM<-read.csv("emxflocM.csv")

V1 <- gsub("PSelect","",myFiles)
V1 <- gsub("hk.csv","",V1)
V1 <-unique(V1)

#Recode hk markers
Mk<- c("Affx.88872747","Affx.88836738","Affx.88827079","Affx.88820302","Affx.88834011","Affx.88894284","Affx.88903226","Affx.88880503","Affx.88855607","Affx.88867030","Affx.88879239","Affx.88815852","Affx.88816867","Affx.88829155","Affx.88891632")
for (marker in 1:length(Mk)){
  Mki<-Mk[marker]
  locM[Mki] <- data.frame(lapply(locM[Mki], as.character), stringsAsFactors=FALSE)
}

for (i in 1:length(locM$Affx.88889211)) {locM$Affx.88889211[i] <- if (locM$Affx.88889211[i] == "hh") locM$ Affx.88889211[i] <- "A" else "B"}
for (i in 1:length(locM$ Affx.88872747)) {locM$Affx.88872747[i] <- if (locM$Affx.88872747[i] == "hh") locM$ Affx.88872747[i] <- "A" else if (locM$Affx.88872747[i] == "kk") locM$ Affx.88872747[i] <- "A" else "B"}
for (i in 1:length(locM$ Affx.88836738)) {locM$Affx.88836738[i] <- if (locM$Affx.88836738[i] == "kh") locM$ Affx.88836738[i] <- "A" else "B"}
for (i in 1:length(locM$ Affx.88827079)) {locM$Affx.88827079[i] <- if (locM$Affx.88827079[i] == "hh") locM$ Affx.88827079[i] <- "A" else if (locM$Affx.88827079[i] == "kk") locM$ Affx.88827079[i] <- "A" else "B"}
for (i in 1:length(locM$ Affx.88820302)) {locM$Affx.88820302[i] <- if (locM$Affx.88820302[i] == "hk") locM$ Affx.88820302[i] <- "C" else if (locM$Affx.88820302[i] == "kh") locM$ Affx.88820302[i] <- "B" else "A"}
for (i in 1:length(locM$ Affx.88834011)) {locM$Affx.88834011[i] <- if (locM$Affx.88834011[i] == "hh") locM$ Affx.88834011[i] <- "A" else if (locM$Affx.88834011[i] == "kk") locM$ Affx.88834011[i] <- "B" else "C"}
for (i in 1:length(locM$ Affx.88894284)) {locM$Affx.88894284[i] <- if (locM$Affx.88894284[i] == "kk") locM$ Affx.88894284[i] <- "B" else "A"}
for (i in 1:length(locM$ Affx.88903226)) {locM$Affx.88903226[i] <- if (locM$Affx.88903226[i] == "hh") locM$ Affx.88903226[i] <- "A" else if (locM$Affx.88903226[i] == "kk") locM$ Affx.88903226[i] <- "A" else "B"}
for (i in 1:length(locM$ Affx.88880503)) {locM$Affx.88880503[i] <- if (locM$Affx.88880503[i] == "kk") locM$ Affx.88880503[i] <- "B" else "A"}
for (i in 1:length(locM$ Affx.88867030)) {locM$Affx.88867030[i] <- if (locM$Affx.88867030[i] == "hk") locM$ Affx.88867030[i] <- "B" else "A"}
for (i in 1:length(locM$ Affx.88879239)) {locM$Affx.88879239[i] <- if (locM$Affx.88879239[i] == "hh") locM$ Affx.88879239[i] <- "B" else if (locM$Affx.88879239[i] == "hk") locM$ Affx.88879239[i] <- "C" else "A"}
for (i in 1:length(locM$ Affx.88815852)) {locM$Affx.88815852[i] <- if (locM$Affx.88815852[i] == "hh") locM$ Affx.88815852[i] <- "A" else if (locM$Affx.88815852[i] == "kk") locM$ Affx.88815852[i] <- "A" else "B"}
for (i in 1:length(locM$ Affx.88816867)) {locM$Affx.88816867[i] <- if (locM$Affx.88816867[i] == "hh") locM$ Affx.88816867[i] <- "C" else if (locM$Affx.88816867[i] == "kk") locM$ Affx.88816867[i] <- "B" else "A"}
for (i in 1:length(locM$ Affx.88829155)) {locM$Affx.88829155[i] <- if (locM$Affx.88829155[i] == "hh") locM$ Affx.88829155[i] <- "C" else if (locM$Affx.88829155[i] == "kk") locM$ Affx.88829155[i] <- "B" else "A"}
for (i in 1:length(locM$ Affx.88891632)) {locM$Affx.88891632[i] <- if (locM$Affx.88891632[i] == "hh") locM$ Affx.88891632[i] <- "C" else if (locM$Affx.88891632[i] == "kk") locM$ Affx.88891632[i] <- "B" else "A"}

for (marker in 1:length(Mk)){
  Mki<-Mk[marker]
  locM[Mki] <- data.frame(lapply(locM[Mki], as.factor), stringsAsFactors=FALSE)
}

i=1
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
  
  # Creat models without each varriable in turn to calculate PRE score
  
  for(d in 1:length(y[, 1]))
  {
    qtl_modelc = paste("qtl_modelc",d,"= lm(EFcombmil~",sep="")
    for(p in 1:length(y[, 1]))
    {
      if(p==d) next 
      else
        qtl_modelc<- paste(qtl_modelc, y[p,1], '+')
      d+1
    }
    file_name1<-"EFcombmilqtl_modelc.R"
    qtl_modelc<-paste(qtl_modelc, ", data=locM)\n")
    qtl_modelc<-gsub("+ ,",",",fixed = TRUE, qtl_modelc)
    write(qtl_modelc, file=file_name1,append=TRUE)
  }
  
  source(file_name1)
  
  # create formulas to calcualte rss
  x2<-NULL 
    #paste("rss <- sum(residuals(qtl_model)^2)\n")
  for (i in 1:length(y[, 1])) {
    file_name3<-"EFcombmilrss.R"
    x<-paste('rss',i,"<- sum(residuals(qtl_modelc",i,")^2)\n", sep="")
    x2<-paste(x2,x,sep = "")
    write(x2, file=file_name3,append=TRUE)
  }
  source(file_name3)
  
  # create formulas to calcualte PRE
  x2=NULL
  for (i in 1:length(y[, 1])) {
    file_name2<-"X.R"
    x<-paste('X',i,"<- (rss",i,"-rss)/rss",i,"*100\n", sep="")
    x2<-paste(x2,x,sep = "")
    write(x2, file=file_name2,append=TRUE)
  }
  
  source(file_name2)
  
  # create list of PRE
  x<-"h=c("
  for (i in 1:length(y[, 1])) {
    x<-paste(x,'X',i,",", sep="")  
  }
  x<-paste(x,")", sep="")
  x<-gsub(",)",")",fixed = TRUE, x)
  write(x, file="h.r")
  source("h.r")
  
  QTLs<- QTL[which(QTL$Rname %in% y[,1]),]
  QTLs
  QTLs$PRE=h
  filename<-paste("Significant_QTL_hk_BF_PRE",V1a,".csv", sep="")
  write.csv(QTLs, file=filename,row.names= F)
  


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


locM$Affx.88846914 <- as.character(locM$Affx.88846914)
locM$Affx.88880944 <- as.character(locM$Affx.88880944)
for (i in 1:length(locM$Affx.88846914)) {
locM$Affx.88846914[i] <- if (locM$Affx.88846914[i] == "lm") locM$ Affx.88846914[i] <- "A" else "B"}
for (i in 1:length(locM$Affx.88880944)) {
locM$Affx.88880944[i] <- if (locM$Affx.88880944[i] == "np") locM$ Affx.88880944[i] <- "A" else "B"}
locM$Affx.88846914 <-as.factor(locM$Affx.88846914)
locM$Affx.88880944 <-as.factor(locM$Affx.88880944)

effect_model<-lm(formula = EFcombmil ~ Affx.88846914 + Affx.88879339 + Affx.88831545 + 
                   Affx.88902178 + Affx.88880944 + Affx.88862333 + Affx.88810959, 
                 data = locM)
summary(effect_model)

sumA<-summary(effect_model)
filename<-paste("Model_SummaryBF",V1a,"switched.txt", sep="")
capture.output(sumA, file = filename)

