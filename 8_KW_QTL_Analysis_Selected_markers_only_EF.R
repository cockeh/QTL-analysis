library(ggplot2)

setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/")
pheno<-read.csv("emxflocM.csv")
basedir<-"./exf_Hap/"
for (j in 2:8){
score <- pheno[j]
df = data.frame(lg=character(),Rname=character(),pos=double(),phase=character(),pval=double(),k=double(),sig=character(),mtype=character(),ll=double(),lm=double(),
                nn=double(),np=double(),hh=double(),hk=double(),kh=double(),kk=double(),
                stringsAsFactors=FALSE)

for (hg in 1:7)
{
  #subgenome
  for(sg in c("A","B","C","D"))
  {
    lg    <- paste0(hg,sg) #ie 1A... 7D
    fname <- paste0(basedir,"MS35",lg,"_genotypes_cult.csv")
    print(lg)
    
    #load genotype data
    all_data <- read.csv(fname,row.names=1)
    firstcol <- 6
    n_cols   <- length(colnames(all_data))
    
    pos   <- as.numeric(all_data$cM) #marker centimorgan positions on this LG
    Rname <- row.names(all_data)
    phase <-all_data$phase
    #for each marker
    i=0
    for (snp in rownames(all_data))
    {
      i <- as.numeric(i + 1)
      mtype <- all_data[snp,2] #lmxll, nnxnp or hkxhk
      genotypes <- as.data.frame(t(all_data[snp,firstcol:n_cols]))
      colnames(genotypes)[1] <- 'call'
      genotypes$call <- as.factor(genotypes$call)
      both <- cbind(genotypes,score)
      colnames(both)[2] <- 'score'
      p <- kruskal.test(score ~ call, data = both)[[3]] #get just the p value from the Kruskal Wallis test
      k <- kruskal.test(score ~ call, data = both)[[1]] 
      bycall <- aggregate(both$score,list(both$call),mean,na.rm=TRUE)
      colnames(bycall) <- c("call","mean_score")
      if(mtype == "<lmxll>") mtype <- "lmxll" #marker type
      if(mtype == "<nnxnp>") mtype <- "nnxnp"
      if(mtype == "<hkxhk>") mtype <- "hkxhk"
      p<-as.numeric(p)
      
      sig <- if(0.1>p && 0.05<p){"."} else if (0.5>p && 0.01<p){"*"} else if (0.01>p && 0.001<p){"**"} else if (0.001>p && 0.0001<p){"***"} else if (0.0001>p && 0.00001<p){"****"} else if (0.00001>p && 0.000001<p){"*****"} else if (0.000001>p && 0.0000001<p){"******"} else if (0.0000001>p){"*******"} else {""}     
      x<-c("A","B")
      y<-c("A","B")
      df2<-rbind(x,y)
      colnames(df2) <- c("call","mean_score")
      bycall<-rbind(bycall,df2)
      
      ll<-if (bycall[1,1] == "ll") {as.numeric(bycall[1,2])} else {NA}
      lm<-if (bycall[2,1] == "lm") {as.numeric(bycall[2,2])} else {NA}
      nn<-if (bycall[1,1] == "nn") {as.numeric(bycall[1,2])} else {NA}
      np<-if (bycall[2,1] == "np") {as.numeric(bycall[2,2])} else {NA}
      hh<-if (bycall[1,1] == "hh") {as.numeric(bycall[1,2])} else {NA}
      hk<-if (bycall[2,1] == "hk") {as.numeric(bycall[2,2])} else {NA}
      kh<-if (bycall[3,1] == "kh") {as.numeric(bycall[3,2])} else {NA}
      kk<-if (bycall[4,1] == "kk") {as.numeric(bycall[4,2])} else {NA}
      
      df[nrow(df) + 1, ] <- c(lg,Rname[[i]],pos[[i]],phase[[i]],p,k,sig,mtype,ll,lm,nn,np,hh,hk,kh,kk)
    }
  }
}


#convert string to factors
df$mtype = as.factor(df$mtype)

df$pos = as.numeric(df$pos)
df$pval = as.numeric(df$pval)
V1<-names(pheno[j])
name<-paste("Sig_values",V1,"MS.csv",sep="")
write.csv(df,file=name)
}
