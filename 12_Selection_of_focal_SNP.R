# Selection of focal SNP

library(data.table)

setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/")
dir()
myFiles <- list.files(pattern="Sig_values*")
notmyFiles <- list.files(pattern="*MS.csv")
myFiles <- myFiles[!myFiles %in% notmyFiles]
i=1
for(i in 1:length(myFiles)){
  name<-read.csv(myFiles[i]) 
  name$group <- paste0(name$lg,name$mtype)
  name<-data.table(name)
  sig<- data.frame(name[ , .SD[which.min(pval)], by = group])
  y<- data.frame(matrix(unlist(rep(myFiles[i],length(z)))))
  sig<- cbind(sig,y)
  sig <- sig[sig$pval <=0.05 , ]
  fn<-gsub("Sig_values","",myFiles[i])
  filenamep<-paste0("PSelect",fn)
  write.csv(sig, file=filenamep)
}
