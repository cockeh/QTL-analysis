
setwd("/Users/helencockerton/Desktop/QTL/Files for git hub")

score<-read.csv("./data/emxflocM.csv")
for (i in 2:5){
  pheno<-colnames(score[i])
for (hg in 1:7)
{
  #subgenome
  for(sg in c("a","b","c","d"))
  {
    lg<- paste0(hg,sg) #ie 1A... 7D
    fname<- paste0("X",lg,pheno,"_rand_Pvalues")
    data<- read.csv(fname)
    assign(fname, data)
    x<-(nrow(data)) *0.05
    round(x)
    print(data[x,4])
  }}
}