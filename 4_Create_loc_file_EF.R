# Import AUDPC data for RGH mildew and change mis called phenotypes to NA
# Then order by population map file
# Bind genotypes together
# Add phenotypes to marker data to make locfile

setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data")
order<-read.csv("Progeny_order_EF.csv",header=FALSE)
my_files<-list.files(pattern="*_M.csv")
for (i in 1:7){
  data<-my_files[i]
  dataR<-read.csv(data)
  order$id  <- 1:nrow(order)
  dataR$X<-gsub("EF","EMxFE",dataR$X)
  dataR<-merge(dataR, order, by.x= "X", by.y= "V1", all.y=T)
  dataR<-dataR[order(dataR$id), ]
  dataR<-dataR[,c(1,2)]
  colnames(dataR)<- c("progeny", "score")
  data<-gsub("A._M.csv","",data)
  file_name<-paste0(data,"T2.csv")
  write.csv(dataR, file_name, row.names=FALSE)
}
x<-NULL
my_files<-list.files(pattern="*T2.csv")
for (i in 1:length(my_files)){
  data<-my_files[i]
  dataR<-read.csv(data)
  x<-append(x,dataR[2])
}

ff<-cbind(dataR[1],as.data.frame(x))

colnames(ff)[c(2:8)]<-c("mil11","mil12D",	"mil12",	"mil13D",	"mil13",	"mil14",	"EFcombmil") 

markers<-read.csv("emxfloc.csv")
ff$id  <- 1:nrow(ff)
dataM<-merge(ff, markers, by.x= "progeny", by.y= "individual", all.y=T)
dataM<-dataM[order(dataM$id), ]
dataM<-dataM[,-9]
file_name<-"emxflocM.csv"
write.csv(dataM, file_name, row.names=FALSE)
