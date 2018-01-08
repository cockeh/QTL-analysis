# Import AUDPC data for RGH mildew and change mis called phenotypes to NA
# Then order by population map file
# Bind genotypes together
# Add phenotypes to marker data to make locfile

setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data")
order<-read.csv("Progeny_order.csv",header=FALSE)
my_files<-list.files(pattern="*_M.csv")
for (i in 8:length(my_files)){
  data<-my_files[i]
  dataR<-read.csv(data)
  dataR$X<- gsub("RG","RGxHA",dataR$X)
  dataR$X<- gsub("RH","RGxHA",dataR$X)
  drops <- c("RGxHA1011","RGxHA023","RGxHA072","RGxHA126","Hapil","Elsanta","Sonata","Redgauntlet","Redgntlt")
  for(i in 1:length(drops)){
    dataR$Mean[dataR$X == drops[i]] <- NA
  }
  order$id  <- 1:nrow(order)
  dataR<-merge(dataR, order, by.x= "X", by.y= "V2", all.y=T)
  dataR<-dataR[order(dataR$id), ]
  dataR<-dataR[,c(3,2)]
  colnames(dataR)<- c("progeny", "score")
  data<-gsub("A._M.csv","",data)
  file_name<-paste0(data,"T.csv")
  write.csv(dataR, file_name, row.names=FALSE)
}
x<-NULL
my_files<-list.files(pattern="*T.csv")
for (i in 1:length(my_files)){
  data<-my_files[i]
  dataR<-read.csv(data)
  x<-append(x,dataR[2])
}

ff<-cbind(dataR[1],as.data.frame(x))

colnames(ff)[c(2:6)]<-c("RGH_mil_12","RGH_mil_13","RGH_mil_14","RGH_mil_16", "RHcombmil") 

markers<-read.csv("rgxhloc.csv")
ff$id  <- 1:nrow(ff)
dataM<-merge(ff, markers, by.x= "progeny", by.y= "individual", all.y=T)
dataM<-dataM[order(dataM$id), ]
dataM<-dataM[,-7]
file_name<-"rgxhlocM.csv"
write.csv(dataM, file_name, row.names=FALSE)
