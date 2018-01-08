setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data")
my_files<-list.files(pattern="*A.csv")
for (i in 1:length(my_files)){
  print(i)
  data<-my_files[i]
  dataR<-na.omit(read.csv(data))
  data<-gsub("csv","",data)
  Mean<-with(dataR,tapply(audpc,Genotype,mean)) 
  Mean<-as.data.frame(Mean)
  name_data<-paste0(data,"_M.csv")
  write.csv(Mean, name_data, row.names=TRUE)
}