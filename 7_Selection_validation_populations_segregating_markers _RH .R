setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/")
decode<-read.csv("found_geno_name.csv")
decode$probeset_id<-gsub("PHR-","AX",decode$probeset_id)
decode$probeset_id<-gsub("NMH-","AX",decode$probeset_id)
decode$probeset_id<-gsub("MHR-","AX",decode$probeset_id)
recode<-read.csv("id_decode.csv")
recode$probe_id_35<-gsub("AX-","AX",recode$probe_id_35)
mergecode<-merge(decode,recode, by.x= "probeset_id", by.y="probe_id_35", all.x =TRUE)
setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/rxh_Hap/")

my_files<-list.files(pattern="35*")
for (i in 1:30){
  file_name<-my_files[i]
  LG<-read.csv(file_name)
  Validation<- LG[which(LG$V1 %in% mergecode$snp_id_90),]
  print(my_files[i])
  print(nrow(Validation))
  out_file_name<-paste0("MS",file_name)
  write.csv(Validation,out_file_name, row.names = F)
}

