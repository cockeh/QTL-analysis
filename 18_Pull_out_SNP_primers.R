setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/")
markersef<-read.csv("Significant_QTL_hk_BFEFcombmil.csv")
EF<-markersef$Rname
markersrg<-read.csv("Significant_QTL_hk_BFRHcombmil.csv")
RH<-markersrg$Rname
All<-unlist(list(EF,RH))
primers<-read.table("../bedtools/istraw90_vesca_v1.1_snp_positions.gff3")
All<- gsub("Affx\\.","ID=Affx-",All)
All<-as.data.frame(All)
Correct_names<- merge(primers, All, by.x = "V9", by.y = "All", all.y=T)
Correct_names<-Correct_names[c(2,3,4,5,6,7,8,9,1)]
write.table(Correct_names,"mil_primers.gff",row.names=FALSE,col.names=FALSE,sep="\t", quote = FALSE)
