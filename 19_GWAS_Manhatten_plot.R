
setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/")
GWAS<- read.csv("No_Mil2.txt_istraw90.out_fix_min10_pheno_0.2_strat.assoc.linear.csv")
position_of_markers<-read.csv("vesca2consensus_map_integerposns_2017-07-09.csv", header=FALSE)

GWAS_man<-merge(GWAS,position_of_markers, by.x="SNP",by.y="V1", all.x=TRUE) 
GWAS_man <- GWAS_man[order(GWAS_man$V3),]
GWAS_man<- GWAS_man[order(GWAS_man$V2),]

GWAS_man<-GWAS_man[!is.na(GWAS_man$V3), ]
GWAS_man_NA<-GWAS_man[is.na(GWAS_man$V2), ]

PosVal<- data.frame(Pv = numeric())
for (i in unique(GWAS_man$V2)) {
  x <- tail(GWAS_man[GWAS_man$V2 == i,], n=1)
  z<- x$V3 
  PosVal[i,] = z   
}
setDT(PosVal, keep.rownames = TRUE)[]
newrow<- NA
PosVal <- PosVal[order(PosVal$rn),]
PosVal = rbind(newrow,PosVal,fill=TRUE)
nc  <- nrow(PosVal)
PosVal[1, 3] <- 0

shift <- function(x, n){
  c(x[-(seq(n))])
}
PosVal<- rbind(PosVal, "29"=NA, fill=TRUE)
PosVal$rn <- shift(PosVal$rn, 1)
#error message is not an issue
PosVal$x <- NULL
PosVal <- PosVal[-30, ]
PosVal$lg <- PosVal$rn
PosVal$Pv2<-cumsum(PosVal$Pv)

PosVal$lty<- rep(3,nrow(PosVal))
thick<- c(1,5,9,13,17,21,25,29)
PosVal$lty[thick]<-gsub("3","1",PosVal$lty[thick])

plyr1 <- merge(GWAS_man, PosVal, by.x = "V2",by.y ="lg")
plyr1$position <- plyr1$V3 + plyr1$Pv2
threshold<-0.05
threshold2<-0.01

for (i in 1:nrow(PosVal)){
  i2<-i+1
  textpos<- PosVal[i,4]+((PosVal[i2,4]-PosVal[i,4])/2)
  print(textpos)
  PosVal$V5[i]<-textpos
}
print(unlist(PosVal$V5))

octoploid <- ggplot(plyr1, aes(x=position, y=-log10(UNADJ))) +
  geom_point(pch=1)+
  #ggtitle("Octoploid")+ 
  geom_vline(data=PosVal,aes(xintercept=Pv2,linetype=lty))+
  # scale_y_log10(limits = c(0.0000001,1))+
  geom_hline(aes(yintercept=3.823909))+
  annotate("text", x = unlist(PosVal$V5[1:28]), y = -0.2, label = c("1A","1B","1C","1D","2A","2B","2C","2D","3A","3B","3C","3D","4A","4B","4C","4D","5A","5B","5C","5D","6A","6B","6C","6D","7A","7B","7C","7D"))+
  scale_x_continuous(expand = c(0, 0))+
  coord_cartesian(xlim=c(0,B29),ylim=c(-0.5, 5.5))+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'), axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x="Position (Mb)",y=expression(-log[10](italic(p))))+
  theme(legend.position='none')

# snp_position_vesca<- read.csv("istraw90_vesca_v1.1_snp_positions.csv",header=FALSE)
# 
# GWAS_man<-merge(GWAS,snp_position_vesca, by.x="SNP",by.y="V9", all.x=TRUE) 
# GWAS_man <- GWAS_man[order(GWAS_man$V4),]
# GWAS_man<- GWAS_man[order(GWAS_man$V1),]
# 
# GWAS_man<-GWAS_man[!is.na(GWAS_man$V4), ]
# 
# PosVal<- data.frame(Pv = numeric())
# for (i in unique(GWAS_man$V1)) {
#   x <- tail(GWAS_man[GWAS_man$V1 == i,], n=1)
#   z<- x$V4 
#   PosVal[i,] = z   
# }
# setDT(PosVal, keep.rownames = TRUE)[]
# newrow<- NA
# PosVal <- PosVal[order(PosVal$rn),]
# PosVal = rbind(newrow,PosVal,fill=TRUE)
# nc  <- nrow(PosVal)
# PosVal[1, 3] <- 0
# 
# shift <- function(x, n){
#   c(x[-(seq(n))])
# }
# PosVal<- rbind(PosVal, "9"=NA, fill=TRUE)
# PosVal$rn <- shift(PosVal$rn, 1)
# #error message is not an issue
# PosVal$x <- NULL
# PosVal <- PosVal[-9, ]
# PosVal$lg <- PosVal$rn
# PosVal$Pv2<-cumsum(PosVal$Pv)
# 
# plyr1 <- merge(GWAS_man, PosVal, by.x = "V1",by.y ="lg")
# plyr1$position <- plyr1$V4 + plyr1$Pv2
# threshold<-0.05
# threshold2<-0.01
# 
# for (i in 1:nrow(PosVal)){
#   i2<-i+1
#   textpos<- PosVal[i,4]+((PosVal[i2,4]-PosVal[i,4])/2)
#   print(textpos)
#   PosVal$V5[i]<-textpos
# }
# print(unlist(PosVal$V5))
# 
# vesca <- ggplot(plyr1, aes(x=position, y=-log10(UNADJ))) +
#   geom_point(pch=1)+
#   ggtitle("Vesca")+ 
#   geom_vline(data=PosVal,aes(xintercept=Pv2))+
#   # scale_y_log10(limits = c(0.0000001,1))+
#   geom_hline(aes(yintercept=3.823909))+
#   annotate("text", x = unlist(PosVal$V5[1:7]), y = -0.2, label = c("1","2","3","4","5","6","7"))+
#   scale_x_continuous(expand = c(0, 0))+
#   coord_cartesian(xlim=c(0,194954649),ylim=c(-0.5, 5.5))+
#   theme(panel.background = element_rect(fill = 'white', colour = 'black'), axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())+
#   labs(x="Position (Mb)",y=expression(-log[10](italic(p))))+
#   theme(legend.position='none')

# Does not make sense to display vesca SNPs as all sig markers on LG 4 in vesca map to 6C in octopliod  
setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/figs/")
name<-"Fig_6.pdf"
pdf(name,width=12,height=4.3,paper='special') 
#multiplot(octoploid, vesca, cols=1)
octoploid
dev.off()