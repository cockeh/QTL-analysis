#!/usr/local/bin/Rscript

# Permutation test 
# cd /Users/helencockerton/Desktop/QTL/Files\ for\ git\ hub/
# for i in $(seq 1 1000); do
# echo "Itteration - $i"
# ./9_Permutation_test.R --V 1a 
# done
# setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/")

# hg=1
# sg="A"

library(optparse)

opt_list = list(
  make_option("--V", type="character", help="linkage_group")
)
opt = parse_args(OptionParser(option_list=opt_list))
lg = opt$V

score<-read.csv("./data/emxflocM.csv")
basedirectory<-"./data/exf_Hap/"
  
for (i in c(2:5)){
#     lg    <- paste0(hg,sg) #ie 1A... 7D
    fname <- paste0(basedirectory,lg,"_genotypes_cult.csv")
    scorei<-unlist((score[i]))
    sc<-colnames(score[i])
    score2<-sample(scorei,replace = TRUE)
    score2<-as.data.frame(score2)
    all_data <- read.csv(fname,row.names=1)
    mod_data<- all_data[-c(1:5)]
    mod_data<-t(mod_data)
    mod_data <- data.frame(mod_data)
    mod_data<-cbind(mod_data,score2)
df <- data.frame()
name.of.function <- function(x) {
  p<-kruskal.test(mod_data$score2~x)[[3]]
}
df<-data.frame(lapply(mod_data,name.of.function))
drops <- "score2"
df<-df[ , !(names(df) %in% drops)]
df<-t(df)
filename<-paste0("X",lg,sc,"s_rand_Pvalues") 
write.table(df,file=filename,col.names=F,append=T,sep=",")
}
