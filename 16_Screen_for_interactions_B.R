setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/")
Srh<-read.csv("Significant_QTL_hk_BFRHcombmil.csv")
Sef<-read.csv("Significant_QTL_hk_BFEFcombmil.csv")
rxhloc<-read.csv("rgxhlocM.csv")
exfloc<-read.csv("emxflocM.csv")

Srh$Rname
Sef$Rname

rh<-aov(RHcombmil ~ ("marker.1","marker.2","marker.3")^2, data=rxhloc)
summary(rh)

# No interactions

ef<-aov(EFcombmil ~ ("marker.1","marker.2","marker.3")^2, data=exfloc)
summary(ef)
        
# No interactions

