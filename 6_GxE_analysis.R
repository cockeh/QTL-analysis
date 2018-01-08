setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/")
library(MASS)
library(lme4)
install.packages("lmerTest")
library(lmerTest)

All <- read.csv("RH_EF_relative_AUDPC.csv")
drops <- c("RH011","RH007","RH023","RH011","RH017","RH023","RH064","RH072","RH126","RH015","RH044","RH046","RH099","RH166","RH171","RH183","RH040","RH052","Emily","Fenella","Hapil","Elsanta","Sonata","Redgauntlet","Redgntlt")
for(i in 1:length(drops)){
  All$Genotype[All$Genotype == drops[i]] <- NA
}
All <-All[!is.na(All$audpc_r),]
write.csv(All,"RH_EF_relative_AUDPC2.csv")
exf<-All[grep("emxf", All$pop), ]
rxh<-All[grep("rgxh", All$pop), ]

# Test to see if the residuals are normally distrubuted

x<-resid(lm(audpc_r~Genotype*year, data=exf))
hist(x)
#too many data point for shaprio.test()
m <- mean(x) 
s <- sd(x) 
ks.test(x,"pnorm",m,s)
# One-sample Kolmogorov-Smirnov test
# 
# data:  x
# D = 0.097428, p-value < 2.2e-16
# alternative hypothesis: two-sided

x<-resid(lm(audpc_r~Genotype*year, data=rxh))
hist(x)
m <- mean(x) 
s <- sd(x) 
ks.test(x,"pnorm",m,s)

# One-sample Kolmogorov-Smirnov test
# 
# data:  x
# D = 0.16145,, p-value = 2.33e-12
# alternative hypothesis: two-sided

### REML and compare varrience
exf <-na.omit(exf)
l1 <- lmer(exf$audpc_r ~  1+ (1|exf$Genotype) + (1|exf$year) + (1|exf$Genotype:exf$year))
l2 <- lmer(exf$audpc_r ~  1+ (1|exf$Genotype) + (1|exf$year))
l3 <- lmer(exf$audpc_r ~  1+ (1|exf$year) + (1|exf$Genotype:exf$year))
l4 <- lmer(exf$audpc_r ~  1+ (1|exf$Genotype) + (1|exf$Genotype:exf$year))
print(anova(l1,l2))
print(anova(l1,l3))
print(anova(l1,l4))
summary(l1)
summary(l2)
summary(l3)
summary(l4)

l1 <- lmer(rxh$audpc_r ~  1+ (1|rxh$Genotype) + (1|rxh$year) + (1|rxh$Genotype:rxh$year))
l2 <- lmer(rxh$audpc_r ~  1+ (1|rxh$Genotype) + (1|rxh$year))
l3 <- lmer(rxh$audpc_r ~  1+ (1|rxh$year) + (1|rxh$Genotype:rxh$year))
l4 <- lmer(rxh$audpc_r ~  1+ (1|rxh$Genotype) + (1|rxh$Genotype:rxh$year))
print(anova(l1,l2))
print(anova(l1,l3))
print(anova(l1,l4))
summary(l1)
summary(l2)
summary(l3)
summary(l4)

