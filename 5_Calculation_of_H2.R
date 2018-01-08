# Broadsense heritability 

library(lme4)

setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/")
myfiles<- list.files(pattern="*A.csv")
for(i in 1:length(myfiles)){
Mildew_pheno<- read.csv(myfiles[i], header=TRUE)
drops <- c("RH011","RH007","RH023","RH011","RH017","RH023","RH064","RH072","RH126","RH015","RH044","RH046","RH099","RH166","RH171","RH183","RH040","RH052","Emily","Fenella","Hapil","Elsanta","Sonata","Redgauntlet","Redgntlt")
for(i2 in 1:length(drops)){
  Mildew_pheno$Genotype[Mildew_pheno$Genotype == drops[i]] <- NA
}
Mildew_pheno<-na.omit(Mildew_pheno)
lmfit =lmer(Mildew_pheno$audpc ~1 +(1|Mildew_pheno$Genotype), na.action = na.exclude)
print(myfiles[i])
print(summary(lmfit))
print("H = Vg^2/ (Vg^2 +Ve^2)")
}

H16<-(57.34^2)/((57.34^2+77.49^2))
H14r<-737.6^2/(737.6^2+657.5^2)
H13r<-66.11^2/(66.11^2+102.83^2)
H12r<-111.6^2/(111.6^2+103.2^2)
H14e<-101.7^2/(101.7^2+112.5^2)
H13De<-575.9^2/(575.9^2+ 699.7^2)
H13e<-52.56^2/(52.56^2+55.35^2)
H12De<-12.88^2/(12.88^2+20.90^2)
H12e<-9.178^2/(9.178^2+18.009^2)
H11e<-32.05 ^2/(32.05 ^2+19.42^2)

H16
H14r
H13r
H12r
H14e
H13De
H13e
H12De
H12e
H11e
# [1] "emxf2011A.csv"
# summary from lme4 is returned
# some computational error has occurred in lmerTest
# Linear mixed model fit by REML ['lmerMod']
# Formula: Mildew_pheno$audpc ~ 1 + (1 | Mildew_pheno$Genotype)
# 
# REML criterion at convergence: 3198.9
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -2.58457 -0.48220  0.08418  0.53898  2.34970 
# 
# Random effects:
#   Groups                Name        Variance Std.Dev.
# Mildew_pheno$Genotype (Intercept) 32.05    5.661   
# Residual                          19.42    4.407   
# Number of obs: 500, groups:  Mildew_pheno$Genotype, 168
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  33.7925     0.4794   70.49
# [1] "H = Vg^2/ (Vg^2 +Ve^2)"
# [1] "emxf2012A.csv"
# summary from lme4 is returned
# some computational error has occurred in lmerTest
# Linear mixed model fit by REML ['lmerMod']
# Formula: Mildew_pheno$audpc ~ 1 + (1 | Mildew_pheno$Genotype)
# 
# REML criterion at convergence: 2448.3
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.9987 -0.4525  0.1825  0.6018  2.3033 
# 
# Random effects:
#   Groups                Name        Variance Std.Dev.
# Mildew_pheno$Genotype (Intercept)  9.178   3.030   
# Residual                          18.009   4.244   
# Number of obs: 404, groups:  Mildew_pheno$Genotype, 180
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  25.0414     0.3138   79.81
# [1] "H = Vg^2/ (Vg^2 +Ve^2)"
# [1] "emxf2012DA.csv"
# summary from lme4 is returned
# some computational error has occurred in lmerTest
# Linear mixed model fit by REML ['lmerMod']
# Formula: Mildew_pheno$audpc ~ 1 + (1 | Mildew_pheno$Genotype)
# 
# REML criterion at convergence: 5882.7
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -2.55304 -0.62072 -0.04567  0.53745  2.96367 
# 
# Random effects:
#   Groups                Name        Variance Std.Dev.
# Mildew_pheno$Genotype (Intercept) 12.88    3.588   
# Residual                          20.90    4.571   
# Number of obs: 955, groups:  Mildew_pheno$Genotype, 194
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  19.5335     0.2973   65.71
# [1] "H = Vg^2/ (Vg^2 +Ve^2)"
# [1] "emxf2013A.csv"
# summary from lme4 is returned
# some computational error has occurred in lmerTest
# Linear mixed model fit by REML ['lmerMod']
# Formula: Mildew_pheno$audpc ~ 1 + (1 | Mildew_pheno$Genotype)
# 
# REML criterion at convergence: 4478.3
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.6625 -0.5946 -0.0360  0.5325  3.3551 
# 
# Random effects:
#   Groups                Name        Variance Std.Dev.
# Mildew_pheno$Genotype (Intercept) 52.56    7.25    
# Residual                          55.35    7.44    
# Number of obs: 618, groups:  Mildew_pheno$Genotype, 160
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  29.2844     0.6472   45.24
# [1] "H = Vg^2/ (Vg^2 +Ve^2)"
# [1] "emxf2013DA.csv"
# summary from lme4 is returned
# some computational error has occurred in lmerTest
# Linear mixed model fit by REML ['lmerMod']
# Formula: Mildew_pheno$audpc ~ 1 + (1 | Mildew_pheno$Genotype)
# 
# REML criterion at convergence: 3820.9
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -2.74377 -0.59786 -0.01747  0.56292  2.49408 
# 
# Random effects:
#   Groups                Name        Variance Std.Dev.
# Mildew_pheno$Genotype (Intercept) 101.7    10.08   
# Residual                          112.5    10.61   
# Number of obs: 478, groups:  Mildew_pheno$Genotype, 179
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  60.6878     0.9203   65.95
# [1] "H = Vg^2/ (Vg^2 +Ve^2)"
# [1] "emxf2014A.csv"
# summary from lme4 is returned
# some computational error has occurred in lmerTest
# Linear mixed model fit by REML ['lmerMod']
# Formula: Mildew_pheno$audpc ~ 1 + (1 | Mildew_pheno$Genotype)
# 
# REML criterion at convergence: 6064.9
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -1.7456 -0.6681 -0.2549  0.4532  3.1840 
# 
# Random effects:
#   Groups                Name        Variance Std.Dev.
# Mildew_pheno$Genotype (Intercept) 575.9    24.00   
# Residual                          699.7    26.45   
# Number of obs: 622, groups:  Mildew_pheno$Genotype, 161
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)   90.987      2.172   41.89
# [1] "H = Vg^2/ (Vg^2 +Ve^2)"
# [1] "rgxh2012A.csv"
# summary from lme4 is returned
# some computational error has occurred in lmerTest
# Linear mixed model fit by REML ['lmerMod']
# Formula: Mildew_pheno$audpc ~ 1 + (1 | Mildew_pheno$Genotype)
# 
# REML criterion at convergence: 3579.1
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.6991 -0.5924  0.1185  0.5655  2.4689 
# 
# Random effects:
#   Groups                Name        Variance Std.Dev.
# Mildew_pheno$Genotype (Intercept) 111.6    10.56   
# Residual                          103.2    10.16   
# Number of obs: 449, groups:  Mildew_pheno$Genotype, 168
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  57.0375     0.9509   59.98
# [1] "H = Vg^2/ (Vg^2 +Ve^2)"
# [1] "rgxh2013A.csv"
# summary from lme4 is returned
# some computational error has occurred in lmerTest
# Linear mixed model fit by REML ['lmerMod']
# Formula: Mildew_pheno$audpc ~ 1 + (1 | Mildew_pheno$Genotype)
# 
# REML criterion at convergence: 5451.2
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -2.37698 -0.63742  0.05288  0.68795  2.70361 
# 
# Random effects:
#   Groups                Name        Variance Std.Dev.
# Mildew_pheno$Genotype (Intercept)  66.11    8.131  
# Residual                          102.83   10.141  
# Number of obs: 700, groups:  Mildew_pheno$Genotype, 177
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  32.5852     0.7216   45.16
# [1] "H = Vg^2/ (Vg^2 +Ve^2)"
# [1] "rgxh2014A.csv"
# summary from lme4 is returned
# some computational error has occurred in lmerTest
# Linear mixed model fit by REML ['lmerMod']
# Formula: Mildew_pheno$audpc ~ 1 + (1 | Mildew_pheno$Genotype)
# 
# REML criterion at convergence: 6804.4
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.5350 -0.5147 -0.2628  0.5721  3.0468 
# 
# Random effects:
#   Groups                Name        Variance Std.Dev.
# Mildew_pheno$Genotype (Intercept) 737.6    27.16   
# Residual                          657.5    25.64   
# Number of obs: 698, groups:  Mildew_pheno$Genotype, 177
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)   93.978      2.261   41.56
# [1] "H = Vg^2/ (Vg^2 +Ve^2)"
# [1] "rgxh2016A.csv"
# summary from lme4 is returned
# some computational error has occurred in lmerTest
# Linear mixed model fit by REML ['lmerMod']
# Formula: Mildew_pheno$audpc ~ 1 + (1 | Mildew_pheno$Genotype)
# 
# REML criterion at convergence: 5481.5
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -4.8757 -0.5827  0.0982  0.6399  2.8086 
# 
# Random effects:
#   Groups                Name        Variance Std.Dev.
# Mildew_pheno$Genotype (Intercept) 57.34    7.572   
# Residual                          77.49    8.803   
# Number of obs: 734, groups:  Mildew_pheno$Genotype, 123
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  49.9514     0.7562   66.06
# [1] "H = Vg^2/ (Vg^2 +Ve^2)"