# Mildew Correlation graphs 

setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/data/")
library("PerformanceAnalytics")
library(corrplot)
dataexf<-read.csv("emxflocM.csv")
datarxh<-read.csv("rgxhlocM.csv")
my_dataexf <- dataexf[, c(2:7)]
my_datarxh <- datarxh[, c(2:5)]
# colnames(my_datarxh)<- c(2012,2013,2014,2016)
# colnames(my_dataexf)<- c(2011,2012,"2012D",2013,"2013D",2014)
# ^ You have to remove the headers or they over lap with the histogram
colnames(my_datarxh)<- c("","","","")
colnames(my_dataexf)<- c("","","","","","")

# Below is altered source code from library("PerformanceAnalytics") or chart.Correlation
Cor_script<- function (R, histogram = TRUE, method = c("pearson", "kendall", 
                                                "spearman"), ...) 
{
  x = checkData(R, method = "matrix")
  if (missing(method)) 
    method = method[1]
  panel.cor <- function(x, y, digits = 2, prefix = "", use = "pairwise.complete.obs", 
                        method, cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = use, method = method)
    r<-r^2
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if (missing(cex.cor)) 
      cex <- 0.8/strwidth(txt)
    test <- cor.test(x, y, method = method)
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", 
                                                                              "**", "*", ".", " "))
    text(0.5, 0.5, txt, cex = cex * (abs(r) + 0.3)/1.3)
    text(0.7, 0.8, Signif, cex = cex, col = 1)
  }
  f <- function(t) {
    dnorm(t, mean = mean(x), sd = sd.xts(x))
  }
  hist.panel = function(x, ...) {
    par(new = TRUE)
    hist(x, main=NULL, col = "azure", probability = TRUE, axes = FALSE, breaks = "FD")
    lines(density(x, na.rm = TRUE), col = "black", lwd = 2)
    for (i in 1:ncol(R)){
      title(names(x)[i],line = -2, adj = 0, cex=2)
    }
    rug(x)
  }
  if (histogram)
    pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor,
          diag.panel = hist.panel, method = method, ...)
  else pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor,
             method = method, ...)
}

setwd("/Users/helencockerton/Desktop/QTL/Files for git hub/figs/")
pdf("exf_cor.pdf",width=6,height=6,paper='special') 
Cor_script(my_dataexf, histogram=TRUE,col.bg="blue", col="grey18", pch=16, method = c("pearson"))
dev.off()

pdf("rxh_cor.pdf",width=6,height=6,paper='special') 
Cor_script(my_datarxh, histogram=TRUE,col.bg="blue", col="grey18", pch=16, method = "pearson")
dev.off()

# Add year labels to PDF manually 


