library(qqman)
library(data.table)
library(lattice)
library(dplyr)

qqplot <- function(pvector, col = "black", add = F, ylim = NULL, shape = 21){  
  expectedP <- -log10(ppoints(length(pvector)))
  observedP <- -log10(sort(pvector, decreasing = F))
  if (add == F) {
    plot(x = expectedP, y = observedP, col = col,
         xlab = expression(Expected ~ ~-log[10](italic(p))),
         ylab = expression(Observed ~ ~-log[10](italic(p))),
         pch = shape, cex = 0.7, ylim=c(-10,260))
    abline(0, 1, col = "black")
  }else{
    points(x = expectedP, y = observedP, col = col,
           xlab = expression(Expected ~ ~-log[10](italic(p))),
           ylab = expression(Observed ~ ~-log[10](italic(p))),
           pch = shape, cex = 0.7) 
  }
}

predixcan = as.data.frame(fread("/Users/chachritkhun/Desktop/TWAS_files/QQ/predixcan_QQ_analysis2.txt"))
predixcan = predixcan %>% filter(V1 > 0)

dapg = as.data.frame(fread("/Users/chachritkhun/Desktop/TWAS_files/QQ/dapg_QQ_analysis2.txt"))
dapg = dapg %>% filter(V1 > 0)

multi = as.data.frame(fread("/Users/chachritkhun/Desktop/TWAS_files/QQ/multi-omics_QQ_analysis2.txt"))
multi = multi %>% filter(V1 > 0)

utmost = as.data.frame(fread("/Users/chachritkhun/Desktop/TWAS_files/QQ/utmost_QQ_analysis2.txt"))
utmost = utmost %>% filter(V1 > 0)

mashr = as.data.frame(fread("/Users/chachritkhun/Desktop/TWAS_files/QQ/mashr_QQ_analysis2.txt"))
mashr = mashr %>% filter(V1 > 0)

predixcan_mashr_cauchy = as.data.frame(fread("/Users/chachritkhun/Desktop/TWAS_files/QQ/predixcan_mashr_cauchy_QQ_analysis2.txt"))
predixcan_mashr_cauchy = predixcan_mashr_cauchy %>% filter(V1 > 0)

predixcan_utmost_cauchy = as.data.frame(fread("/Users/chachritkhun/Desktop/TWAS_files/QQ/predixcan_utmost_cauchy_QQ_analysis2.txt"))
predixcan_utmost_cauchy = predixcan_utmost_cauchy %>% filter(V1 > 0)

dapg_mashr_cauchy = as.data.frame(fread("/Users/chachritkhun/Desktop/TWAS_files/QQ/dapg_mashr_cauchy_QQ_analysis2.txt"))
dapg_mashr_cauchy = dapg_mashr_cauchy %>% filter(V1 > 0)

dapg_utmost_cauchy = as.data.frame(fread("/Users/chachritkhun/Desktop/TWAS_files/QQ/dapg_utmost_cauchy_QQ_analysis2.txt"))
dapg_utmost_cauchy = dapg_utmost_cauchy %>% filter(V1 > 0)

combined_cauchy = as.data.frame(fread("/Users/chachritkhun/Desktop/TWAS_files/QQ/combined_cauchy_QQ_analysis2.txt"))
combined_cauchy = combined_cauchy %>% filter(V1 > 0)

combined_mashr_cauchy = as.data.frame(fread("/Users/chachritkhun/Desktop/TWAS_files/QQ/combined_mashr_cauchy_QQ_analysis2.txt"))
combined_mashr_cauchy = combined_mashr_cauchy %>% filter(V1 > 0)


tiff("/Users/chachritkhun/Desktop/plot1.tiff", units="in",  width=5, height=5, res=300)
qqplot(predixcan_utmost_cauchy$V1, col = "#59268F", shape = 24, add = F)
qqplot(dapg_utmost_cauchy$V1, col = "#24A860", shape = 24, add = T)
qqplot(combined_cauchy$V1, col = "#F18032", shape = 24, add = T)
qqplot(combined_mashr_cauchy$V1, col = "#F18032", shape = 23, add = T)
dev.off()

tiff("/Users/chachritkhun/Desktop/plot2.tiff", units="in",  width=5, height=5, res=300)
qqplot(predixcan_mashr_cauchy$V1, col = "#59268F", shape = 23, add = F)
qqplot(dapg_mashr_cauchy$V1, col = "#24A860", shape = 23, add = T)
qqplot(combined_mashr_cauchy$V1, col = "#F18032", shape = 23, add = T)
dev.off()

tiff("/Users/chachritkhun/Desktop/plot3.tiff", units="in",  width=5, height=5, res=300)
qqplot(mashr$V1, col = "#FEF236", shape = 21, add = F)
qqplot(dapg$V1, col = "#24A860", shape = 21, add = T)
qqplot(utmost$V1, col = "#1265B2", shape = 21, add = T)
qqplot(predixcan$V1, col = "#59268F", shape = 21, add = T)
qqplot(multi$V1, col = "#F18032", shape = 21, add = T)
qqplot(combined_mashr_cauchy$V1, col = "#F18032", shape = 23, add = T)
dev.off()




