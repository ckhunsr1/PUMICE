library(ggplot2)
library(reshape2)
library(dplyr)
library(colortools)

############################################################################################################################################################
predixcan = read.table("/Users/chachritkhun/Desktop/TWAS_files/predicted_exp/ENSG00000152684_DGN_predixcan.txt", header = TRUE, na.string ='.', as.is=TRUE)
colnames(predixcan) = c("observed", "predicted")
predixcan$observed = scale(predixcan$observed)
predixcan$predicted = scale(predixcan$predicted)
predixcan$method = as.factor("PrediXcan")

utmost = read.table("/Users/chachritkhun/Desktop/TWAS_files/predicted_exp/ENSG00000152684_DGN_utmost.txt", header = TRUE, na.string ='.', as.is=TRUE)
colnames(utmost) = c("observed", "predicted")
utmost$observed = scale(utmost$observed)
utmost$predicted = scale(utmost$predicted)
utmost$method = as.factor("UTMOST")

multi = read.table("/Users/chachritkhun/Desktop/TWAS_files/predicted_exp/ENSG00000152684_DGN_multi-omics.txt", header = TRUE, na.string ='.', as.is=TRUE)
colnames(multi) = c("observed", "predicted")
multi$observed = scale(multi$observed)
multi$predicted = scale(multi$predicted)
multi$method = as.factor("PUMICE")

dapg = read.table("/Users/chachritkhun/Desktop/TWAS_files/predicted_exp/ENSG00000152684_DGN_dapg.txt", header = TRUE, na.string ='.', as.is=TRUE)
colnames(dapg) = c("observed", "predicted")
dapg$observed = scale(dapg$observed)
dapg$predicted = scale(dapg$predicted)
dapg$method = as.factor("DAPGW")

blup = read.table("/Users/chachritkhun/Desktop/TWAS_files/predicted_exp/ENSG00000152684_DGN_blup.txt", header = TRUE, na.string ='.', as.is=TRUE)
colnames(blup) = c("observed", "predicted")
blup$observed = scale(blup$observed)
blup$predicted = scale(blup$predicted)
blup$method = as.factor("BLUP")

bslmm = read.table("/Users/chachritkhun/Desktop/TWAS_files/predicted_exp/ENSG00000152684_DGN_bslmm.txt", header = TRUE, na.string ='.', as.is=TRUE)
colnames(bslmm) = c("observed", "predicted")
bslmm$observed = scale(bslmm$observed)
bslmm$predicted = scale(bslmm$predicted)
bslmm$method = as.factor("BSLMM")

dpr_add = read.table("/Users/chachritkhun/Desktop/TWAS_files/predicted_exp/ENSG00000152684_DGN_dpr_add.txt", header = TRUE, na.string ='.', as.is=TRUE)
colnames(dpr_add) = c("observed", "predicted")
dpr_add$observed = scale(dpr_add$observed)
dpr_add$predicted = scale(dpr_add$predicted)
dpr_add$method = as.factor("DPR(VB)")

mcmc_add = read.table("/Users/chachritkhun/Desktop/TWAS_files/predicted_exp/ENSG00000152684_DGN_mcmc_add.txt", header = TRUE, na.string ='.', as.is=TRUE)
colnames(mcmc_add) = c("observed", "predicted")
mcmc_add$observed = scale(mcmc_add$observed)
mcmc_add$predicted = scale(mcmc_add$predicted)
mcmc_add$method = as.factor("DPR(MCMC)")

epixcan = read.table("/Users/chachritkhun/Desktop/TWAS_files/predicted_exp/ENSG00000152684_DGN_epixcan.txt", header = TRUE, na.string ='.', as.is=TRUE)
colnames(epixcan) = c("observed", "predicted")
epixcan$observed = scale(epixcan$observed)
epixcan$predicted = scale(epixcan$predicted)
epixcan$method = as.factor("EpiXcan")

mashr = read.table("/Users/chachritkhun/Desktop/TWAS_files/predicted_exp/ENSG00000152684_DGN_mashr.txt", header = TRUE, na.string ='.', as.is=TRUE)
colnames(mashr) = c("observed", "predicted")
mashr$observed = scale(mashr$observed)
mashr$predicted = scale(mashr$predicted)
mashr$method = as.factor("MASHR")

merge = rbind(multi, bslmm)
merge$method = factor(merge$method, levels = c("PUMICE", "BSLMM"))

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=6, height=6, res=300)
ggplot(merge, aes(predicted, observed, fill = method)) +
  theme_bw() +
  geom_point(shape = 21, alpha = 0.4) +
  geom_smooth(method=lm, color = "black", size = 0.5) +
  scale_fill_manual(values = c("#F18032", "gray"), name = "Method", labels = c(bquote("PUMICE (window = TAD, "*Phi*" = 1/6"*")"), "BSLMM")) +
  #ylim(-2.5, 2.5) +
  #xlim(-2.5, 2.5) +
  theme(legend.position = c(.01, .99), legend.justification = c("left", "top"), legend.box.just = "left", legend.margin = margin(6, 6, 6, 6)) +
  labs(x = "Predicted expression", y = "Observed expression", fill = "Method") +
  #ggtitle("PELO prediction in DGN") +
  #theme(plot.title = element_text(face="bold", size = 15, family = "Helvetica", hjust = 0.5)) +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 14, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 14, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 12, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 12, family = "Helvetica")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title=element_text(face="bold", size=10, family = "Helvetica"), legend.text=element_text(size=8, family = "Helvetica")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) 
  #scale_fill_discrete(name = "Method", labels = c(bquote("PERSIMMON (window = Domain, "*Phi*" = 1/6"*")"), "PrediXcan", "UTMOST"))
  #theme(legend.position = "none") 
dev.off()