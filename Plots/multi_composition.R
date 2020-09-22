library(dplyr)
library(ggplot2)
library(reshape2)

data1 = read.table("/Users/chachritkhun/Desktop/TWAS_files/multi-omics_window_composition.txt", header = TRUE, na.string ='.', as.is=TRUE)
rownames(data1) = data1$tissue
data1 = data1 %>% select(-c("tissue"))

data1$sum = apply(data1, 1, function(x)(sum(x)))
data1$per_1000 = 100*data1$X1000/data1$sum
data1$per_250 = 100*data1$X250/data1$sum
data1$per_Domain = 100*data1$Domain/data1$sum
data1$per_Loop = 100*data1$Loop/data1$sum
data1$per_pcHiC = 100*data1$pcHiC/data1$sum
data1$per_TAD = 100*data1$TAD/data1$sum

data1 = data1 %>% select(per_1000, per_250, per_Domain, per_Loop, per_pcHiC, per_TAD)
data1 = rbind(data1, "gg" = 1:6)
data1 = t(data1)
data1 = as.data.frame(data1)
data1 <- reshape2::melt(data1, "gg")
data1$gg = factor(data1$gg)

tiff("/Users/chachritkhun/Desktop/plot1.tiff", units="in",  width=4, height=4, res=300)
ggplot(data1, aes(x=gg, y=value, color = gg)) + 
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.1), cex=1.2) +
  scale_color_brewer(palette="Dark2") +
  theme_classic() +
  ylab("% Composition") +
  xlab("Window") +
  scale_x_discrete(labels=c("1" = "1000kb", "2" = "250kb", "3" = "Domain", "4" = "Loop", "5" = "pcHiC", "6" = "TAD")) +
  theme(axis.title.x = element_text(face="bold", size = 10, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2.5, face="bold", size = 10, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 8, family = "Helvetica", vjust = 1)) +
  theme(axis.text.y = element_text(face="bold", size = 8, family = "Helvetica")) +
  theme(legend.position = "none")
dev.off()

####################################################################################################################################################
data1 = read.table("/Users/chachritkhun/Desktop/TWAS_files/multi-omics_penalty_composition.txt", header = TRUE, na.string ='.', as.is=TRUE)
rownames(data1) = data1$tissue
data1 = data1 %>% select(-c("tissue"))

data1$sum = apply(data1, 1, function(x)(sum(x)))
data1$per_1 = 100*data1$X1/data1$sum
data1$per_2 = 100*data1$X2/data1$sum
data1$per_3 = 100*data1$X3/data1$sum
data1$per_4 = 100*data1$X4/data1$sum
data1$per_5 = 100*data1$X5/data1$sum
data1$per_6 = 100*data1$X6/data1$sum
data1$per_7 = 100*data1$X7/data1$sum

data1 = data1 %>% select(per_1, per_2, per_3, per_4, per_5, per_6, per_7)
data1 = rbind(data1, "gg" = 1:7)
data1 = t(data1)
data1 = as.data.frame(data1)
data1 <- reshape2::melt(data1, "gg")
data1$gg = factor(data1$gg)

tiff("/Users/chachritkhun/Desktop/plot1.tiff", units="in",  width=4, height=4, res=300)
ggplot(data1, aes(x=gg, y=value, color = gg)) + 
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.1), cex=1.2) +
  scale_color_brewer(palette="Dark2") +
  theme_classic() +
  ylab("% Composition") +
  xlab("Penalty factor for established predictor") +
  scale_x_discrete(labels=c("1" = "0", "2" = "1/6", "3" = "1/3", "4" = "1/2", "5" = "2/3", "6" = "5/6", "7" = "1")) +
  theme(axis.title.x = element_text(face="bold", size = 10, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2.5, face="bold", size = 10, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 8, family = "Helvetica", vjust = 1)) +
  theme(axis.text.y = element_text(face="bold", size = 8, family = "Helvetica")) +
  theme(legend.position = "none")
dev.off()



