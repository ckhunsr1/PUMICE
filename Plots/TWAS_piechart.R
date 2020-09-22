library(dplyr)
library(ggplot2)
library(reshape2)

#data1 = read.table("/Users/chachritkhun/Desktop/TWAS_files/snp_count/Whole_Blood_multi-omics_snps_count.txt", header = TRUE, na.string ='.', as.is=TRUE)
data1 = read.table("/Users/chachritkhun/Desktop/TWAS_files/snp_count/multi-omics_snps_count.txt", header = TRUE, na.string ='.', as.is=TRUE)
df1 = as.data.frame(cbind("method" = "PUMICE", "num" = data1$n))
df1$num = as.numeric(as.character(df1$num))

#data2 = read.table("/Users/chachritkhun/Desktop/TWAS_files/snp_count/Whole_Blood_predixcan_snps_count.txt", header = TRUE, na.string ='.', as.is=TRUE)
data2 = read.table("/Users/chachritkhun/Desktop/TWAS_files/snp_count/predixcan_snps_count.txt", header = TRUE, na.string ='.', as.is=TRUE)
df2 = as.data.frame(cbind("method" = "PrediXcan", "num" = data2$n))
df2$num = as.numeric(as.character(df2$num))

#data3 = read.table("/Users/chachritkhun/Desktop/TWAS_files/snp_count/Whole_Blood_utmost_snps_count.txt", header = TRUE, na.string ='.', as.is=TRUE)
data3 = read.table("/Users/chachritkhun/Desktop/TWAS_files/snp_count/utmost_snps_count.txt", header = TRUE, na.string ='.', as.is=TRUE)
df3 = as.data.frame(cbind("method" = "UTMOST", "num" = data3$n))
df3$num = as.numeric(as.character(df3$num))

#data4 = read.table("/Users/chachritkhun/Desktop/TWAS_files/snp_count/Whole_Blood_dapg_snps_count.txt", header = TRUE, na.string ='.', as.is=TRUE)
data4 = read.table("/Users/chachritkhun/Desktop/TWAS_files/snp_count/dapg_snps_count.txt", header = TRUE, na.string ='.', as.is=TRUE)
df4 = as.data.frame(cbind("method" = "DAPGW", "num" = data4$n))
df4$num = as.numeric(as.character(df4$num))

#data5 = read.table("/Users/chachritkhun/Desktop/TWAS_files/snp_count/Whole_Blood_mashr_snps_count.txt", header = TRUE, na.string ='.', as.is=TRUE)
data5 = read.table("/Users/chachritkhun/Desktop/TWAS_files/snp_count/mashr_snps_count.txt", header = TRUE, na.string ='.', as.is=TRUE)
df5 = as.data.frame(cbind("method" = "MASHR", "num" = data5$n))
df5$num = as.numeric(as.character(df5$num))

#data6 = read.table("/Users/chachritkhun/Desktop/TWAS_files/snp_count/Whole_Blood_epixcan_snps_count.txt", header = TRUE, na.string ='.', as.is=TRUE)
#df6 = as.data.frame(cbind("method" = "EpiXcan", "num" = data6$n))
#df6$num = as.numeric(as.character(df6$num))

df = rbind(df2, df4, df1, df5, df3)
medians = cbind(median(df2$num), median(df4$num), median(df1$num), median(df5$num), median(df3$num))
means = cbind(mean(df2$num), mean(df4$num), mean(df1$num), mean(df5$num), mean(df3$num))
  
#df = rbind(df2, df4, df6, df1, df5, df3)
#medians = cbind(median(df2$num), median(df4$num), median(df6$num), median(df1$num), median(df5$num), median(df3$num))
#means = cbind(mean(df2$num), mean(df4$num), mean(df6$num), mean(df1$num), mean(df5$num), mean(df3$num))

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=10, height=5, res=200)
ggplot(df, aes(x = num, fill = method)) +
  geom_density(alpha=0.3, bw = 10) +
  geom_vline(xintercept = medians[1,1], linetype="dashed", color = "#59268F", size=1) +
  geom_vline(xintercept = medians[1,2], linetype="dashed", color = "#24A860", size=1) +
  #geom_vline(xintercept = medians[1,3], linetype="dashed", color = "gray", size=1) +
  geom_vline(xintercept = medians[1,3], linetype="dashed", color = "#F18032", size=1) +
  geom_vline(xintercept = medians[1,4], linetype="dashed", color = "#FEF236", size=1) +
  geom_vline(xintercept = medians[1,5], linetype="dashed", color = "#1265B2", size=1) +
  theme_classic() +
  scale_x_continuous(limits = c(0,100), breaks = c(0, medians[1,5], medians[1,4], medians[1,1], medians[1,2], 50, medians[1,3], 100)) +
  xlab("Number of non-zero SNP weight") +
  ylab("Density") +
  scale_fill_manual(values=c("#59268F", "#24A860", "#F18032", "#FEF236", "#1265B2")) +
  #scale_fill_manual(values=c("#59268F", "#24A860", "gray", "#F18032", "#FEF236", "#1265B2")) +
  labs(fill = "Method") +
  theme(legend.text = element_text(colour="black", size=8, family = "Helvetica")) +
  theme(legend.title = element_text(colour="black", size=10, face="bold", family = "Helvetica")) +
  theme(legend.position=c(0.98,0.5), legend.justification=c(0.98,0.5)) +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 14, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 10, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 10, family = "Helvetica")) 
dev.off()

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=10, height=5, res=200)
ggplot(df, aes(x = num, fill = method)) +
  geom_density(alpha=0.3, bw = 10) +
  geom_vline(xintercept = medians[1,1], linetype="dashed", color = "#F18032", size=1) +
  geom_vline(xintercept = medians[1,2], linetype="dashed", color = "#59268F", size=1) +
  geom_vline(xintercept = medians[1,3], linetype="dashed", color = "#1265B2", size=1) +
  geom_vline(xintercept = medians[1,4], linetype="dashed", color = "#24A860", size=1) +
  geom_vline(xintercept = medians[1,5], linetype="dashed", color = "#FEF236", size=1) +
  geom_vline(xintercept = medians[1,6], linetype="dashed", color = "gray", size=1) +
  theme_classic() +
  scale_x_continuous(limits = c(0,100), breaks = c(0, medians[1,5], medians[1,4], medians[1,1], medians[1,6], medians[1,2], 50, medians[1,3], 100)) +
  xlab("Number of non-zero SNP weight") +
  ylab("Density") +
  scale_fill_manual(values=c("#F18032", "#59268F", "#1265B2", "#24A860", "#FEF236", "gray")) +
  labs(fill = "Method") +
  theme(legend.text = element_text(colour="black", size=8, family = "Helvetica")) +
  theme(legend.title = element_text(colour="black", size=10, face="bold", family = "Helvetica")) +
  theme(legend.position=c(0.98,0.5), legend.justification=c(0.98,0.5)) +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 14, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 10, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 10, family = "Helvetica")) 
dev.off()

##Barplot for func type##
tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=5, height=5, res=200)
composition = as.data.frame(table(data1$V3))
ggplot(composition, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values=c("#3399FF", "#3399FF", "#FF9900", "#FF9900", "#FF9900", "#FF9900"))
dev.off()

##Barplot for penalty factor##
tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=5, height=5, res=200)
composition = as.data.frame(table(data1$V7))
ggplot(composition, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()
dev.off()





