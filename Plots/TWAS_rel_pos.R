library(dplyr)
library(ggplot2)
library(reshape2)

##Containing 98.25% of non-zero weight SNPs from Multi-omics##
#data1 = read.table("/Users/chachritkhun/Desktop/TWAS_files/rel_pos/Whole_Blood_multi-omics_rel_pos.txt", header = TRUE, na.string ='.', as.is=TRUE)
data1 = read.table("/Users/chachritkhun/Desktop/TWAS_files/rel_pos/multi-omics_rel_pos.txt", header = TRUE, na.string ='.', as.is=TRUE)
data1$rel_bin = as.numeric(as.character(data1$rel_bin))

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=5, height=5, res=200)
ggplot(data1) +
  geom_histogram(aes(x = rel_bin, y = ..density..), binwidth = 1, fill = "#F18032", color = "black") +
  geom_vline(xintercept=c(50,100), colour="red", size = 0.8, alpha = 0.6) +
  theme_classic() +
  xlab("Position relative to gene (kb)") +
  ylab("Density") +
  ggtitle("PUMICE") +
  ylim(0, 0.08) +
  scale_x_continuous(labels=  c("-1000kb", "TSS", "TES", "+1000kb")) +
  theme(plot.title = element_text(face="bold", size = 10, family = "Helvetica", hjust = 0.5)) +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 14, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 10, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 10, family = "Helvetica"))
dev.off()

##Containing100% of non-zero weight SNPs from PrediXcan##
data1 = read.table("/Users/chachritkhun/Desktop/TWAS_files/rel_pos/Whole_Blood_predixcan_rel_pos.txt", header = TRUE, na.string ='.', as.is=TRUE)
#data1 = read.table("/Users/chachritkhun/Desktop/TWAS_files/rel_pos/predixcan_rel_pos.txt", header = TRUE, na.string ='.', as.is=TRUE)
data1$rel_bin = as.numeric(as.character(data1$rel_bin))

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=5, height=5, res=200)
ggplot(data1) +
  geom_histogram(aes(x = rel_bin, y = ..density..), binwidth = 1, fill = "#59268F", color = "black") +
  geom_vline(xintercept=c(50,100), colour="red", size = 0.8, alpha = 0.6) +
  theme_classic() +
  xlab("Position relative to gene (kb)") +
  ylab("Density") +
  ggtitle("PrediXcan") +
  ylim(0, 0.08) +
  scale_x_continuous(labels=  c("-1000kb", "TSS", "TES", "+1000kb")) +
  theme(plot.title = element_text(face="bold", size = 10, family = "Helvetica", hjust = 0.5)) +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 14, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 10, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 10, family = "Helvetica"))
dev.off()


##Containing100% of non-zero weight SNPs from UTMOST##
data1 = read.table("/Users/chachritkhun/Desktop/TWAS_files/rel_pos/Whole_Blood_utmost_rel_pos.txt", header = TRUE, na.string ='.', as.is=TRUE)
#data1 = read.table("/Users/chachritkhun/Desktop/TWAS_files/rel_pos/utmost_rel_pos.txt", header = TRUE, na.string ='.', as.is=TRUE)
data1$rel_bin = as.numeric(as.character(data1$rel_bin))

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=5, height=5, res=200)
ggplot(data1) +
  geom_histogram(aes(x = rel_bin, y = ..density..), binwidth = 1, fill = "#1265B2", color = "black") +
  geom_vline(xintercept=c(50,100), colour="red", size = 0.8, alpha = 0.6) +
  theme_classic() +
  xlab("Position relative to gene (kb)") +
  ylab("Density") +
  ggtitle("UTMOST") +
  ylim(0, 0.08) +
  scale_x_continuous(labels=  c("-1000kb", "TSS", "TES", "+1000kb")) +
  theme(plot.title = element_text(face="bold", size = 10, family = "Helvetica", hjust = 0.5)) +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 14, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 10, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 10, family = "Helvetica"))
dev.off()

##Containing100% of non-zero weight SNPs from DAPGW##
data1 = read.table("/Users/chachritkhun/Desktop/TWAS_files/rel_pos/Whole_Blood_dapg_rel_pos.txt", header = TRUE, na.string ='.', as.is=TRUE)
#data1 = read.table("/Users/chachritkhun/Desktop/TWAS_files/rel_pos/dapg_rel_pos.txt", header = TRUE, na.string ='.', as.is=TRUE)
data1$rel_bin = as.numeric(as.character(data1$rel_bin))

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=5, height=5, res=200)
ggplot(data1) +
  geom_histogram(aes(x = rel_bin, y = ..density..), binwidth = 1, fill = "#24A860", color = "black") +
  geom_vline(xintercept=c(50,100), colour="red", size = 0.8, alpha = 0.6) +
  theme_classic() +
  xlab("Position relative to gene (kb)") +
  ylab("Density") +
  ggtitle("DAPGW") +
  ylim(0, 0.08) +
  scale_x_continuous(labels=  c("-1000kb", "TSS", "TES", "+1000kb")) +
  theme(plot.title = element_text(face="bold", size = 10, family = "Helvetica", hjust = 0.5)) +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 14, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 10, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 10, family = "Helvetica"))
dev.off()

##Containing100% of non-zero weight SNPs from MASHR##
data1 = read.table("/Users/chachritkhun/Desktop/TWAS_files/rel_pos/Whole_Blood_mashr_rel_pos.txt", header = TRUE, na.string ='.', as.is=TRUE)
#data1 = read.table("/Users/chachritkhun/Desktop/TWAS_files/rel_pos/mashr_rel_pos.txt", header = TRUE, na.string ='.', as.is=TRUE)
data1$rel_bin = as.numeric(as.character(data1$rel_bin))

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=5, height=5, res=200)
ggplot(data1) +
  geom_histogram(aes(x = rel_bin, y = ..density..), binwidth = 1, fill = "#FEF236", color = "black") +
  geom_vline(xintercept=c(50,100), colour="red", size = 0.8, alpha = 0.6) +
  theme_classic() +
  xlab("Position relative to gene (kb)") +
  ylab("Density") +
  ggtitle("MASHR") +
  ylim(0, 0.08) +
  scale_x_continuous(labels=  c("-1000kb", "TSS", "TES", "+1000kb")) +
  theme(plot.title = element_text(face="bold", size = 10, family = "Helvetica", hjust = 0.5)) +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 14, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 10, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 10, family = "Helvetica"))
dev.off()


##Containing 100% of non-zero weight SNPs from EpiXcan##
data1 = read.table("/Users/chachritkhun/Desktop/TWAS_files/rel_pos/Whole_Blood_epixcan_rel_pos.txt", header = TRUE, na.string ='.', as.is=TRUE)
data1$rel_bin = as.numeric(as.character(data1$rel_bin))

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=5, height=5, res=200)
ggplot(data1) +
  geom_histogram(aes(x = rel_bin, y = ..density..), binwidth = 1, fill = "gray", color = "black") +
  geom_vline(xintercept=c(50,100), colour="red", size = 0.8, alpha = 0.6) +
  theme_classic() +
  xlab("Position relative to gene (kb)") +
  ylab("Density") +
  ggtitle("EpiXcan") +
  ylim(0, 0.08) +
  scale_x_continuous(labels=  c("-1000kb", "TSS", "TES", "+1000kb")) +
  theme(plot.title = element_text(face="bold", size = 10, family = "Helvetica", hjust = 0.5)) +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 14, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 10, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 10, family = "Helvetica"))
dev.off()
