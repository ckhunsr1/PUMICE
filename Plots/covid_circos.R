library(dplyr)
library(data.table)
library(ggplot2)
library(CMplot)

df = read.table("/Users/chachritkhun/Desktop/TWAS_files/covid_pval.txt", header = TRUE)
df = df %>% select(SNP, Chromosome, Position, utmost, predixcan, mashr, dapg, multi.omics, combined, combined_mashr)
df = df %>% select(SNP, Chromosome, Position,combined, combined_mashr)

CMplot(df,type="p",plot.type="c",chr.labels=paste("Chr",c(1:22),sep=""),r=0.4,cir.legend=TRUE,
       outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",file="tiff",  ylim = c(0,10),
       memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10)


df = read.table("/Users/chachritkhun/Desktop/TWAS_files/covid_pval.txt", header = TRUE)
df = df %>% select(SNP, Chromosome, Position, combined)
colnames(df) = c("SNP", "CHR", "BP", "P")
gwasResults = df

don <- gwasResults %>% group_by(CHR) %>% summarise(chr_len=max(BP)) %>% mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>% select(-chr_len) %>%
       left_join(gwasResults, ., by=c("CHR"="CHR")) %>% arrange(CHR, BP) %>% mutate( BPcum=BP+tot)
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=10, height=2, res=300)
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  geom_point( aes(color=as.factor(CHR)), alpha=1, size=0.5) +
  #geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2, nudge_y = c(7, 6, 5, 8, 8)) +
  geom_hline(yintercept = -log10(0.05/21000), color = "red", size = 0.3) + 
  scale_color_manual(values = rep(c("#14234A", "#DC5E26"), 22 )) +
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  ylim(0,10) + 
  theme_classic() +
  theme(legend.position="none", panel.border = element_blank(), panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +
  labs(x = "Chromosome", y = "-log(P)") +
  theme(axis.title.x = element_text(vjust=1.5, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=1.5, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.text.x = element_text(angle = 30, size = 10, family = "Helvetica", hjust = 1)) +
  theme(axis.text.y = element_text(size = 10, family = "Helvetica"))
dev.off()

#############################################

df = as.data.frame(fread("/Users/chachritkhun/Desktop/TWAS_files/ss_clean_covid.txt"))
colnames(df) = c("SNP", "CHR", "BP", "P")
gwasResults = df %>% filter(-log10(P)>1)

don <- gwasResults %>% group_by(CHR) %>% summarise(chr_len=max(BP)) %>% mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>% select(-chr_len) %>%
  left_join(gwasResults, ., by=c("CHR"="CHR")) %>% arrange(CHR, BP) %>% mutate( BPcum=BP+tot)
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=10, height=2, res=300)
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  geom_point( aes(color=as.factor(CHR)), alpha=1, size=0.5) +
  geom_hline(yintercept = -log10(0.05/1000000), color = "red", size = 0.3) + 
  scale_color_manual(values = rep(c("#14234A", "#DC5E26"), 22 )) +
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  theme_classic() +
  theme(legend.position="none", panel.border = element_blank(), panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) +
  labs(x = "Chromosome", y = "-log(P)") +
  theme(axis.title.x = element_text(vjust=1.5, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=1.5, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.text.x = element_text(angle = 30, size = 10, family = "Helvetica", hjust = 1)) +
  theme(axis.text.y = element_text(size = 10, family = "Helvetica"))
dev.off()



