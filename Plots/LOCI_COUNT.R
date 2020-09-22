library(ggplot2)
library(ggpattern)
library(reshape2)

##Loci count##
options(scipen=10000)

##Create bar plot for number of models##
method = c("PrediXcan", "DAPGW", "PUMICE", "UTMOST", "MASHR")
value = c(262426, 517973, 383465, 341283, 509642)

data <- data.frame(method, value)
data$method = factor(data$method, levels = method)

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=4, height=4, res=300)
ggplot(data=data, aes(x=method, y=value, fill = method)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("PrediXcan" = "#59268F", "DAPGW" = "#24A860", "MASHR" = "#FEF236", "UTMOST" = "#1265B2", "PUMICE" = "#F18032")) +
  theme_bw() +
  xlab("") +
  ylab("Number of models") +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.text.x = element_text(size = 8, family = "Helvetica", angle = 30, vjust = 0.6)) +
  theme(axis.text.y = element_text(size = 8, family = "Helvetica"))
dev.off()

###################################################################################################################
##Create bar plot for GTA counts##
method = c("PrediXcan", "DAPGW", "PUMICE", "UTMOST", "MASHR", 
           "PrediXcan+(UTMOST)", "DAPGW+(UTMOST)", "PUMICE+(UTMOST)",
           "PrediXcan+(MASHR)", "DAPGW+(MASHR)", "PUMICE+(MASHR)")
value = c(4.4001, 4.3955, 5.5335, 7.1523, 6.0509, 8.2151, 8.9257, 8.7462, 7.4562, 6.5010, 8.3885)

data <- data.frame(method, value)
data$method = factor(data$method, levels = method)

col1 = c("#59268F", "#24A860", "#F18032", "#1265B2", "#FEF236", "#1265B2", "#1265B2", "#1265B2", "#FEF236", "#FEF236", "#FEF236")
col2 = c("#59268F", "#24A860", "#F18032", "#1265B2", "#FEF236", "#59268F", "#24A860", "#F18032", "#59268F", "#24A860", "#F18032")
tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=5, height=4, res=300)
ggplot(data, aes(method, value)) +
  geom_col_pattern(
    aes(pattern = method, pattern_fill = method, pattern_colour = method),
    pattern = c('stripe'),
    pattern_alpha = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1),
    fill = col1,
    colour = 'black', pattern_density = 0.5, pattern_angle = 30) +
  scale_pattern_color_manual(values = col1) +
  scale_pattern_fill_manual(values = col2) +
    theme(legend.position = 'none') +
  #coord_fixed(ratio = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  theme(legend.position = "none") +
  scale_y_continuous(breaks=c(0, 2.5, 5, 7.5), labels=c("0", "25000", "50000", "75000")) +
  labs(x = "", y = "Number of gene-trait associations") +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.text.x = element_text(size = 8, family = "Helvetica")) +
  theme(axis.text.y = element_text(size = 8, family = "Helvetica"))
dev.off()

##Create bar plot for GTA counts##
method = c("PrediXcan", "DAPGW", "PUMICE", "UTMOST", "MASHR","PUMICE+(MASHR)")
value = c(4.4001, 4.3955, 5.5335, 7.1523, 6.0509, 8.3885)

data <- data.frame(method, value)
data$method = factor(data$method, levels = method)

col1 = c("#59268F", "#24A860", "#F18032", "#1265B2", "#FEF236", "#FEF236")
col2 = c("#59268F", "#24A860", "#F18032", "#1265B2", "#FEF236", "#F18032")
tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=4, height=4, res=300)
ggplot(data, aes(method, value)) +
  geom_col_pattern(
    aes(pattern = method, pattern_fill = method, pattern_colour = method),
    pattern = c('stripe'),
    pattern_alpha = c(0, 0, 0, 0, 0, 1),
    fill = col1,
    colour = 'black', pattern_density = 0.5, pattern_angle = 30) +
  scale_pattern_color_manual(values = col1) +
  scale_pattern_fill_manual(values = col2) +
  theme(legend.position = 'none') +
  #coord_fixed(ratio = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  theme(legend.position = "none") +
  scale_y_continuous(breaks=c(0, 2.5, 5, 7.5), labels=c("0", "25000", "50000", "75000")) +
  labs(x = "", y = "Number of gene-trait associations") +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.text.x = element_text(size = 8, family = "Helvetica")) +
  theme(axis.text.y = element_text(size = 8, family = "Helvetica"))
dev.off()

###################################################################################################################
##Create bar plot for unique gene counts##
method = c("PrediXcan", "DAPGW", "PUMICE", "UTMOST", "MASHR", 
           "PrediXcan+(UTMOST)", "DAPGW+(UTMOST)", "PUMICE+(UTMOST)",
           "PrediXcan+(MASHR)", "DAPGW+(MASHR)", "PUMICE+(MASHR)")
value = c(11.511, 16.207, 15.935, 8.017, 16.695, 13.587, 18.550, 17.978, 20.171, 18.295, 22.968)

data <- data.frame(method, value)
data$method = factor(data$method, levels = method)

col1 = c("#59268F", "#24A860", "#F18032", "#1265B2", "#FEF236", "#1265B2", "#1265B2", "#1265B2", "#FEF236", "#FEF236", "#FEF236")
col2 = c("#59268F", "#24A860", "#F18032", "#1265B2", "#FEF236", "#59268F", "#24A860", "#F18032", "#59268F", "#24A860", "#F18032")
tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=5, height=4, res=300)
ggplot(data, aes(method, value)) +
  geom_col_pattern(
    aes(pattern = method, pattern_fill = method, pattern_colour = method),
    pattern = c('stripe'),
    pattern_alpha = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1),
    fill = col1,
    colour = 'black', pattern_density = 0.5, pattern_angle = 30) +
  scale_pattern_color_manual(values = col1) +
  scale_pattern_fill_manual(values = col2) +
  theme(legend.position = 'none') +
  #coord_fixed(ratio = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  theme(legend.position = "none") +
  scale_y_continuous(breaks=c(0, 5, 10, 15, 20), labels=c("0", "5000", "10000", "15000", "20000")) +
  labs(x = "", y = "Unique gene count") +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.text.x = element_text(size = 8, family = "Helvetica")) +
  theme(axis.text.y = element_text(size = 8, family = "Helvetica"))
dev.off()

method = c("PrediXcan", "DAPGW", "PUMICE", "UTMOST", "MASHR", "PUMICE+(MASHR)")
value = c(11.511, 16.207, 15.935, 8.017, 16.695, 22.968)

data <- data.frame(method, value)
data$method = factor(data$method, levels = method)

col1 = c("#59268F", "#24A860", "#F18032", "#1265B2", "#FEF236", "#FEF236")
col2 = c("#59268F", "#24A860", "#F18032", "#1265B2", "#FEF236", "#F18032")
tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=4, height=4, res=300)
ggplot(data, aes(method, value)) +
  geom_col_pattern(
    aes(pattern = method, pattern_fill = method, pattern_colour = method),
    pattern = c('stripe'),
    pattern_alpha = c(0, 0, 0, 0, 0, 1),
    fill = col1,
    colour = 'black', pattern_density = 0.5, pattern_angle = 30) +
  scale_pattern_color_manual(values = col1) +
  scale_pattern_fill_manual(values = col2) +
  theme(legend.position = 'none') +
  #coord_fixed(ratio = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  theme(legend.position = "none") +
  scale_y_continuous(breaks=c(0, 5, 10, 15, 20), labels=c("0", "5000", "10000", "15000", "20000")) +
  labs(x = "", y = "Unique gene count") +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.text.x = element_text(size = 8, family = "Helvetica")) +
  theme(axis.text.y = element_text(size = 8, family = "Helvetica"))
dev.off()

###################################################################################################################
##Create gene-trait pair plot (MAGMA)##
method = c("PrediXcan", "DAPGW", "PUMICE", "UTMOST", "MASHR", 
           "PrediXcan+(UTMOST)", "DAPGW+(UTMOST)", "PUMICE+(UTMOST)",
           "PrediXcan+(MASHR)", "DAPGW+(MASHR)", "PUMICE+(MASHR)")
value = c(6.937, 10.566, 9.982, 4.549, 10.839, 8.393, 12.309, 11.544, 13.622, 12.104, 15.751)

data <- data.frame(method, value)
data$method = factor(data$method, levels = method)

col1 = c("#59268F", "#24A860", "#F18032", "#1265B2", "#FEF236", "#1265B2", "#1265B2", "#1265B2", "#FEF236", "#FEF236", "#FEF236")
col2 = c("#59268F", "#24A860", "#F18032", "#1265B2", "#FEF236", "#59268F", "#24A860", "#F18032", "#59268F", "#24A860", "#F18032")
tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=5, height=4, res=300)
ggplot(data, aes(method, value)) +
  geom_col_pattern(
    aes(pattern = method, pattern_fill = method, pattern_colour = method),
    pattern = c('stripe'),
    pattern_alpha = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1),
    fill = col1,
    colour = 'black', pattern_density = 0.5, pattern_angle = 30) +
  scale_pattern_color_manual(values = col1) +
  scale_pattern_fill_manual(values = col2) +
  theme(legend.position = 'none') +
  #coord_fixed(ratio = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  theme(legend.position = "none") +
  scale_y_continuous(breaks=c(0, 5, 10, 15), labels=c("0", "5000", "10000", "15000")) +
  labs(x = "", y = "Novel gene count (MAGMA)") +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.text.x = element_text(size = 8, family = "Helvetica")) +
  theme(axis.text.y = element_text(size = 8, family = "Helvetica"))
dev.off()

method = c("PrediXcan", "DAPGW", "PUMICE", "UTMOST", "MASHR", "PUMICE+(MASHR)")
value = c(6.937, 10.566, 9.982, 4.549, 10.839, 15.751)

data <- data.frame(method, value)
data$method = factor(data$method, levels = method)

col1 = c("#59268F", "#24A860", "#F18032", "#1265B2", "#FEF236", "#FEF236")
col2 = c("#59268F", "#24A860", "#F18032", "#1265B2", "#FEF236", "#F18032")
tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=4, height=4, res=300)
ggplot(data, aes(method, value)) +
  geom_col_pattern(
    aes(pattern = method, pattern_fill = method, pattern_colour = method),
    pattern = c('stripe'),
    pattern_alpha = c(0, 0, 0, 0, 0, 1),
    fill = col1,
    colour = 'black', pattern_density = 0.5, pattern_angle = 30) +
  scale_pattern_color_manual(values = col1) +
  scale_pattern_fill_manual(values = col2) +
  theme(legend.position = 'none') +
  #coord_fixed(ratio = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  theme(legend.position = "none") +
  scale_y_continuous(breaks=c(0, 5, 10, 15), labels=c("0", "5000", "10000", "15000")) +
  labs(x = "", y = "Novel gene count (MAGMA)") +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.text.x = element_text(size = 8, family = "Helvetica")) +
  theme(axis.text.y = element_text(size = 8, family = "Helvetica"))
dev.off()


###################################################################################################################
##Create gene-trait pair plot (GWAS)##
method = c("PrediXcan", "DAPGW", "PUMICE", "UTMOST", "MASHR", 
           "PrediXcan+(UTMOST)", "DAPGW+(UTMOST)", "PUMICE+(UTMOST)",
           "PrediXcan+(MASHR)", "DAPGW+(MASHR)", "PUMICE+(MASHR)")
value = c(2.711, 2.279, 3.744, 1.834, 2.397, 3.204, 3.179, 4.166, 3.654, 2.693, 4.511)

data <- data.frame(method, value)
data$method = factor(data$method, levels = method)

col1 = c("#59268F", "#24A860", "#F18032", "#1265B2", "#FEF236", "#1265B2", "#1265B2", "#1265B2", "#FEF236", "#FEF236", "#FEF236")
col2 = c("#59268F", "#24A860", "#F18032", "#1265B2", "#FEF236", "#59268F", "#24A860", "#F18032", "#59268F", "#24A860", "#F18032")
tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=5, height=4, res=300)
ggplot(data, aes(method, value)) +
  geom_col_pattern(
    aes(pattern = method, pattern_fill = method, pattern_colour = method),
    pattern = c('stripe'),
    pattern_alpha = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1),
    fill = col1,
    colour = 'black', pattern_density = 0.5, pattern_angle = 30) +
  scale_pattern_color_manual(values = col1) +
  scale_pattern_fill_manual(values = col2) +
  theme(legend.position = 'none') +
  #coord_fixed(ratio = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  theme(legend.position = "none") +
  scale_y_continuous(breaks=c(0, 1, 2, 3, 4), labels=c("0", "1000", "2000", "3000" , "4000")) +
  labs(x = "", y = "Novel gene count (GWAS)") +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.text.x = element_text(size = 8, family = "Helvetica")) +
  theme(axis.text.y = element_text(size = 8, family = "Helvetica"))
dev.off()

method = c("PrediXcan", "DAPGW", "PUMICE", "UTMOST", "MASHR","PUMICE+(MASHR)")
value = c(2.711, 2.279, 3.744, 1.834, 2.397,4.511)

data <- data.frame(method, value)
data$method = factor(data$method, levels = method)

col1 = c("#59268F", "#24A860", "#F18032", "#1265B2", "#FEF236", "#FEF236")
col2 = c("#59268F", "#24A860", "#F18032", "#1265B2", "#FEF236", "#F18032")
tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=4, height=4, res=300)
ggplot(data, aes(method, value)) +
  geom_col_pattern(
    aes(pattern = method, pattern_fill = method, pattern_colour = method),
    pattern = c('stripe'),
    pattern_alpha = c(0, 0, 0, 0, 0, 1),
    fill = col1,
    colour = 'black', pattern_density = 0.5, pattern_angle = 30) +
  scale_pattern_color_manual(values = col1) +
  scale_pattern_fill_manual(values = col2) +
  theme(legend.position = 'none') +
  #coord_fixed(ratio = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  theme(legend.position = "none") +
  scale_y_continuous(breaks=c(0, 1, 2, 3, 4), labels=c("0", "1000", "2000", "3000" , "4000")) +
  labs(x = "", y = "Novel gene count (GWAS)") +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.text.x = element_text(size = 8, family = "Helvetica")) +
  theme(axis.text.y = element_text(size = 8, family = "Helvetica"))
dev.off()

############################################################################################################
##COVID result##
##Create bar plot for GTA counts##
method = c("PrediXcan", "DAPGW", "MASHR", "UTMOST", "PERSIMMON", "Cauchy")
value = c(2, 3, 2, 1, 187, 177)

data <- data.frame(method, value)
data$method = factor(data$method, levels = method)

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=4, height=4, res=300)
ggplot(data=data, aes(x=method, y=value, fill = method)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("PrediXcan" = "#59268F", "DAPGW" = "#24A860", "MASHR" = "#FEF236", "UTMOST" = "#1265B2", "PERSIMMON" = "#F18032", "Cauchy" = "#E7002D")) +
  theme_classic() +
  xlab("") +
  ylab("Gene-trait associations (COVID19)") +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 14, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 10, family = "Helvetica", angle = 30, vjust = 0.6)) +
  theme(axis.text.y = element_text(face="bold", size = 10, family = "Helvetica"))
dev.off()


##Create gene-trait pair plot (MAGMA)##
B = c(2-1, 3-2, 2-1, 1-1, 173-172, 164-163)
A = c(1, 2, 1, 1, 172, 163)
data = as.data.frame(cbind(A, B))
rownames(data) = c("PrediXcan", "DAPGW", "MASHR", "UTMOST", "PERSIMMON", "Cauchy")
data = rbind(data, "gg" = 1:2)
data = t(data)
data = as.data.frame(data)
data <- reshape2::melt(data, "gg")

col = c("#59268F", "white", "#24A860", "white", "#FEF236", "white", "#1265B2", "white", "#F18032", "white", "#E7002D", "white")
tiff("/Users/chachritkhun/Desktop/plot1.tiff", units="in",  width=4, height=4, res=300)
ggplot(data=data, aes(x=variable, y=value)) +
  geom_bar(stat="identity", color="black", fill = col, size=0.25) +
  theme_classic() +
  xlab("") +
  ylab("Gene count") +
  #ggtitle("MAGMA") +
  #theme(plot.title = element_text(face="bold", size = 10, family = "Helvetica", hjust = 0.5)) +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 14, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 10, family = "Helvetica", angle = 30, vjust = 0.6)) +
  theme(axis.text.y = element_text(face="bold", size = 10, family = "Helvetica"))
dev.off()

##Create gene-trait pair plot (GWAS)##
B = c(2-1, 3-1, 2-0, 1-1, 173-171, 164-162)
A = c(1, 1, 0, 1, 171, 162)
data = as.data.frame(cbind(A, B))
rownames(data) = c("PrediXcan", "DAPGW", "MASHR", "UTMOST", "PERSIMMON", "Cauchy")
data = rbind(data, "gg" = 1:2)
data = t(data)
data = as.data.frame(data)
data <- reshape2::melt(data, "gg")

col = c("#59268F", "white", "#24A860", "white", "#FEF236", "white", "#1265B2", "white", "#F18032", "white", "#E7002D", "white")
tiff("/Users/chachritkhun/Desktop/plot2.tiff", units="in",  width=4, height=4, res=300)
ggplot(data=data, aes(x=variable, y=value)) +
  geom_bar(stat="identity", color="black", fill = col, size=0.25) +
  theme_classic() +
  xlab("") +
  ylab("Gene count") +
  #ggtitle("GWAS Clumping") +
  #theme(plot.title = element_text(face="bold", size = 10, family = "Helvetica", hjust = 0.5)) +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 14, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 12, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 10, family = "Helvetica", angle = 30, vjust = 0.6)) +
  theme(axis.text.y = element_text(face="bold", size = 10, family = "Helvetica"))
dev.off()


