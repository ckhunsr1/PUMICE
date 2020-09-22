library(dplyr)
library(ggrepel)
library(reshape2)

#######################################################################################################################################
df = read.table("/Users/chachritkhun/Desktop/TWAS_files/espresso_finemap_result2.txt", header = TRUE, na.string ='.', as.is=TRUE)
#######################################################################################################################################
##Plot boxplot##
method = c("PrediXcan", "DAPGW", "PUMICE", "UTMOST", "MASHR", "PUMICE+(MASHR)")
ex = df %>% select(pheno, predixcan, dapg, multi, utmost, mashr, cauchymam)
ex = melt(data = ex, id.vars = "pheno", measure.vars = colnames(ex)[2:7])
levels(ex$variable) = method

col1 = c("#59268F", "#24A860", "#F18032", "#1265B2", "#FEF236", "#FEF236")
col2 = c("#59268F", "#24A860", "#F18032", "#1265B2", "#FEF236", "#F18032")
tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=6, height=8, res=300)
ggplot(ex, aes(x=variable, y=value)) + 
  geom_boxplot_pattern(
    aes(pattern = variable, pattern_fill = variable),
    pattern = c('stripe'),
    pattern_alpha = c(0, 0, 0, 0, 0, 1),
    fill = col1,
    colour = 'black', pattern_density = 0.5, pattern_angle = 30, pattern_key_scale_factor = 2,
    lwd=0.5, outlier.shape  = NA) +
  scale_pattern_color_manual(values = col1) +
  scale_pattern_fill_manual(values = col2) +
  geom_point(data = data.frame(x = levels(ex$variable), y = c(1.99, 2.63, 2.05, 2.10, 2.25, 1.96)), aes(x=x, y=y), color = 'red', size = 2) +
  ylim(0.8, 4.5) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Average size of 90% credible set") +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 14, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 14, family = "Helvetica"))
dev.off()
#######################################################################################################################################
##MU vs. MMA##
df$col = "A"
df[df$cauchymu < df$cauchymam, "col"] = "B"
df[df$cauchymu > df$cauchymam, "col"] = "C"
df$col = as.factor(df$col)

##creat nudge vector##
diff = as.data.frame(cbind(df$pheno, df$cauchymu, df$cauchymam, df$cauchymam - df$cauchymu))
diff$V2 = as.numeric(as.character(diff$V2))
diff$V3 = as.numeric(as.character(diff$V3))
diff$V4 = as.numeric(as.character(diff$V4))
diff = diff[order(diff$V4), ]
row_select = c(1:3, (nrow(diff)-2):nrow(diff))
diff = diff[row_select,]

nudge_vec = c()
for (i in 1:nrow(diff)){
  if ( diff$V4[i] > 0 ) {
    temp = diff$V3[i] - 1
    nudge_vec = c(nudge_vec, temp)
  } else {
    temp = diff$V3[i] - 6
    nudge_vec = c(nudge_vec, temp)
  }
}

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=6, height=6, res=300)
ggplot(df, aes(x=cauchymam, y=cauchymu, label = pheno, color = col)) +
  geom_point(aes(color=col, shape=col, size = col)) +
  scale_color_manual(values = c("black", "#F18032", "#F18032")) +
  scale_shape_manual(values=c(20, 17, 18)) +
  scale_size_manual(values=c(3,3,4))+
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  xlim(0, 7) +
  ylim(0, 7) +
  geom_label_repel(data = diff %>% filter(V4 > 0), aes(x=V3, y=V2, label = V1), color = "black", nudge_x = nudge_vec[which(diff$V4 > 0)], size = 4,
                   box.padding = 0.5, point.padding = 0.5, force = 10,
                   segment.size  = 0.4, fontface = 'bold', direction = "y") +
  geom_label_repel(data = diff %>% filter(V4 < 0), aes(x=V3, y=V2, label = V1), color = "black", nudge_x = nudge_vec[which(diff$V4 < 0)],size = 4,
                   box.padding = 0.5, point.padding = 0.5, force = 20,
                   segment.size  = 0.4, fontface = 'bold', direction = "y") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "PUMICE+(MASHR)", y = "PUMICE+(UTMOST)") +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 14, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 14, family = "Helvetica"))
dev.off()

#######################################################################################################################################
##MMA vs. DMA##
df$col = "A"
df[df$cauchymam < df$cauchymad, "col"] = "B"
df[df$cauchymam > df$cauchymad, "col"] = "C"
df$col = as.factor(df$col)

##creat nudge vector##
diff = as.data.frame(cbind(df$pheno, df$cauchymam, df$cauchymad, df$cauchymad - df$cauchymam))
diff$V2 = as.numeric(as.character(diff$V2))
diff$V3 = as.numeric(as.character(diff$V3))
diff$V4 = as.numeric(as.character(diff$V4))
diff = diff[order(diff$V4), ]
row_select = c(1:3, (nrow(diff)-2):nrow(diff))
diff = diff[row_select,]

nudge_vec = c()
for (i in 1:nrow(diff)){
  if ( diff$V4[i] > 0 ) {
    temp = diff$V2[i] - 6
    nudge_vec = c(nudge_vec, temp)
  } else {
    temp = diff$V2[i] - 1
    nudge_vec = c(nudge_vec, temp)
  }
}

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=6, height=6, res=300)
ggplot(df, aes(x=cauchymam, y=cauchymad, label = pheno, color = col)) +
  geom_point(aes(color=col, shape = col, size = col)) +
  scale_color_manual(values = c("black", "#F18032", "#24A860")) +
  scale_shape_manual(values=c(20, 18, 18)) +
  scale_size_manual(values=c(3,4,4))+
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  xlim(0, 7) +
  ylim(0, 7) +
  geom_label_repel(data = diff %>% filter(V4 > 0), aes(x=V2, y=V3, label = V1), color = "black", nudge_x = nudge_vec[which(diff$V4 > 0)], size = 4,
                   box.padding = 0.5, point.padding = 0.5, force = 10,
                   segment.size  = 0.4, fontface = 'bold', direction = "y") +
  geom_label_repel(data = diff %>% filter(V4 < 0), aes(x=V2, y=V3, label = V1), color = "black", nudge_x = nudge_vec[which(diff$V4 < 0)],size = 4,
                   box.padding = 0.5, point.padding = 0.5, force = 20,
                   segment.size  = 0.4, fontface = 'bold', direction = "y") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "PUMICE+(MASHR)", y = "DAPGW+(MASHR)") +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 14, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 14, family = "Helvetica"))
dev.off()

#######################################################################################################################################
##MMA vs. PMA##
df$col = "A"
df[df$cauchymam < df$cauchymap, "col"] = "B"
df[df$cauchymam > df$cauchymap, "col"] = "C"
df$col = as.factor(df$col)

##creat nudge vector##
diff = as.data.frame(cbind(df$pheno, df$cauchymam, df$cauchymap, df$cauchymap - df$cauchymam))
diff$V2 = as.numeric(as.character(diff$V2))
diff$V3 = as.numeric(as.character(diff$V3))
diff$V4 = as.numeric(as.character(diff$V4))
diff = diff[order(diff$V4), ]
row_select = c(1:3, (nrow(diff)-2):nrow(diff))
diff = diff[row_select,]

nudge_vec = c()
for (i in 1:nrow(diff)){
  if ( diff$V4[i] > 0 ) {
    temp = diff$V2[i] - 4
    nudge_vec = c(nudge_vec, temp)
  } else {
    temp = diff$V2[i] - 1
    nudge_vec = c(nudge_vec, temp)
  }
}

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=6, height=6, res=300)
ggplot(df, aes(x=cauchymam, y=cauchymap, label = pheno, color = col)) +
  geom_point(aes(color=col, size = col, shape = col)) +
  scale_color_manual(values = c("black", "#F18032", "#59268F")) +
  scale_shape_manual(values=c(20, 18, 18)) +
  scale_size_manual(values=c(3,4,4))+
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  xlim(0, 7) +
  ylim(0, 7) +
  geom_label_repel(data = diff %>% filter(V4 > 0), aes(x=V2, y=V3, label = V1), color = "black", nudge_x = nudge_vec[which(diff$V4 > 0)], size = 4,
                   box.padding = 0.5, point.padding = 0.5, force = 10,
                   segment.size  = 0.4, fontface = 'bold', direction = "y") +
  geom_label_repel(data = diff %>% filter(V4 < 0), aes(x=V2, y=V3, label = V1), color = "black", nudge_x = nudge_vec[which(diff$V4 < 0)],size = 4,
                   box.padding = 0.5, point.padding = 0.5, force = 20,
                   segment.size  = 0.4, fontface = 'bold', direction = "y") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "PUMICE+(MASHR)", y = "PrediXcan+(MASHR)") +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 14, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 14, family = "Helvetica"))
dev.off()

#######################################################################################################################################
##MMA vs. PU##
df$col = "A"
df[df$cauchymam < df$cauchypu, "col"] = "B"
df[df$cauchymam > df$cauchypu, "col"] = "C"
df$col = as.factor(df$col)

##creat nudge vector##
diff = as.data.frame(cbind(df$pheno, df$cauchymam, df$cauchypu, df$cauchypu - df$cauchymam))
diff$V2 = as.numeric(as.character(diff$V2))
diff$V3 = as.numeric(as.character(diff$V3))
diff$V4 = as.numeric(as.character(diff$V4))
diff = diff[order(diff$V4), ]
row_select = c(1:3, (nrow(diff)-2):nrow(diff))
diff = diff[row_select,]

nudge_vec = c()
for (i in 1:nrow(diff)){
  if ( diff$V4[i] > 0 ) {
    temp = diff$V2[i] - 4
    nudge_vec = c(nudge_vec, temp)
  } else {
    temp = diff$V2[i] - 1
    nudge_vec = c(nudge_vec, temp)
  }
}

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=6, height=6, res=300)
ggplot(df, aes(x=cauchymam, y=cauchypu, label = pheno, color = col)) +
  geom_point(aes(color=col, size = col, shape = col)) +
  scale_color_manual(values = c("black", "#F18032", "#59268F")) +
  scale_shape_manual(values=c(20, 18, 17)) +
  scale_size_manual(values=c(3,4,3))+
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  xlim(0, 7) +
  ylim(0, 7) +
  geom_label_repel(data = diff %>% filter(V4 > 0), aes(x=V2, y=V3, label = V1), color = "black", nudge_x = nudge_vec[which(diff$V4 > 0)], size = 4,
                   box.padding = 0.5, point.padding = 0.5, force = 10,
                   segment.size  = 0.4, fontface = 'bold', direction = "y") +
  geom_label_repel(data = diff %>% filter(V4 < 0), aes(x=V2, y=V3, label = V1), color = "black", nudge_x = nudge_vec[which(diff$V4 < 0)],size = 4,
                   box.padding = 0.5, point.padding = 0.5, force = 20,
                   segment.size  = 0.4, fontface = 'bold', direction = "y") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "PUMICE+(MASHR)", y = "PrediXcan+(UTMOST)") +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 14, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 14, family = "Helvetica"))
dev.off()
#######################################################################################################################################
##MMA vs. DU##
df$col = "A"
df[df$cauchymam < df$cauchydu, "col"] = "B"
df[df$cauchymam > df$cauchydu, "col"] = "C"

df$col = as.factor(df$col)

##creat nudge vector##
diff = as.data.frame(cbind(df$pheno, df$cauchymam, df$cauchydu, df$cauchydu - df$cauchymam))
diff$V2 = as.numeric(as.character(diff$V2))
diff$V3 = as.numeric(as.character(diff$V3))
diff$V4 = as.numeric(as.character(diff$V4))
diff = diff[order(diff$V4), ]
row_select = c(1:3, (nrow(diff)-2):nrow(diff))
diff = diff[row_select,]

nudge_vec = c()
for (i in 1:nrow(diff)){
  if ( diff$V4[i] > 0 ) {
    temp = diff$V2[i] - 3
    nudge_vec = c(nudge_vec, temp)
  } else {
    temp = diff$V2[i] -1
    nudge_vec = c(nudge_vec, temp)
  }
}

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=6, height=6, res=300)
ggplot(df, aes(x=cauchymam, y=cauchydu, label = pheno, color = col)) +
  geom_point(aes(color=col, size = col, shape = col)) +
  scale_color_manual(values = c("black", "#F18032", "#24A860")) +
  scale_shape_manual(values=c(20, 18, 17)) +
  scale_size_manual(values=c(3,4,3))+
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  xlim(0, 7) +
  ylim(0, 7) +
  geom_label_repel(data = diff %>% filter(V4 > 0), aes(x=V2, y=V3, label = V1), color = "black", nudge_x = nudge_vec[which(diff$V4 > 0)], size = 4,
                   box.padding = 0.5, point.padding = 0.5, force = 10,
                   segment.size  = 0.4, fontface = 'bold', direction = "y") +
  geom_label_repel(data = diff %>% filter(V4 < 0), aes(x=V2, y=V3, label = V1), color = "black", nudge_x = nudge_vec[which(diff$V4 < 0)],size = 4,
                   box.padding = 0.5, point.padding = 0.5, force = 10,
                   segment.size  = 0.4, fontface = 'bold', direction = "y") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "PUMICE+(MASHR)", y = "DAPGW+(UTMOST)") +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 14, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 14, family = "Helvetica"))
dev.off()
#######################################################################################################################################
##MMA vs. PERSIMMON##
df$col = "A"
df[df$cauchymam < df$multi, "col"] = "B"
df[df$cauchymam > df$multi, "col"] = "C"
df$col = as.factor(df$col)

##creat nudge vector##
diff = as.data.frame(cbind(df$pheno, df$cauchymam, df$multi, df$multi - df$cauchymam))
diff$V2 = as.numeric(as.character(diff$V2))
diff$V3 = as.numeric(as.character(diff$V3))
diff$V4 = as.numeric(as.character(diff$V4))
diff = diff[order(diff$V4), ]
row_select = c(1:3, (nrow(diff)-2):nrow(diff))
diff = diff[row_select,]

nudge_vec = c()
for (i in 1:nrow(diff)){
  if ( diff$V4[i] > 0 ) {
    temp = diff$V2[i] - 5
    nudge_vec = c(nudge_vec, temp)
  } else {
    temp = diff$V2[i]
    nudge_vec = c(nudge_vec, temp)
  }
}

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=6, height=6, res=300)
ggplot(df, aes(x=cauchymam, y=multi, label = pheno, color = col)) +
  geom_point(aes(color=col, size = col, shape = col)) +
  scale_color_manual(values = c("black", "#F18032", "#F18032")) +
  scale_shape_manual(values=c(20, 18, 20)) +
  scale_size_manual(values=c(3,4,3))+
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  xlim(0, 7) +
  ylim(0, 7) +
  geom_label_repel(data = diff %>% filter(V4 > 0), aes(x=V2, y=V3, label = V1), color = "black", nudge_x = nudge_vec[which(diff$V4 > 0)], size = 4,
                   box.padding = 0.5, point.padding = 0.5, force = 10,
                   segment.size  = 0.4, fontface = 'bold', direction = "y") +
  geom_label_repel(data = diff %>% filter(V4 < 0), aes(x=V2, y=V3, label = V1), color = "black", nudge_x = nudge_vec[which(diff$V4 < 0)],size = 4,
                   box.padding = 0.5, point.padding = 0.5, force = 10,
                   segment.size  = 0.4, fontface = 'bold', direction = "y") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "PUMICE+(MASHR)", y = "PUMICE") +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 14, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 14, family = "Helvetica"))
dev.off()
#######################################################################################################################################
##MMA vs. PrediXcan##
df$col = "A"
df[df$cauchymam < df$predixcan, "col"] = "B"
df[df$cauchymam > df$predixcan, "col"] = "C"
df$col = as.factor(df$col)

##creat nudge vector##
diff = as.data.frame(cbind(df$pheno, df$cauchymam, df$predixcan, df$predixcan - df$cauchymam))
diff$V2 = as.numeric(as.character(diff$V2))
diff$V3 = as.numeric(as.character(diff$V3))
diff$V4 = as.numeric(as.character(diff$V4))
diff = diff[order(diff$V4), ]
row_select = c(1:3, (nrow(diff)-2):nrow(diff))
diff = diff[row_select,]

nudge_vec = c()
for (i in 1:nrow(diff)){
  if ( diff$V4[i] > 0 ) {
    temp = diff$V2[i] - 5
    nudge_vec = c(nudge_vec, temp)
  } else {
    temp = diff$V2[i] - 2
    nudge_vec = c(nudge_vec, temp)
  }
}

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=6, height=6, res=300)
ggplot(df, aes(x=cauchymam, y=predixcan, label = pheno, color = col)) +
  geom_point(aes(color=col, size = col, shape = col)) +
  scale_color_manual(values = c("black", "#F18032", "#59268F")) +
  scale_shape_manual(values=c(20, 18, 20)) +
  scale_size_manual(values=c(3,4,3))+
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  xlim(0, 7) +
  ylim(0, 7) +
  geom_label_repel(data = diff %>% filter(V4 > 0), aes(x=V2, y=V3, label = V1), color = "black", nudge_x = nudge_vec[which(diff$V4 > 0)], size = 4,
                   box.padding = 0.5, point.padding = 0.5, force = 10,
                   segment.size  = 0.4, fontface = 'bold', direction = "y") +
  geom_label_repel(data = diff %>% filter(V4 < 0), aes(x=V2, y=V3, label = V1), color = "black", nudge_x = nudge_vec[which(diff$V4 < 0)],size = 4,
                   box.padding = 0.5, point.padding = 0.5, force = 10,
                   segment.size  = 0.4, fontface = 'bold', direction = "y") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "PUMICE+(MASHR)", y = "PrediXcan") +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 14, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 14, family = "Helvetica"))
dev.off()
#######################################################################################################################################
##MMA vs. UTMOST##
df$col = "A"
df[df$cauchymam < df$utmost, "col"] = "B"
df[df$cauchymam > df$utmost, "col"] = "C"
df$col = as.factor(df$col)

##creat nudge vector##
diff = as.data.frame(cbind(df$pheno, df$cauchymam, df$utmost, df$utmost - df$cauchymam))
diff$V2 = as.numeric(as.character(diff$V2))
diff$V3 = as.numeric(as.character(diff$V3))
diff$V4 = as.numeric(as.character(diff$V4))
diff = diff[order(diff$V4), ]
row_select = c(1:3, (nrow(diff)-2):nrow(diff))
diff = diff[row_select,]

nudge_vec = c()
for (i in 1:nrow(diff)){
  if ( diff$V4[i] > 0 ) {
    temp = diff$V2[i] - 4
    nudge_vec = c(nudge_vec, temp)
  } else {
    temp = diff$V2[i] 
    nudge_vec = c(nudge_vec, temp)
  }
}

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=6, height=6, res=300)
ggplot(df, aes(x=cauchymam, y=utmost, label = pheno, color = col)) +
  geom_point(aes(color=col, shape = col, size = col)) +
  scale_color_manual(values = c("black", "#F18032", "#1265B2")) +
  scale_shape_manual(values=c(20, 18, 20)) +
  scale_size_manual(values=c(3,4,3))+
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  xlim(0, 7) +
  ylim(0, 7) +
  geom_label_repel(data = diff %>% filter(V4 > 0), aes(x=V2, y=V3, label = V1), color = "black", nudge_x = nudge_vec[which(diff$V4 > 0)],size = 4,
                   box.padding = 0.5, point.padding = 0.5, force = 10,
                   segment.size  = 0.4, fontface = 'bold', direction = "y") +
  geom_label_repel(data = diff %>% filter(V4 < 0), aes(x=V2, y=V3, label = V1), color = "black", nudge_x = nudge_vec[which(diff$V4 < 0)],size = 4,
                   box.padding = 0.5, point.padding = 0.5, force = 10,
                   segment.size  = 0.4, fontface = 'bold', direction = "y") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "PUMICE+(MASHR)", y = "UTMOST") +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 14, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 14, family = "Helvetica"))
dev.off()
#######################################################################################################################################
##MMA vs. dapgw##
df$col = "A"
df[df$cauchymam < df$dapg, "col"] = "B"
df[df$cauchymam > df$dapg, "col"] = "C"
df$col = as.factor(df$col)

##creat nudge vector##
diff = as.data.frame(cbind(df$pheno, df$cauchymam, df$dapg, df$dapg - df$cauchymam))
diff$V2 = as.numeric(as.character(diff$V2))
diff$V3 = as.numeric(as.character(diff$V3))
diff$V4 = as.numeric(as.character(diff$V4))
diff = diff[order(diff$V4), ]
row_select = c(1:3, (nrow(diff)-2):nrow(diff))
diff = diff[row_select,]

nudge_vec = c()
for (i in 1:nrow(diff)){
  if ( diff$V4[i] > 0 ) {
    temp = diff$V2[i] - 4
    nudge_vec = c(nudge_vec, temp)
  } else {
    temp = diff$V2[i] - 1
    nudge_vec = c(nudge_vec, temp)
  }
}

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=6, height=6, res=300)
ggplot(df, aes(x=cauchymam, y=dapg, label = pheno, color = col)) +
  geom_point(aes(color=col, size = col, shape = col)) +
  scale_color_manual(values = c("black", "#F18032", "#24A860")) +
  scale_shape_manual(values=c(20, 18, 20)) +
  scale_size_manual(values=c(3,4,3))+
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  xlim(0, 7) +
  ylim(0, 7) +
  geom_label_repel(data = diff %>% filter(V4 > 0), aes(x=V2, y=V3, label = V1), color = "black", nudge_x = nudge_vec[which(diff$V4 > 0)], size = 4,
                   box.padding = 0.5, point.padding = 0.5, force = 10,
                   segment.size  = 0.4, fontface = 'bold', direction = "y") +
  geom_label_repel(data = diff %>% filter(V4 < 0), aes(x=V2, y=V3, label = V1), color = "black", nudge_x = nudge_vec[which(diff$V4 < 0)],size = 4,
                   box.padding = 0.5, point.padding = 0.5, force = 10,
                   segment.size  = 0.4, fontface = 'bold', direction = "y") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "PUMICE+(MASHR)", y = "DAPGW") +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 14, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 14, family = "Helvetica"))
dev.off()
#######################################################################################################################################
##MMA vs. mashr##
df$col = "A"
df[df$cauchymam < df$mashr, "col"] = "B"
df[df$cauchymam > df$mashr, "col"] = "C"
df$col = as.factor(df$col)

##creat nudge vector##
diff = as.data.frame(cbind(df$pheno, df$cauchymam, df$mashr, df$mashr - df$cauchymam))
diff$V2 = as.numeric(as.character(diff$V2))
diff$V3 = as.numeric(as.character(diff$V3))
diff$V4 = as.numeric(as.character(diff$V4))
diff = diff[order(diff$V4), ]
row_select = c(1:3, (nrow(diff)-2):nrow(diff))
diff = diff[row_select,]


nudge_vec = c()
for (i in 1:nrow(diff)){
  if ( diff$V4[i] > 0 ) {
    temp = diff$V2[i] - 5
    nudge_vec = c(nudge_vec, temp)
  } else {
    temp = diff$V2[i] 
    nudge_vec = c(nudge_vec, temp)
  }
}

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=6, height=6, res=300)
ggplot(df, aes(x=cauchymam, y=mashr, label = pheno, color = col)) +
  geom_point(aes(color=col, size = col, shape = col)) +
  scale_color_manual(values = c("black", "#F18032", "#FEF236")) +
  scale_shape_manual(values=c(20, 18, 20)) +
  scale_size_manual(values=c(3,4,3))+
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  xlim(0, 7) +
  ylim(0, 7) +
  geom_label_repel(data = diff %>% filter(V4 > 0), aes(x=V2, y=V3, label = V1), color = "black", nudge_x = nudge_vec[which(diff$V4 > 0)], size = 4,
                   box.padding = 0.5, point.padding = 0.5, force = 10,
                   segment.size  = 0.4, fontface = 'bold', direction = "y") +
  geom_label_repel(data = diff %>% filter(V4 < 0), aes(x=V2, y=V3, label = V1), color = "black", nudge_x = nudge_vec[which(diff$V4 < 0)],size = 4,
                   box.padding = 0.5, point.padding = 0.5, force = 10,
                   segment.size  = 0.4, fontface = 'bold', direction = "y") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "PUMICE+(MASHR)", y = "MASHR") +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 16, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 14, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 14, family = "Helvetica"))
dev.off()
#######################################################################################################################################
