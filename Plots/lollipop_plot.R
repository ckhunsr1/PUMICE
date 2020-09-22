library(ggplot2)
library(dplyr)

##Plot for delta R2##
data = read.table("/Users/chachritkhun/Desktop/TWAS_files/internal_validation/iv_performance.txt", header = TRUE)

##Non-significant tissue when compared to utmost##
#idx = c(17, 19, 28, 29, 30, 31, 34, 35, 36, 37, 38, 39, 40, 41, 43, 44, 45, 46, 47, 48)
tissue_list = as.character(data$tissue)
#tissue_list[idx] = paste("# ", tissue_list[idx], paste = "")
data$tissue = factor(tissue_list, levels = tissue_list)

##Create lollipop plot##
col= c(  rep(c("#7CC10D"), 12), rep(c("#C10D7C"), 11), rep(c("#0D7CC1"), 25))
tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=6, height=6, res=300)
ggplot(data, aes(x=tissue, y=mean_mp)) +
  geom_segment( aes(x=tissue, xend=tissue, y=0, yend=mean_mp), color= col) +
  geom_point( color= col, size=2) +
  theme_classic() +
  coord_flip() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab(bquote(bold("Average"~Delta*"R"^"2"))) +
  ggtitle("PUMICE - PrediXcan") +
  theme(plot.title = element_text(face="bold", size = 10, family = "Helvetica", hjust = 0.4)) +
  theme(axis.text.x = element_text(size = 8, family = "Helvetica")) +
  theme(axis.text.y = element_text(size = 8, family = "Helvetica"))  +
  theme(axis.title.x = element_text(face="bold", size = 10, family = "Helvetica")) 
dev.off()

##Plot for Percent change of mean R2##
data = read.table("/Users/chachritkhun/Desktop/TWAS_files/internal_validation/iv_per_performance.txt", header = TRUE)
tissue_list = as.character(data$tissue)
data$tissue = factor(tissue_list, levels = tissue_list)
data$gain_mp = 100*data$gain_mp
data$gain_mu = 100*data$gain_mu
data$gain_up = 100*data$gain_up

##Create lollipop plot##
col= c(  rep(c("#7CC10D"), 12), rep(c("#C10D7C"), 11), rep(c("#0D7CC1"), 25))
tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=6, height=6, res=300)
ggplot(data, aes(x=tissue, y=gain_mp)) +
  geom_segment( aes(x=tissue, xend=tissue, y=0, yend=gain_mp), color= col) +
  geom_point( color= col, size=2) +
  theme_classic() +
  coord_flip() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab(bquote(bold("% change of"~Delta*"(mean R"^"2"*")"))) +
  ggtitle("PUMICE - PrediXcan") +
  theme(plot.title = element_text(face="bold", size = 10, family = "Helvetica", hjust = 0.4)) +
  theme(axis.text.x = element_text(size = 8, family = "Helvetica")) +
  theme(axis.text.y = element_text(size = 8, family = "Helvetica"))  +
  theme(axis.title.x = element_text(face="bold", size = 10, family = "Helvetica")) 
dev.off()

##Create Cleaveland dot plot##
data = read.table("/Users/chachritkhun/Desktop/TWAS_files/internal_validation/iv_sig_count.txt", header = TRUE)

##Non-significant tissue when compared to utmost##
#idx = c(17, 19, 28, 29, 30, 31, 34, 35, 36, 37, 38, 39, 40, 41, 43, 44, 45, 46, 47, 48)
tissue_list = as.character(data$tissue)
#tissue_list[idx] = paste("# ", tissue_list[idx], paste = "")
data$tissue = factor(tissue_list, levels = tissue_list)

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=6, height=6, res=300)
ggplot(data) +
  geom_segment( aes(x=tissue, xend=tissue, y=sig_count_p, yend=sig_count_m), color=col, ) +
  geom_point( aes(x=tissue, y=sig_count_p), color=col, size=2, shape = 19 ) +
  geom_point( aes(x=tissue, y=sig_count_m), color=col, size=2, shape = 17) +
  coord_flip()+
  theme_classic() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  xlab("") +
  ylab("Number of significant models") +
  ggtitle("PUMICE - PrediXcan") +
  theme(plot.title = element_text(face="bold", size = 10, family = "Helvetica", hjust = 0.4)) +
  theme(axis.text.x = element_text(size = 8, family = "Helvetica")) +
  theme(axis.text.y = element_text(size = 8, family = "Helvetica"))  +
  theme(axis.title.x = element_text(face="bold", size = 10, family = "Helvetica"))
dev.off()
