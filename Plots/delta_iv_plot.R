library(dplyr)
library(ggplot2)
library(ggforce)
library(reshape2)
library(ggpubr)
###########################################################################################################
tissue = "Cells_EBV-transformed_lymphocytes"
df = read.table(paste("/Users/chachritkhun/Desktop/TWAS_files/prediction_performance/", tissue, "_GTEx_cv_result.txt", sep = ""), header = TRUE, na.string ='.', as.is=TRUE)
#df$delta = (df$pear_avg_t.x)^2 - (df$pear_avg_t.y)^2
method_list = c("mp", "me", "mbl", "mbs", "md", "mmc", "mu", "mdap")

df$type = factor(df$type, levels = method_list)
#x = df %>% filter(pear_avg_t.x == 0 & pear_avg_t.y == 0)
#df = setdiff(df, x)
###########################################################################################################
df_iv = df %>% group_by(type) %>% mutate(MEANEXP = mean(delta, na.rm = T))

iv_names <- c(
  `mp` = "PrediXcan",
  `me` = "EpiXcan",
  `mbl` = "BLUP",
  `mbs` = "BSLMM",
  `md` = "DPR(VB)",
  `mmc` = "DPR(MCMC)",
  `mu` = "UTMOST",
  `mdap` = "DAPGW"
)
##Drawing delta
tiff("/Users/chachritkhun/Desktop/calc_r1.tiff", units="in",  width=10, height=12, res=300)
ggplot(data = df_iv, aes(x=delta)) + 
  geom_density(alpha = 0.8,fill="#87CFD6",colour="#87CFD6") +
  facet_wrap( ~ type, nrow = 8, strip.position="right", labeller = as_labeller(iv_names)) +
  xlim(-1, 1) +
  geom_vline(aes(xintercept=MEANEXP), color="blue") +
  geom_vline(xintercept=0, linetype = "dashed") +
  theme_classic() +
  theme(strip.background = element_rect(color="black", fill="#87CFD6", size=1.5, linetype="solid")) +
  theme(strip.text = element_text(face="bold", size = 12, family = "Helvetica")) +
  labs(title="Imputation performance in cross-validation dataset (GTEx)") +
  xlab(expression(paste(Delta, "R", sep = ""))) +
  ylab("Density") +
  theme(plot.title = element_text(hjust =0.5, face="bold", size = 24, family = "Helvetica", margin = margin(b = 10))) +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 20, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 20, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 12, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 12, family = "Helvetica")) 
dev.off()

for (method in method_list) {
  print((wilcox.test( (df_iv %>% filter(type == method))$delta, mu = 0, alternative = "greater")))
}
table(df_iv$type)
###########################################################################################################
tissue = "Brain_Frontal_Cortex_BA9"
ext = "CMC"
df = read.table(paste("/Users/chachritkhun/Desktop/TWAS_files/prediction_performance/", tissue, "_GTEx_", ext, "_result.txt", sep = ""), header = TRUE, na.string ='.', as.is=TRUE)
#df$delta = (df$r.x) - (df$r.y)
method_list = c("mp", "me", "mbl", "mbs", "md", "mmc", "mu", "mdap", "mmash")

df$type = factor(df$type, levels = method_list)
#x = df %>% filter(r.x == 0 & r.y == 0)
#df = setdiff(df, x)
###########################################################################################################

df_ext = df %>% group_by(type) %>% mutate(MEANEXP = mean(delta, na.rm = T))

ext_names <- c(
  `mp` = "PrediXcan",
  `me` = "EpiXcan",
  `mbl` = "BLUP",
  `mbs` = "BSLMM",
  `md` = "DPR(VB)",
  `mmc` = "DPR(MCMC)",
  `mu` = "UTMOST",
  `mdap` = "DAPGW",
  `mmash` = "MASHR"
)
##Drawing delta instead
tiff("/Users/chachritkhun/Desktop/calc_r2.tiff", units="in",  width=10, height=12, res=300)
ggplot(data = df_ext, aes(x=delta)) + 
  geom_density(alpha = 0.8,fill="#F0D483",colour="#F0D483") +
  facet_wrap( ~ type, nrow = 9, strip.position="right", labeller = as_labeller(ext_names)) +
  xlim(-1, 1) +
  geom_vline(aes(xintercept=MEANEXP), color="red") +
  geom_vline(xintercept=0, linetype = "dashed") +
  theme_classic() +
  theme(strip.background = element_rect(color="black", fill="#F0D483", size=1.5, linetype="solid")) +
  theme(strip.text = element_text(face="bold", size = 12, family = "Helvetica")) +
  labs(title= paste("Imputation performance in ", ext," cohort (union)", sep = "")) +
  xlab(expression(paste(Delta, "R", sep = ""))) +
  ylab("Density") +
  theme(plot.title = element_text(hjust =0.5, face="bold", size = 21, family = "Helvetica", margin = margin(b = 10))) +
  theme(axis.title.x = element_text(vjust=-0.7, face="bold", size = 20, family = "Helvetica")) +
  theme(axis.title.y = element_text(vjust=2, face="bold", size = 20, family = "Helvetica")) +
  theme(axis.text.x = element_text(face="bold", size = 12, family = "Helvetica")) +
  theme(axis.text.y = element_text(face="bold", size = 12, family = "Helvetica")) 
dev.off()

for (method in method_list) {
  print((wilcox.test( (df_ext %>% filter(type == method))$delta, mu = 0, alternative = "greater")))
}
table(df_ext$type)

###########################################################################################################
