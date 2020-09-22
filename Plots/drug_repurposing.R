library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

##################################################################################################################
df = read.table("/Users/chachritkhun/Desktop/TWAS_files/result_all_tissues_new.txt", header = TRUE)
df$name_all = paste(df$method, df$name_d, sep = "_")
x = df %>% select(method, name_d, name_all, mean) %>% group_by(name_all) %>% dplyr::slice(which.min(mean))
inc = (x %>% group_by(name_d) %>% tally() %>% filter(n == length(unique(x$method))))$name_d
x = as.data.frame(x[x$name_d %in% inc, ]) %>% select(method, name_d, mean)
colnames(x)[2:3] = c("variable", "value")

##Determine the order of drug based on hierachical clustering##
mat = matrix(0, length(unique(x$variable)), length(unique(x$method)))
for ( i in 1:length(unique(x$method)) ) {
  meth = unique(x$method)[i]
  for ( j in 1:length(unique(x$variable)) ) {
    drug = unique(x$variable)[j]
    temp = x %>% filter(method == meth) %>% filter(variable == drug)
    mat[j, i] = temp$value
  }
}
y = as.data.frame(mat)
colnames(y) = unique(x$method)
rownames(y) = unique(x$variable)

mat = scale(mat)
res.dist <- dist(mat, method = "euclidean")
#res.dist = as.matrix(res.dist)
res.hc <- hclust(d = res.dist, method = "ward.D2")

x$variable = factor(x$variable, levels = unique(x$variable)[res.hc$order])

z = x %>% group_by(variable) %>% mutate(neg = sum(value < 0)) %>% mutate(pos = sum(value > 0)) 
z$value = abs(z$value)
z = z %>% filter(value > 0) %>% filter(neg >= 0)%>% group_by(variable) %>% tally() %>% filter(n > 0)
#z = z %>% filter(neg > 0) %>% group_by(variable) %>% dplyr::slice(which.max(value)) %>% filter(value > 10)
#z = z %>% filter(neg >= 0) %>% select(-c("neg")) %>% group_by(variable) %>% dplyr::slice(which.max(value)) %>% filter(value > 20)
y = y[rownames(y) %in% unique(z$variable), ]
x = x[x$variable %in% unique(z$variable), ]


m = x %>% group_by(variable) %>% mutate(mean = mean(value))
nrow(m %>% filter(mean < 0) %>% filter(value < mean) %>% filter(method == "cauchy_combined"))
nrow(m %>% filter(mean < 0) %>% filter(value < mean) %>% filter(method == "cauchy_combined_mashr"))
nrow(m %>% filter(mean < 0) %>% filter(value < mean) %>% filter(method == "cauchy_dapg_mashr"))
nrow(m %>% filter(mean < 0) %>% filter(value < mean) %>% filter(method == "cauchy_predixcan_mashr"))
nrow(m %>% filter(mean < 0) %>% filter(value < mean) %>% filter(method == "cauchy_dapg_utmost"))
nrow(m %>% filter(mean < 0) %>% filter(value < mean) %>% filter(method == "cauchy_predixcan_utmost"))
nrow(m %>% filter(mean < 0) %>% filter(value < mean) %>% filter(method == "utmost"))
nrow(m %>% filter(mean < 0) %>% filter(value < mean) %>% filter(method == "multi-omics"))
nrow(m %>% filter(mean < 0) %>% filter(value < mean) %>% filter(method == "dapg"))
nrow(m %>% filter(mean < 0) %>% filter(value < mean) %>% filter(method == "mashr"))
nrow(m %>% filter(mean < 0) %>% filter(value < mean) %>% filter(method == "predixcan"))

mean((x %>% filter(method == "cauchy_combined"))$value)
mean((x %>% filter(method == "cauchy_combined_mashr"))$value)
mean((x %>% filter(method == "cauchy_dapg_mashr"))$value)
mean((x %>% filter(method == "cauchy_predixcan_mashr"))$value)
mean((x %>% filter(method == "cauchy_dapg_utmost"))$value)
mean((x %>% filter(method == "cauchy_predixcan_utmost"))$value)
mean((x %>% filter(method == "utmost"))$value)
mean((x %>% filter(method == "multi-omics"))$value)
mean((x %>% filter(method == "dapg"))$value)
mean((x %>% filter(method == "mashr"))$value)
mean((x %>% filter(method == "predixcan"))$value)


##################################################################################################################
x$cut = as.factor(as.numeric(cut(x$value, c(-100, -80, -60, -40, -20, 0, 20, 40, 60))))
#x$cut = as.factor(as.numeric(cut(x$value, c(-50, -40, -30, -20, -10, 0, 10, 20, 30))))
levels(x$cut) <- c("A","B","C", "D", "E", "F", "G", "H")

idx = c(1:5,7:9)
jcolors = brewer.pal(n = 11, name = "RdBu")[idx]
names(jcolors) = c("A","B","C", "D", "E", "F", "G", "H")

x$method = factor(x$method, levels = rev(c("predixcan", "dapg", "multi-omics", "mashr", "utmost", "cauchy_predixcan_utmost", "cauchy_dapg_utmost", "cauchy_combined", "cauchy_predixcan_mashr","cauchy_dapg_mashr", "cauchy_combined_mashr")))
levels(x$method) <- rev(c("PrediXcan", "DAPGW", "PUMICE", "MASHR", "UTMOST", "PrediXcan+(UTMOST)", "DAPGW+(UTMOST)", "PUMICE+(UTMOST)", "PrediXcan+(MASHR)", "DAPGW+(MASHR)", "PUMICE+(MASHR)"))
##################################################################################################################

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=12, height=3, res=300)
ggplot(x, aes(variable, method)) +
  geom_point(pch=15, aes(color=cut, size = abs(value))) +
  scale_size_continuous(range = c(0, 4)) +
  #scale_colour_gradient2(low = "#B72333", mid = "#F8F8F8", high = "#306EAF") +
  geom_vline(xintercept=seq(0, nrow(x) - 1, 1)+.5,color="black") +
  geom_hline(yintercept=seq(0, length(unique(x$method)) + 1, 1)+.5,color="black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.45, hjust=1)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  ggtitle("Computational drug repurposing") +
  theme(plot.title = element_text(face="bold", size = 12, family = "Helvetica", hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 8, family = "Helvetica")) +
  theme(axis.text.y = element_text(size = 8, family = "Helvetica")) +
  theme(legend.key.size = unit(0.1, "cm")) +
  #scale_color_discrete(name = "Drug effect", labels = c("Inducing", "", "", "", "", "Reversing")) +
  scale_color_manual(values=jcolors)
dev.off()

