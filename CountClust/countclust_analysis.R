library(maptpx)
library(CountClust)
library(singleCellRNASeqMouseDeng2014)
library(RColorBrewer)
library(data.table)
library(tibble)

gtex.meta_data <- read.table("/Users/poom/Desktop/countcluster/gtex_meta_data_final.txt", header = FALSE, sep = "\t")
gtex.meta_data$V3 = as.character(gtex.meta_data$V3)
gtex.meta_data <- gtex.meta_data %>% filter(V3 != "Bladder") %>% filter(V3 != "Cervix - Ectocervix") %>% filter(V3 != "Kidney - Cortex") %>% filter(V3 != "Cervix - Endocervix") %>% filter(V3 != "Fallopian Tube")

omega <- read.table("/Users/poom/Desktop/countcluster/omega_k20.txt", header = TRUE, sep = "\t")
rownames(omega) = omega$sample
omega <- omega %>% select(-sample)
colnames(omega) = 1:ncol(omega)
omega = omega[rownames(omega) %in% gtex.meta_data$V2,]

tissue_labels = gtex.meta_data$V3
docweights_per_tissue_mean <- apply(omega, 2,
                                    function(x) { tapply(x, tissue_labels, mean) })


##Draw structure plot##
annotation <- data.frame(
  sample_id = paste0("X", 1:length(tissue_labels)),
  tissue_label = factor(tissue_labels,
                        levels = rev(unique(tissue_labels) ) ) );

tiff("/Users/poom/Desktop/countclust_splot.tiff", units="in",  width=11, height=8, res=300)
StructureGGplot(omega = omega,
                annotation= annotation,
                yaxis_label = "",
                palette = c(RColorBrewer::brewer.pal(5, "Dark2"), RColorBrewer::brewer.pal(5, "Accent"), RColorBrewer::brewer.pal(5, "Greens"), RColorBrewer::brewer.pal(5, "Reds")),
                order_sample = TRUE,
                split_line = list(split_lwd = .4,
                                  split_col = "white"),
                axis_tick = list(axis_ticks_length = .1,
                                 axis_ticks_lwd_y = .1,
                                 axis_ticks_lwd_x = .1,
                                 axis_label_size = 7,
                                 axis_label_face = "bold"))
dev.off()

##Draw heat-map##
tiff("/Users/poom/Desktop/countclust_hmap.tiff", units="in",  width=13, height=8, res=300)
ordering <- heatmap(docweights_per_tissue_mean, distfun = function(x) dist(x, method = "manhattan"))$rowInd
dev.off()

dist = as.matrix(dist(docweights_per_tissue_mean, diag = TRUE, upper = TRUE, method = "euclidean"))
feng_list = c("Adrenal Gland", "Artery - Aorta", "Brain - Frontal Cortex (BA9)", "Brain - Hippocampus", "Breast - Mammary Tissue",
              "Cells - EBV-transformed lymphocytes", "Heart - Left Ventricle", "Liver", "Lung", "Muscle - Skeletal",
              "Ovary", "Pancreas", "Skin - Not Sun Exposed (Suprapubic)", "Small Intestine - Terminal Ileum", "Spleen")

lieberman_list = c("Breast - Mammary Tissue", "Cells - EBV-transformed lymphocytes", "Lung", "Skin - Not Sun Exposed (Suprapubic)")

pchic_list = c("Adipose - Subcutaneous", "Adrenal Gland", "Artery - Aorta", "Brain - Frontal Cortex (BA9)", "Brain - Hippocampus", 
               "Cells - EBV-transformed lymphocytes", "Colon - Sigmoid", "Esophagus - Gastroesophageal Junction", "Heart - Atrial Appendage",
               "Heart - Left Ventricle", "Liver", "Lung", "Muscle - Skeletal", "Ovary", "Pancreas", "Small Intestine - Terminal Ileum", "Spleen",
               "Stomach")

result = data.frame()
for (i in 1:nrow(dist)) {
  temp = t(as.matrix(colnames(dist)[order(dist[i,], decreasing = FALSE)]))
  
  index1 = vector()
  for (j in 1:length(feng_list)){
    index1_temp = which(temp == feng_list[j])
    index1 = append(index1, index1_temp)
  }
  
  index2 = vector()
  for (j in 1:length(lieberman_list)){
    index2_temp = which(temp == lieberman_list[j])
    index2 = append(index2, index2_temp)
  }
  
  index3 = vector()
  for (j in 1:length(pchic_list)){
    index3_temp = which(temp == pchic_list[j])
    index3 = append(index3, index3_temp)
  }
  
  result_temp = cbind("tissue" = temp[1], "feng_proxy" = feng_list[which(index1 == min(index1))], 
                      "lieberman_proxy" = lieberman_list[which(index2 == min(index2))],
                      "pchic_proxy" = pchic_list[which(index3 == min(index3))])
  result = rbind(result, result_temp)
}
write.table(result, "/Users/poom/Desktop/countcluster/result_final.txt", quote= FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

