library(maptpx)
library(CountClust)
library(singleCellRNASeqMouseDeng2014)
library(RColorBrewer)
library(data.table)
library(tibble)

#############################################################################################
##GTEx V7 data##
gtex.counts <- read.table("/Users/poom/Desktop/countcluster/gtex_gene_count_brain_final.txt", header = TRUE, sep = "\t")
rownames(gtex.counts) = gtex.counts$gene_id
gtex.counts <- gtex.counts %>% select(-gene_id)
colnames(gtex.counts) <- str_replace_all(colnames(gtex.counts), "[.]", "-")
gtex.meta_data <- read.table("/Users/poom/Desktop/countcluster/gtex_meta_data_brain_final.txt", header = FALSE, sep = "\t")

##Start clustering##
a=FitGoM(t(gtex.counts),K=6,tol=1);
omega=a$fit$omega;
tissue_labels_brain = gtex.meta_data$V3
docweights_per_tissue_mean <- apply(omega, 2,
                                    function(x) { tapply(x, tissue_labels_brain, mean) })


##Draw structure plot##
annotation <- data.frame(
  sample_id = paste0("X", 1:length(tissue_labels_brain)),
  tissue_label = factor(tissue_labels_brain,
                        levels = rev(unique(tissue_labels_brain) ) ) );

tiff("/Users/poom/Desktop/countclust_brain_splot.tiff", units="in",  width=11, height=8, res=300)
StructureGGplot(omega = omega,
                annotation= annotation,
                yaxis_label = "",
                palette = RColorBrewer::brewer.pal(6, "Dark2"),
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
tiff("/Users/poom/Desktop/countclust_brain_hmap.tiff", units="in",  width=13, height=8, res=300)
ordering <- heatmap(docweights_per_tissue_mean, distfun = function(x) dist(x, method = "manhattan"))$rowInd
dev.off()

dist = as.matrix(dist(docweights_per_tissue_mean, diag = TRUE, upper = TRUE, method = "euclidean"))
result = data.frame()
for (i in 1:nrow(dist)) {
  temp = t(as.matrix(colnames(dist)[order(dist[i,], decreasing = FALSE)]))
  result = rbind(result, temp)
}

result$summ = paste(result$V1, result$V2, result$V3, result$V4, result$V5, result$V6, result$V7, result$V8,
                    result$V9, result$V10, result$V11, result$V12, result$V13, sep = " > ")
result = result %>% select(V1, summ)
write.table(result, "/Users/poom/Desktop/countcluster/result_brain_final.txt", quote= FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
