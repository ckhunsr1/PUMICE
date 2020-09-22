library(dplyr)
library(tidyverse)

type = "multi-omics"
process = "coef_total"

tissue_list = c("Muscle_Skeletal", "Skin_Sun_Exposed_Lower_leg", "Thyroid", "Adipose_Subcutaneous", "Lung", "Artery_Tibial",
                        "Whole_Blood", "Nerve_Tibial", "Esophagus_Mucosa", "Skin_Not_Sun_Exposed_Suprapubic", "Esophagus_Muscularis", "Adipose_Visceral_Omentum",
                        "Cells_Transformed_fibroblasts", "Artery_Aorta", "Heart_Left_Ventricle", "Heart_Atrial_Appendage", "Breast_Mammary_Tissue", "Colon_Transverse",
                        "Stomach", "Testis", "Esophagus_Gastroesophageal_Junction", "Colon_Sigmoid", "Pancreas", "Pituitary",
                        "Adrenal_Gland", "Brain_Cerebellum", "Liver", "Brain_Caudate_basal_ganglia", "Artery_Coronary", "Brain_Cortex",
                        "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Spleen", "Prostate", "Brain_Frontal_Cortex_BA9", "Brain_Anterior_cingulate_cortex_BA24",
                        "Small_Intestine_Terminal_Ileum", "Brain_Hippocampus", "Brain_Putamen_basal_ganglia", "Brain_Hypothalamus", "Ovary", "Cells_EBV-transformed_lymphocytes",
                        "Vagina", "Uterus", "Brain_Amygdala", "Brain_Spinal_cord_cervical_c-1", "Minor_Salivary_Gland", "Brain_Substantia_nigra")

window = as.data.frame(matrix(0, length(tissue_list), 6))
penalty = as.data.frame(matrix(0, length(tissue_list), 7)) 
for (idx in 1:length(tissue_list)){
	tissue = tissue_list[idx]

	##Import internal validation result##
	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/", process, "/", type, "/internal_validation/", tissue, sep = "")
	filename = dir(path, pattern = "chr")

	sig = data.frame()
	for (i in 1:length(filename)){
        	print(i)
        	try({
		temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE) %>%
                	filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05)
        	sig = rbind(sig, temp)
        	})
	}
	window[idx,] = t(sig %>% group_by(type) %>% tally)[2,]
	penalty[idx,] = t(sig %>% group_by(penalty) %>% tally)[2,]
}
colnames(window) = c("1000", "250", "Domain", "Loop", "pcHiC", "TAD")
colnames(penalty) = 1:7
window$tissue = tissue_list
penalty$tissue = tissue_list

index = c(7, 1:6)
window = window[, index]
index = c(8, 1:7)
penalty = penalty[, index]

write.table(window, paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", type, "_window_composition.txt", sep =""),
               col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(penalty, paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", type, "_penalty_composition.txt", sep =""),
               col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
