args = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(stringr)

tissue_list = c("Muscle_Skeletal", "Skin_Sun_Exposed_Lower_leg", "Thyroid", "Adipose_Subcutaneous", "Lung", "Artery_Tibial",
                        "Whole_Blood", "Nerve_Tibial", "Esophagus_Mucosa", "Skin_Not_Sun_Exposed_Suprapubic", "Esophagus_Muscularis", "Adipose_Visceral_Omentum",
                        "Cells_Transformed_fibroblasts", "Artery_Aorta", "Heart_Left_Ventricle", "Heart_Atrial_Appendage", "Breast_Mammary_Tissue", "Colon_Transverse",
                        "Stomach", "Testis", "Esophagus_Gastroesophageal_Junction", "Colon_Sigmoid", "Pancreas", "Pituitary",
                        "Adrenal_Gland", "Brain_Cerebellum", "Liver", "Brain_Caudate_basal_ganglia", "Artery_Coronary", "Brain_Cortex",
                        "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Spleen", "Prostate", "Brain_Frontal_Cortex_BA9", "Brain_Anterior_cingulate_cortex_BA24",
                        "Small_Intestine_Terminal_Ileum", "Brain_Hippocampus", "Brain_Putamen_basal_ganglia", "Brain_Hypothalamus", "Ovary", "Cells_EBV-transformed_lymphocytes",
                        "Vagina", "Uterus", "Brain_Amygdala", "Brain_Spinal_cord_cervical_c-1", "Minor_Salivary_Gland", "Brain_Substantia_nigra")

pheno_list = setdiff(dir("/gpfs/group/dxl46/default/private/poom/espresso/sumstats_clean_new"), "check.R")

method_list = args[1]

##Import p-values##
result = list()

if (method_list[1] != "combined" && method_list[1] != "combined_mashr" && method_list[1] != "dapg_mashr" && method_list[1] != "dapg_utmost" && method_list[1] != "predixcan_mashr" && method_list[1] != "predixcan_utmost" ){
	for (idx in 1:length(method_list)){
		method = method_list[idx]
		result_temp = data.frame()

		for (pheno in pheno_list){
			print(pheno)
			for (tissue in tissue_list){
				print(tissue)		
				path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/twas_assoc_new/", method, "/", pheno, "/", tissue, "/output", sep = "")
				filename = dir(path)
			
				df = data.frame()
				for (i in 1:length(filename)){
					temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE) %>% select(ID, TWAS.P)
					temp = cbind("pheno" = pheno, "tisue" = tissue, temp)
					df = rbind(df, temp)
				}
				result_temp = rbind(result_temp, df)				
			}
			result[[idx]] = result_temp   
		}
	}
} else {
	for (idx in 1:length(method_list)){
		method = method_list[idx]
		result_temp = data.frame()

		for (pheno in pheno_list){
                        print(pheno)
                        for (tissue in tissue_list){
                                print(tissue)
                                path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/twas_assoc_new/", method, "/", pheno, "/", tissue, "/output", sep = "")
                                filename = dir(path)

                                df = data.frame()
                                for (i in 1:length(filename)){
					print(i)
                                        temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE) %>% select(gene.vec, pval.minp, pval.cauchy)
					temp = cbind("pheno" = pheno, "tisue" = tissue, temp)
                                        df = rbind(df, temp)
                                }
                                result_temp = rbind(result_temp, df)
                        }
                        result[[idx]] = result_temp
                }
        }
}
write.table(result[[1]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[1], "_QQ_analysis.txt", sep =""), 
	    col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

