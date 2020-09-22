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
#pheno_list = c("anx_cc", "anx_fs", "covid", "ds", "ext", "fa_bmd")

#method_list = c("multi-omics", "predixcan", "utmost", "dapg", "mashr", "combined", "combined_mashr")
method_list = c("dapg_mashr", "dapg_utmost", "predixcan_mashr", "predixcan_utmost")

ref = read.table("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/id_conversion/ensembl_entrez_id.txt", sep = "\t", header = TRUE)
all_count = read.table("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/iv_all_count.txt", sep = "\t", header = TRUE)
all_count = apply(all_count[,2:4], 1, max)

##Import result##
result = list()

for (idx in 1:length(method_list)){
	method = method_list[idx]

	if ( method != "combined" && method != "combined_mashr" && method != "dapg_mashr" && method != "dapg_utmost" && method != "predixcan_mashr" && method != "predixcan_utmost" ){
		result_temp = data.frame()
		for (pheno in pheno_list){
			print(pheno)

			count = 0
			for (tissue in tissue_list){
				print(tissue)		
				path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/twas_assoc_new/", method, "/", pheno, "/", tissue, "/output", sep = "")
				filename = dir(path)
			
				df = data.frame()
				for (i in 1:length(filename)){
					temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE) %>% select(ID, TWAS.P)
					df = rbind(df, temp)
				}

				##Subset for protein coding genes##
				#df = df[df$ID %in% ref$ensembl_gene_id, ]
				df = df %>% filter(TWAS.P < 0.05/all_count[which(tissue_list == tissue)])
				count = count + nrow(df)
			}
			temp = cbind("pheno" = pheno, "total_count"= count)   
			result_temp = rbind(result_temp, temp)
		}
		result[[idx]] = result_temp
	} else {
		result_temp = data.frame()
                for (pheno in pheno_list){
                        print(pheno)

                        count1 = 0
			count2 = 0
                        for (tissue in tissue_list){
                                print(tissue)
                                path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/twas_assoc_new/", method, "/", pheno, "/", tissue, "/output", sep = "")
                                filename = dir(path)

                                df = data.frame()
                                for (i in 1:length(filename)){
                                        temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE) %>% select(gene.vec, pval.minp, pval.cauchy)
                                        df = rbind(df, temp)
                                }

                                ##Subset for protein coding genes##
                                #df = df[df$gene.vec %in% ref$ensembl_gene_id, ]
				threshold = all_count[which(tissue_list == tissue)]
                                df1 = df %>% filter(pval.minp < 0.05/threshold)
				df2 = df %>% filter(pval.cauchy < 0.05/threshold)
                                count1 = count1 + nrow(df1)
				count2 = count2 + nrow(df2)
                        }
                        temp = cbind("pheno" = pheno, "total_count_minp"= count1, "total_count_cauchy" = count2)
                        result_temp = rbind(result_temp, temp)
                }
                result[[idx]] = result_temp
	}

}

write.table(result[[1]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[1], "_gta_count.txt", sep =""), 
	    col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(result[[2]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[2], "_gta_count.txt", sep =""),  
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(result[[3]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[3], "_gta_count.txt", sep =""),  
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(result[[4]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[4], "_gta_count.txt", sep =""),
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
#write.table(result[[5]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[5], "_gta_count.txt", sep =""),
#            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
#write.table(result[[6]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[6], "_gta_count.txt", sep =""),
#            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
#write.table(result[[7]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[7], "_gta_count.txt", sep =""),
#            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
