library(dplyr)
threshold = 0.0005

#method_list = c("multi-omics", "predixcan", "utmost", "dapg", "mashr", "combined", "combined_mashr", "dapg_mashr", "dapg_utmost")
method_list = c("combined", "combined_mashr", "dapg_mashr", "dapg_utmost", "predixcan_mashr", "predixcan_utmost")

tissue_list = c("Muscle_Skeletal", "Skin_Sun_Exposed_Lower_leg", "Thyroid", "Adipose_Subcutaneous", "Lung", "Artery_Tibial",
                        "Whole_Blood", "Nerve_Tibial", "Esophagus_Mucosa", "Skin_Not_Sun_Exposed_Suprapubic", "Esophagus_Muscularis", "Adipose_Visceral_Omentum",
                        "Cells_Transformed_fibroblasts", "Artery_Aorta", "Heart_Left_Ventricle", "Heart_Atrial_Appendage", "Breast_Mammary_Tissue", "Colon_Transverse",
                        "Stomach", "Testis", "Esophagus_Gastroesophageal_Junction", "Colon_Sigmoid", "Pancreas", "Pituitary",
                        "Adrenal_Gland", "Brain_Cerebellum", "Liver", "Brain_Caudate_basal_ganglia", "Artery_Coronary", "Brain_Cortex",
                        "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Spleen", "Prostate", "Brain_Frontal_Cortex_BA9", "Brain_Anterior_cingulate_cortex_BA24",
                        "Small_Intestine_Terminal_Ileum", "Brain_Hippocampus", "Brain_Putamen_basal_ganglia", "Brain_Hypothalamus", "Ovary", "Cells_EBV-transformed_lymphocytes",
                        "Vagina", "Uterus", "Brain_Amygdala", "Brain_Spinal_cord_cervical_c-1", "Minor_Salivary_Gland", "Brain_Substantia_nigra")

#pheno_list = c("ast", "cd", "ecz", "ibd", "ra", "sle", "uc", "vit")
pheno_list = c("covid")

all_count = read.table("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/iv_all_count.txt", sep = "\t", header = TRUE)
all_count = apply(all_count[,2:4], 1, max)

##Convert ensembl gene id to entrez id##
conv = read.table("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/id_conversion/ensembl_entrez_id.txt",
              header = TRUE, sep= "\t", na.strings = "")
conv = conv[!(is.na(conv$entrez_id)), ]

##Make sure that it is one-to-one mapping##
dup1 = conv[duplicated(conv[,2]), 2]
conv = conv[!(conv$entrez_id %in% dup1), ]
dup2 = conv[duplicated(conv[,1]), 1]
conv = conv[!(conv$ensembl_id %in% dup2), ]

for (pheno in pheno_list) {

	for (method in method_list){
	if ( method != "combined" && method != "combined_mashr" && method != "dapg_mashr" && method != "dapg_utmost" && method != "predixcan_mashr" && method != "predixcan_utmost" ){
		result = data.frame()
		for (tissue in tissue_list){

			df = data.frame()
			for (i in 1:22){
				path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/twas_assoc_new/", method, "/", pheno, "/", tissue, "/output/", 
						"chr", i, "_result.txt", sep = "")
				temp = read.table(path, header = TRUE, fill=TRUE)
				df = rbind(df, temp)
			}
			df = df[!(is.na(df$TWAS.Z)),]
			#threshold = 0.05/all_count[which(tissue_list == tissue)]
			df = df %>% filter(TWAS.P < threshold) %>% select(ID, TWAS.Z)
			result = rbind(result, df)
		}
		result = as.data.frame(result %>% group_by(ID) %>% dplyr::summarize(Mean = mean(TWAS.Z, na.rm=TRUE)))

        	merge = merge(result, conv, by.x = "ID", by.y = "ensembl_id")

        	up = merge %>% filter(Mean > 0)
        	up = up[order(-up$Mean),] %>% select(entrez_id)
		up = as.data.frame(cbind(paste(pheno, method, sep = "_"), "up", t(up)))

        	down = merge %>% filter(Mean < 0)
        	down = down[order(down$Mean),] %>% select(entrez_id)
		down = as.data.frame(cbind(paste(pheno, method, sep = "_"), "down", t(down)))

        	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/clue_sig_new/", method, "/", pheno, "_", method, "_clue_up_input_new.txt", sep = "")
        	write.table(up, path, col.names =  FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

        	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/clue_sig_new/", method, "/", pheno, "_", method, "_clue_down_input_new.txt", sep = "")
        	write.table(down, path, col.names =  FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

	} else {
		result1 = data.frame()
                result2 = data.frame()
                for (tissue in tissue_list){

                        df = data.frame()
                        for (i in 1:22){
                                path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/twas_assoc_new/", method, "/", pheno, "/", tissue, "/output/",
                                                "chr", i, "_result.txt", sep = "")
                                temp = read.table(path, header = TRUE, fill=TRUE)
                                df = rbind(df, temp)
                        }

                        #threshold = 0.05/all_count[which(tissue_list == tissue)]
                        df1 = df[!(is.na(df$twas.z.minp)),]
                        df1 = df1 %>% filter(pval.minp < threshold) %>% select(gene.vec, twas.z.minp)

                        df2 = df[!(is.na(df$twas.z.cauchy)),]
                        df2 = df2 %>% filter(pval.cauchy < threshold) %>% select(gene.vec, twas.z.cauchy)

                        result1 = rbind(result1, df1)
                        result2 = rbind(result2, df2)
                }

		result2 = as.data.frame(result2 %>% group_by(gene.vec) %>% dplyr::summarize(Mean = mean(twas.z.cauchy, na.rm=TRUE)))
                merge2 = merge(result2, conv, by.x = "gene.vec", by.y = "ensembl_id")
                up2 = merge2 %>% filter(Mean > 0)
                up2 = up2[order(-up2$Mean),] %>% select(entrez_id)
                up2 = as.data.frame(cbind(paste(pheno, "cauchy", method, sep = "_"), "up", t(up2)))
                down2 = merge2 %>% filter(Mean < 0)
                down2 = down2[order(down2$Mean),] %>% select(entrez_id)
                down2 = as.data.frame(cbind(paste(pheno, "cauchy", method, sep = "_"), "down", t(down2)))
                path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/clue_sig_new/cauchy_", method, "/", pheno, "_cauchy_", method, "_clue_up_input_new.txt", sep = "")
                write.table(up2, path, col.names =  FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
                path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/clue_sig_new/cauchy_", method, "/", pheno, "_cauchy_", method, "_clue_down_input_new.txt", sep = "")
                write.table(down2, path, col.names =  FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

	}
	}
}
