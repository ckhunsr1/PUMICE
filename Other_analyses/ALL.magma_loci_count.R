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

#method_list = c("multi-omics", "predixcan", "utmost", "dapg", "mashr", "combined", "combined_mashr")
method_list = c("dapg_utmost", "predixcan_mashr", "predixcan_utmost")

ref = read.table("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/id_conversion/ensembl_entrez_id.txt", sep = "\t", header = TRUE)
ref = ref[!is.na(ref$entrez_id), ]
colnames(ref)[1] = "ensembl_gene_id"

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

			##Import MAGMA result##
			path = paste("/gpfs/group/dxl46/default/private/poom/magma/magma_output/", pheno, "_ss.35kbup_10_down.genes.out", sep = "")
			magma = read.table(path, header = TRUE)

			##Analyze magma result##
			magma = merge(magma, ref, by.x = "GENE", by.y = "entrez_id")		
			magma = magma[!duplicated(magma$GENE),]
			magma = magma %>% select(ensembl_gene_id, P)		
			magma = magma %>% filter(P < 0.05/nrow(magma))		

			sig_gene = c()
			for (tissue in tissue_list){
				print(tissue)		
				path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/twas_assoc_new/", method, "/", pheno, "/", tissue, "/output", sep = "")
				filename = dir(path)
			
				df = data.frame()
				for (i in 1:length(filename)){
					temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE) %>% select(ID, TWAS.P)
					df = rbind(df, temp)
				}

				df = df %>% filter(TWAS.P < 0.05/all_count[which(tissue_list == tissue)])
				sig_gene_temp = as.character(df$ID)
				sig_gene = unique(c(sig_gene, sig_gene_temp))
			}
			uni_gene = setdiff(sig_gene, as.character(magma$ensembl_gene_id))
			temp = cbind("pheno" = pheno, "total_count" = length(sig_gene), "novel_count" = length(uni_gene))   
			result_temp = rbind(result_temp, temp)
		}
		result[[idx]] = result_temp
	} else {
		result_temp = data.frame()
                for (pheno in pheno_list){
                        print(pheno)

                        ##Import MAGMA result##
                        path = paste("/gpfs/group/dxl46/default/private/poom/magma/magma_output/", pheno, "_ss.35kbup_10_down.genes.out", sep = "")
                        magma = read.table(path, header = TRUE)

                        ##Analyze magma result##
                        magma = merge(magma, ref, by.x = "GENE", by.y = "entrez_id")
                        magma = magma[!duplicated(magma$GENE),]
                        magma = magma %>% select(ensembl_gene_id, P)
                        magma = magma %>% filter(P < 0.05/nrow(magma))

                        sig_gene1 = c()
			sig_gene2 = c()
                        for (tissue in tissue_list){
                                print(tissue)
                                path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/twas_assoc_new/", method, "/", pheno, "/", tissue, "/output", sep = "")
                                filename = dir(path)

                                df = data.frame()
                                for (i in 1:length(filename)){
                                        temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE) %>% select(gene.vec, pval.minp, pval.cauchy)
                                        df = rbind(df, temp)
                                }

				threshold = all_count[which(tissue_list == tissue)]
				df1 = df %>% filter(pval.minp < 0.05/threshold)
                                sig_gene1_temp = as.character(df1$gene.vec)
                                sig_gene1 = unique(c(sig_gene1, sig_gene1_temp))
				
				df2 = df %>% filter(pval.cauchy < 0.05/threshold)
                                sig_gene2_temp = as.character(df2$gene.vec)
                                sig_gene2 = unique(c(sig_gene2, sig_gene2_temp))

                        }
                        uni_gene1 = setdiff(sig_gene1, as.character(magma$ensembl_gene_id))
                        uni_gene2 = setdiff(sig_gene2, as.character(magma$ensembl_gene_id))
			temp = cbind("pheno" = pheno, "total_count_minp" = length(sig_gene1), "novel_count_minp" = length(uni_gene1),
				     "total_count_cauchy" = length(sig_gene2), "novel_count_cauchy" = length(uni_gene2))
                        result_temp = rbind(result_temp, temp)
                }
                result[[idx]] = result_temp
	}
}

write.table(result[[1]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[1], "_magma_loci_count.txt", sep =""), 
	    col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(result[[2]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[2], "_magma_loci_count.txt", sep =""),  
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(result[[3]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[3], "_magma_loci_count.txt", sep =""),  
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
#write.table(result[[4]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[4], "_magma_loci_count.txt", sep =""),
#            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
#write.table(result[[5]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[5], "_magma_loci_count.txt", sep =""),
#            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
#write.table(result[[6]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[6], "_magma_loci_count.txt", sep =""),
#            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
#write.table(result[[7]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[7], "_magma_loci_count.txt", sep =""),
#            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
