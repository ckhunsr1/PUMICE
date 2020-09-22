args = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(tidyverse)

type = args[1]
process = "coef_total"
cut_off = 1000000 ##Window size to look at##
n = 50 ##Numbers of bin##

tissue_list = c("Muscle_Skeletal", "Skin_Sun_Exposed_Lower_leg", "Thyroid", "Adipose_Subcutaneous", "Lung", "Artery_Tibial",
                        "Whole_Blood", "Nerve_Tibial", "Esophagus_Mucosa", "Skin_Not_Sun_Exposed_Suprapubic", "Esophagus_Muscularis", "Adipose_Visceral_Omentum",
                        "Cells_Transformed_fibroblasts", "Artery_Aorta", "Heart_Left_Ventricle", "Heart_Atrial_Appendage", "Breast_Mammary_Tissue", "Colon_Transverse",
                        "Stomach", "Testis", "Esophagus_Gastroesophageal_Junction", "Colon_Sigmoid", "Pancreas", "Pituitary",
                        "Adrenal_Gland", "Brain_Cerebellum", "Liver", "Brain_Caudate_basal_ganglia", "Artery_Coronary", "Brain_Cortex",
                        "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Spleen", "Prostate", "Brain_Frontal_Cortex_BA9", "Brain_Anterior_cingulate_cortex_BA24",
                        "Small_Intestine_Terminal_Ileum", "Brain_Hippocampus", "Brain_Putamen_basal_ganglia", "Brain_Hypothalamus", "Ovary", "Cells_EBV-transformed_lymphocytes",
                        "Vagina", "Uterus", "Brain_Amygdala", "Brain_Spinal_cord_cervical_c-1", "Minor_Salivary_Gland", "Brain_Substantia_nigra")
#tissue_list = c("Whole_Blood")

result = data.frame()
count = data.frame()
for (tissue in tissue_list){
###########################################################################################################################
	##Import coef##
	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/", process, "/", type, "/output/", tissue, sep = "")
	filename = dir(path)
	filename = unique(str_remove(filename, "_new"))

	setwd(path)
	coef = data.frame()
	for (i in 1:length(filename)) {
		tryCatch({
		print(i)
		temp = read.table(filename[i], header = TRUE, na.string ='.', as.is=TRUE, check.names=FALSE)
		coef = rbind(coef, temp)
		}, error=function(e){})
	}
	colnames(coef)[2] = "gene_id" 

###########################################################################################################################
	##Import internal validation result##
	##Filter in only significant genes##
	##Import internal validation result##

	if ( type != "mashr" ) {
		path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/", process, "/", type, "/internal_validation/", tissue, sep = "")
	        filename = dir(path, pattern = "chr")
		sig = c()
		for (i in 1:length(filename)){
        		print(i)
        		try({
			temp = as.character((read.table(paste(path, "/", filename[i], sep = ""), header = TRUE) %>%
        	        	filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05) %>% select(gene_id))[,1])
			sig = c(sig, temp)
        		})
		}
	} else {
		path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/", process, "/", type, "/sig/", tissue, sep = "")
		filename = dir(path, pattern = "chr")
		sig = as.character((read.table(paste(path, "/", filename, sep = ""), header = TRUE))[,1])		
	}

	coef = coef[coef$gene_id %in% sig, ]
	count_temp = coef %>% group_by(gene_id) %>% tally %>% select(n)
	count = rbind(count, count_temp)
###########################################################################################################################
	##Import gene position##
	gene = read.table("/gpfs/group/dxl46/default/private/poom/GTEx/input/ensemble_to_genomic_chr1-23_total.txt", header = TRUE, na.string ='.', as.is=TRUE, check.names=FALSE)

	merge = merge(x = coef, y = gene, by.x = "gene_id", by.y = "ensembl_gene_id") %>% select(snp, start, end)
	merge$snp_pos = sapply(strsplit(merge$snp, "_"), function(x){as.character(x[2])})
	merge = merge %>% select(snp_pos, start, end)
	colnames(merge) = c("snp_pos", "gene_start_pos", "gene_end_pos")
	merge$snp_pos = as.numeric(as.character(merge$snp_pos))
	merge$gene_start_pos = as.numeric(as.character(merge$gene_start_pos))
	merge$gene_end_pos = as.numeric(as.character(merge$gene_end_pos))
	merge$gene_length = merge$gene_end_pos - merge$gene_start_pos + 1

	##Find SNPs within gene regions##
	intra_snp = merge[which(merge$snp_pos >= merge$gene_start_pos & merge$snp_pos <= merge$gene_end_pos), ]
	inter_snp_bf = merge[which(merge$snp_pos < merge$gene_start_pos),]
	inter_snp_af = merge[which(merge$snp_pos > merge$gene_end_pos),]

	##Find relative position##
	intra_snp$rel_pos = (intra_snp$snp_pos - intra_snp$gene_start_pos)*(cut_off)/(intra_snp$gene_length)
	inter_snp_bf$rel_pos = inter_snp_bf$gene_start_pos - inter_snp_bf$snp_pos
	inter_snp_af$rel_pos = inter_snp_af$snp_pos - inter_snp_af$gene_end_pos

	##Filter relative position##
	inter_snp_bf = inter_snp_bf %>% filter(rel_pos <= cut_off)
	inter_snp_af = inter_snp_af %>%	filter(rel_pos <= cut_off)

	##Bin relative position##
	b = seq(0, cut_off, length = n+1)
	bin_name = c(1:n)

	inter_snp_bf$rel_bin = as.numeric(cut(-inter_snp_bf$rel_pos, breaks = -b, labels = bin_name))
	intra_snp$rel_bin = as.numeric(cut(intra_snp$rel_pos, breaks = b, labels = bin_name)) + n
	inter_snp_af$rel_bin = as.numeric(cut(inter_snp_af$rel_pos, breaks = b, labels = bin_name)) + (2*n)

	##Final output##
	df = rbind(inter_snp_bf, intra_snp, inter_snp_af) %>% select(rel_bin)

	result = rbind(result, df)

}

write.table(result, paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", type, "_rel_pos.txt", sep =""),
               col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(count, paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", type, "_snps_count.txt", sep =""),
               col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
