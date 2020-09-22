args = commandArgs(trailingOnly=TRUE)

##Function to import internal validation##
import_iv <- function(path, filename){
        iv = data.frame()
        for (i in 1:length(filename)){
                print(i)
                try({
                temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE)
                iv = rbind(iv, temp)
                })
        }
        iv 
}

##Function to import external prediction in significant genes##
import_ext <- function(path, filename){
	ext_temp = data.frame()
	for (i in 1:length(filename)){
	        print(i)
	        try({
		temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE)
	        ext_temp = rbind(ext_temp, temp)
	        })
	}
	ext_temp
}

##Function to import coef of all genes##
import_coef <- function(path, filename){
        coef_temp = data.frame()
        for (i in 1:length(filename)){
                print(i)
                try({
                temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE)
                coef_temp = rbind(coef_temp, temp)
                })
        }
	coef_temp
}


######################################################################################################################################################################################################
library(dplyr)
library(tidyverse)
tissue_list = dir("/gpfs/group/dxl46/default/private/poom/mashr/fastqtl_input/phenotype")

######################################################################################################################################################################################################
count = data.frame()
for (tissue in tissue_list){

	##Import internal validation result for multi-omics and predixcan##
	type = "multi-omics"
	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/internal_validation/", tissue, sep = "")
	filename = dir(path)
	multi = import_iv(path, filename)
	multi_sig = multi %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05)

	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/output/", tissue, sep = "")
        filename = dir(path)
	multi_coef = import_coef(path, filename)
	multi_sig = multi_sig[multi_sig$gene_id %in% unique(multi_coef$gene_id), ]

	##Import internal validation result for multi-omics and predixcan##
	type = "predixcan"
	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/internal_validation/", tissue, sep = "")
	filename = dir(path)
	predixcan = import_iv(path, filename)
	predixcan_sig = predixcan %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05)

	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/output/", tissue, sep = "")
        filename = dir(path)
        predixcan_coef = import_coef(path, filename)
        predixcan_sig = predixcan_sig[predixcan_sig$gene_id %in% unique(predixcan_coef$gene_id), ]

	##Import internal validation result for utmost##
	type = "utmost"
	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/internal_validation/", tissue, sep = "")
	filename = dir(path)
	utmost = import_iv(path, filename)
	utmost_sig = utmost %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05)
	
	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/output/", tissue, sep = "")
        filename = dir(path)
        utmost_coef = import_coef(path, filename)
        utmost_sig = utmost_sig[utmost_sig$gene_id %in% unique(utmost_coef$gene_id), ]

	##Import internal validation result for epixcan##
	type = "epixcan"
	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/internal_validation/", tissue, sep = "")
	filename = dir(path)
	epixcan = import_iv(path, filename)
	epixcan_sig = epixcan %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05)
	
	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/output/", tissue, sep = "")
        filename = dir(path)
        epixcan_coef = import_coef(path, filename)
        epixcan_sig = epixcan_sig[epixcan_sig$gene_id %in% unique(epixcan_coef$gene_id), ]

	##Import internal validation result for dpr_add##
	type = "dpr_add"
	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/internal_validation/", tissue, sep = "")
	filename = dir(path)
	dpr_add = import_iv(path, filename)
	dpr_add_sig = dpr_add %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05)

	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/output/", tissue, sep = "")
        filename = dir(path)
        dpr_add_coef = import_coef(path, filename)
        dpr_add_sig = dpr_add_sig[dpr_add_sig$gene_id %in% unique(dpr_add_coef$gene_id), ]

	##Import internal validation result for blup##
	type = "blup"
	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/internal_validation/", tissue, sep = "")
	filename = dir(path)
	blup = import_iv(path, filename)
	blup_sig = blup %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05)
	
	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/output/", tissue, sep = "")
        filename = dir(path)
        blup_coef = import_coef(path, filename)
        blup_sig = blup_sig[blup_sig$gene_id %in% unique(blup_coef$gene_id), ]

	##Import internal validation result for bslmm##
	type = "bslmm"
	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/internal_validation/", tissue, sep = "")
	filename = dir(path)
	bslmm = import_iv(path, filename)
	bslmm_sig = bslmm %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05)

	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/output/", tissue, sep = "")
        filename = dir(path)
        bslmm_coef = import_coef(path, filename)
        bslmm_sig = bslmm_sig[bslmm_sig$gene_id %in% unique(bslmm_coef$gene_id), ]

	#Import internal validation result for mcmc_add##
	type = "mcmc_add"
	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/internal_validation/", tissue, sep = "")
	filename = dir(path)
	mcmc_add = import_iv(path, filename)
	mcmc_add_sig = mcmc_add %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05)

	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/output/", tissue, sep = "")
        filename = dir(path)
        mcmc_add_coef = import_coef(path, filename)
        mcmc_add_sig = mcmc_add_sig[mcmc_add_sig$gene_id %in% unique(mcmc_add_coef$gene_id), ]

	#Import internal validation result for dapg##
	type = "dapg"
	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/internal_validation/", tissue, sep = "")
	filename = dir(path)
	dapg = import_iv(path, filename)
	dapg_sig = dapg %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05)

	##Only for dapg, we also need to filter out genes with no variants with PIP > 0.01##
	dapg_pip = c()
	for (i in 1:22){
		load(paste("/gpfs/group/dxl46/default/private/poom/dapg/prior/", tissue, "/", "prior_chr", i, ".RData", sep = ""))
		final_df = final_df %>% filter(pip > 0.01)
		temp = unique(final_df$gene)
		dapg_pip = c(dapg_pip, temp)
	}
	dapg_sig = dapg_sig[dapg_sig$gene_id %in% dapg_pip, ]

	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/output/", tissue, sep = "")
        filename = dir(path)
        dapg_coef = import_coef(path, filename)
        dapg_sig = dapg_sig[dapg_sig$gene_id %in% unique(dapg_coef$gene_id), ]

	##Import mashr sig gene##
	type = "mashr"
	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/sig/", tissue, "/", tissue, "_allchr_GTEx_mashr_sig.txt", sep = "")
	mashr_sig = read.table(path, header = TRUE)
	
	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/output/", tissue, sep = "")
        filename = dir(path)
        mashr_coef = import_coef(path, filename)
        mashr_sig = mashr_sig[mashr_sig$sig_gene %in% unique(mashr_coef$gene_id), , drop = FALSE]

	count_temp = cbind(tissue, nrow(multi_sig), nrow(predixcan_sig), nrow(utmost_sig), nrow(epixcan_sig), 
			   nrow(dpr_add_sig), nrow(blup_sig), nrow(bslmm_sig), nrow(mcmc_add_sig), nrow(dapg_sig), nrow(mashr_sig))
	count = rbind(count, count_temp)
}
colnames(count) = c("tissue", "multi-omics", "predixcan", "utmost", "epixcan", "dpr_add", "blup", "bslmm", "mcmc_add", "dapg", "mashr")

write.table(count, "/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/sig_model_count.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
######################################################################################################################################################################################################
