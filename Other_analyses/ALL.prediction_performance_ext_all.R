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

##Function to determine the performance of cross-validation cohort across all genes##
iv_comp <- function(multi, other, gene_list){
	multi_temp = multi[multi$gene_id %in% gene_list, ] %>% select(gene_id, pear_avg_t)
	if (nrow(multi_temp) < length(gene_list)) {
        	temp = setdiff(gene_list, multi_temp$gene_id)
        	temp = as.data.frame(cbind("gene_id" = temp, "pear_avg_t" = 0))
        	temp$pear_avg_t = as.numeric(as.character(temp$pear_avg_t))
        	multi_temp = rbind(multi_temp, temp)
	}

	other_temp = other[other$gene_id %in% gene_list, ] %>% select(gene_id, pear_avg_t)
	if (nrow(other_temp) < length(gene_list)) {
	        temp = setdiff(gene_list, other_temp$gene_id)
	        temp = as.data.frame(cbind("gene_id" = temp, "pear_avg_t" = 0))
	        temp$pear_avg_t = as.numeric(as.character(temp$pear_avg_t))
	        other_temp = rbind(other_temp, temp)
	}

	merge = merge(multi_temp, other_temp, by.x = "gene_id", by.y = "gene_id") %>% select(gene_id, pear_avg_t.x, pear_avg_t.y)
	merge[merge$pear_avg_t.x < 0, "pear_avg_t.x"] = 0
	merge[merge$pear_avg_t.y < 0, "pear_avg_t.y"] = 0
	merge$delta = merge$pear_avg_t.x - merge$pear_avg_t.y
	merge
}


##Function to determine the performance of independent cohort across all significant genes##
ext_comp <- function(multi_ext, other_ext, gene_list){
        multi_temp = multi_ext[multi_ext$gene_id %in% gene_list, ] %>% select(gene_id, r)
        if (nrow(multi_temp) < length(gene_list)) {
                temp = setdiff(gene_list, multi_temp$gene_id)
                temp = as.data.frame(cbind("gene_id" = temp, "r" = 0))
                temp$r = as.numeric(as.character(temp$r))
                multi_temp = rbind(multi_temp, temp)
        }

	other_temp = other_ext[other_ext$gene_id %in% gene_list, ] %>% select(gene_id, r)
        if (nrow(other_temp) < length(gene_list)) {
                temp = setdiff(gene_list, other_temp$gene_id)
                temp = as.data.frame(cbind("gene_id" = temp, "r" = 0))
                temp$r = as.numeric(as.character(temp$r))
                other_temp = rbind(other_temp, temp)
        }

	merge = merge(multi_temp, other_temp, by.x = "gene_id", by.y = "gene_id") %>% select(gene_id, r.x, r.y)
        merge[merge$r.x < 0, "r.x"] = 0
        merge[merge$r.y < 0, "r.y"] = 0
        merge$delta = merge$r.x - merge$r.y
        merge
}

######################################################################################################################################################################################################
options(stringsAsFactors=F)
library(dplyr)
library(tidyverse)
tissue = args[1]
ext = args[2]

##Import testing datasets##
expression_path = paste("/gpfs/group/dxl46/default/private/poom/", ext, "/input/expression_processed/", ext, "_normalized_expression_final.txt", sep = "")
expression = read.table(expression_path, header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
external_gene = expression$gene_id

######################################################################################################################################################################################################
##Import internal validation result for multi-omics and predixcan##
type = "multi-omics"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/internal_validation/", tissue, sep = "")
filename = dir(path)
multi = import_iv(path, filename)

##Import internal validation result for multi-omics and predixcan##
type = "predixcan"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/internal_validation/", tissue, sep = "")
filename = dir(path)
predixcan = import_iv(path, filename)

##Import internal validation result for utmost##
type = "utmost"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/internal_validation/", tissue, sep = "")
filename = dir(path)
utmost = import_iv(path, filename)

##Import internal validation result for epixcan##
type = "epixcan"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/internal_validation/", tissue, sep = "")
filename = dir(path)
epixcan = import_iv(path, filename)

##Import internal validation result for dpr_add##
type = "dpr_add"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/internal_validation/", tissue, sep = "")
filename = dir(path)
dpr_add = import_iv(path, filename)

##Import internal validation result for blup##
type = "blup"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/internal_validation/", tissue, sep = "")
filename = dir(path)
blup = import_iv(path, filename)

##Import internal validation result for bslmm##
type = "bslmm"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/internal_validation/", tissue, sep = "")
filename = dir(path)
bslmm = import_iv(path, filename)

#Import internal validation result for mcmc_add##
type = "mcmc_add"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/internal_validation/", tissue, sep = "")
filename = dir(path)
mcmc_add = import_iv(path, filename)

#Import internal validation result for dapg##
type = "dapg"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/internal_validation/", tissue, sep = "")
filename = dir(path)
dapg = import_iv(path, filename)

##Import mashr sig gene##
type = "mashr"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/sig/", tissue, "/", tissue, "_allchr_GTEx_mashr_sig.txt", sep = "")
mashr = read.table(path, header = TRUE)

######################################################################################################################################################################################################
##Import external prediction for multi-omics##
type = "multi-omics"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/testing/", ext, "/", type, "/", tissue, sep = "")
filename = dir(path, pattern = tissue)
multi_sig = intersect(external_gene, (multi %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05))$gene_id)
multi_ext = import_ext(path, filename)

##Import external prediction for predixcan##
type = "predixcan"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/testing/", ext, "/", type, "/", tissue, sep = "")
filename = dir(path, pattern = tissue)
predixcan_sig = intersect(external_gene, (get(eval(type)) %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05))$gene_id)
predixcan_ext = import_ext(path, filename)

##Import external prediction for utmost##
type = "utmost"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/testing/", ext, "/", type, "/", tissue, sep = "")
filename = dir(path, pattern = tissue)
utmost_sig = intersect(external_gene, (get(eval(type)) %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05))$gene_id)
utmost_ext = import_ext(path, filename)

##Import external prediction for epixcan##
type = "epixcan"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/testing/", ext, "/", type, "/", tissue, sep = "")
filename = dir(path, pattern = tissue)
epixcan_sig = intersect(external_gene, (get(eval(type)) %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05))$gene_id)
epixcan_ext = import_ext(path, filename)

##Import external prediction for dpr_add##
type = "dpr_add"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/testing/", ext, "/", type, "/", tissue, sep = "")
filename = dir(path, pattern = tissue)
dpr_add_sig = intersect(external_gene, (get(eval(type)) %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05))$gene_id)
dpr_add_ext = import_ext(path, filename)

##Import external prediction for blup##
type = "blup"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/testing/", ext, "/", type, "/", tissue, sep = "")
filename = dir(path, pattern = tissue)
blup_sig = intersect(external_gene, (get(eval(type)) %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05))$gene_id)
blup_ext = import_ext(path, filename)

##Import external prediction for bslmm##
type = "bslmm"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/testing/", ext, "/", type, "/", tissue, sep = "")
filename = dir(path, pattern = tissue)
bslmm_sig = intersect(external_gene, (get(eval(type)) %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05))$gene_id)
bslmm_ext = import_ext(path, filename)

##Import external prediction for mcmc_add##
type = "mcmc_add"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/testing/", ext, "/", type, "/", tissue, sep = "")
filename = dir(path, pattern = tissue)
mcmc_add_sig = intersect(external_gene, (get(eval(type)) %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05))$gene_id)
mcmc_add_ext = import_ext(path, filename)

##Import external prediction for dapg##
type = "dapg"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/testing/", ext, "/", type, "/", tissue, sep = "")
filename = dir(path, pattern = tissue)
dapg_sig = intersect(external_gene, (get(eval(type)) %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05))$gene_id)

##Only for dapg, we also need to filter out genes with no variants with PIP > 0.01##
dapg_pip = c()
for (i in 1:22){
	load(paste("/gpfs/group/dxl46/default/private/poom/dapg/prior/", tissue, "/", "prior_chr", i, ".RData", sep = ""))
	final_df = final_df %>% filter(pip > 0.01)
	temp = unique(final_df$gene)
	dapg_pip = c(dapg_pip, temp)
}
dapg_sig = intersect(dapg_sig, dapg_pip)
dapg_ext = import_ext(path, filename)

##Import external prediction for mashr##
type = "mashr"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/testing/", ext, "/", type, "/", tissue, sep = "")
filename = dir(path, pattern = tissue)
mashr_sig = intersect(external_gene, mashr$sig_gene)
mashr_ext = import_ext(path, filename)

######################################################################################################################################################################################################
##Make comparison for independent cohorts##
ext_result = data.frame()

##Multi vs. PrediXcan##
#gene_list = union(multi_sig, predixcan_sig)
gene_list = intersect(multi_sig, predixcan_sig)
merge = ext_comp(multi_ext, predixcan_ext, gene_list)
merge$type = "mp"
ext_result = rbind(ext_result, merge)

##Multi vs. UTMOST##
#gene_list = union(multi_sig, utmost_sig)
gene_list = intersect(multi_sig, utmost_sig)
merge = ext_comp(multi_ext, utmost_ext, gene_list)
merge$type = "mu"
ext_result = rbind(ext_result, merge)

##Multi vs. EpiXcan##
#gene_list = union(multi_sig, epixcan_sig)
gene_list = intersect(multi_sig, epixcan_sig)
merge = ext_comp(multi_ext, epixcan_ext, gene_list)
merge$type = "me"
ext_result = rbind(ext_result, merge)

##Multi vs. DPR_add##
#gene_list = union(multi_sig, dpr_add_sig)
gene_list = intersect(multi_sig, dpr_add_sig)
merge = ext_comp(multi_ext, dpr_add_ext, gene_list)
merge$type = "md"
ext_result = rbind(ext_result, merge)

##Multi vs. BLUP##
#gene_list = union(multi_sig, blup_sig)
gene_list = intersect(multi_sig, blup_sig)
merge = ext_comp(multi_ext, blup_ext, gene_list)
merge$type = "mbl"
ext_result = rbind(ext_result, merge)

##Multi vs. BSLMM##
#gene_list = union(multi_sig, bslmm_sig)
gene_list = intersect(multi_sig, bslmm_sig)
merge = ext_comp(multi_ext, bslmm_ext, gene_list)
merge$type = "mbs"
ext_result = rbind(ext_result, merge)

##Multi vs. MCMC_add##
#gene_list = union(multi_sig, mcmc_add_sig)
gene_list = intersect(multi_sig, mcmc_add_sig)
merge = ext_comp(multi_ext, mcmc_add_ext, gene_list)
merge$type = "mmc"
ext_result = rbind(ext_result, merge)

##Multi vs. dapg##
#gene_list = union(multi_sig, dapg_sig)
gene_list = intersect(multi_sig, dapg_sig)
merge = ext_comp(multi_ext, dapg_ext, gene_list)
merge$type = "mdap"
ext_result = rbind(ext_result, merge)

##Multi vs. mashr##
#gene_list = union(multi_sig, mashr_sig)
gene_list = intersect(multi_sig, mashr_sig)
merge = ext_comp(multi_ext, mashr_ext, gene_list)
merge$type = "mmash"
ext_result = rbind(ext_result, merge)

##Write output##
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/summary/", tissue, "/", tissue, "_GTEx_", ext, "_result_intersect.txt", sep = "")
write.table(ext_result, path, col.names = TRUE, row.names = FALSE, quote=  FALSE, sep = "\t")

######################################################################################################################################################################################################

multi_ext = multi_ext %>% select(gene_id, r)
colnames(multi_ext)[2] = "multi"
predixcan_ext = predixcan_ext %>% select(gene_id, r)
colnames(predixcan_ext)[2] = "predixcan"
utmost_ext = utmost_ext %>% select(gene_id, r)
colnames(utmost_ext)[2] = "utmost"
dapg_ext = dapg_ext %>% select(gene_id, r)
colnames(dapg_ext)[2] = "dapg"
mashr_ext = mashr_ext %>% select(gene_id, r)
colnames(mashr_ext)[2] = "mashr"
epixcan_ext = epixcan_ext %>% select(gene_id, r)
colnames(epixcan_ext)[2] = "epixcan"
dpr_add_ext = dpr_add_ext %>% select(gene_id, r)
colnames(dpr_add_ext)[2] = "dpr_add"
mcmc_add_ext = mcmc_add_ext %>% select(gene_id, r)
colnames(mcmc_add_ext)[2] = "mcmc_add"
blup_ext = blup_ext %>% select(gene_id, r)
colnames(blup_ext)[2] = "blup"
bslmm_ext = bslmm_ext %>% select(gene_id, r)
colnames(bslmm_ext)[2] = "bslmm"

merge = merge(multi_ext, predixcan_ext, by.x = "gene_id", by.y = "gene_id")
merge = merge(merge, utmost_ext, by.x = "gene_id", by.y = "gene_id")
merge =	merge(merge, dapg_ext, by.x = "gene_id", by.y = "gene_id")
merge =	merge(merge, mashr_ext, by.x = "gene_id", by.y = "gene_id")
merge =	merge(merge, epixcan_ext, by.x = "gene_id", by.y = "gene_id")
merge =	merge(merge, dpr_add_ext, by.x = "gene_id", by.y = "gene_id")
merge =	merge(merge, mcmc_add_ext, by.x = "gene_id", by.y = "gene_id")
merge =	merge(merge, blup_ext, by.x = "gene_id", by.y = "gene_id")
merge = merge(merge, bslmm_ext, by.x = "gene_id", by.y = "gene_id")

x = merge %>% filter(multi > predixcan) %>% filter(multi > utmost) %>% filter(multi > dapg) %>% filter(multi > mashr) %>% filter(multi > epixcan) %>% filter(multi > dpr_add) %>% filter(multi > mcmc_add) %>% filter(multi > blup) %>% filter(multi > bslmm)

sig_list = intersect(multi_sig, intersect(predixcan_sig, intersect(utmost_sig, intersect(dapg_sig, intersect(mashr_sig, intersect(epixcan_sig, intersect(dpr_add_sig, intersect(mcmc_add_sig, intersect(blup_sig, bslmm_sig)))))))))

y = x[x$gene_id %in% sig_list, ]
y %>% filter(multi - predixcan > 0.1)
