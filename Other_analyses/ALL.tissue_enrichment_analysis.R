args = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(stringr)

meth = args[1]
tissue_list = c("Muscle_Skeletal", "Skin_Sun_Exposed_Lower_leg", "Thyroid", "Adipose_Subcutaneous", "Lung", "Artery_Tibial",
                        "Whole_Blood", "Nerve_Tibial", "Esophagus_Mucosa", "Skin_Not_Sun_Exposed_Suprapubic", "Esophagus_Muscularis", "Adipose_Visceral_Omentum",
                        "Cells_Transformed_fibroblasts", "Artery_Aorta", "Heart_Left_Ventricle", "Heart_Atrial_Appendage", "Breast_Mammary_Tissue", "Colon_Transverse",
                        "Stomach", "Testis", "Esophagus_Gastroesophageal_Junction", "Colon_Sigmoid", "Pancreas", "Pituitary",
                        "Adrenal_Gland", "Brain_Cerebellum", "Liver", "Brain_Caudate_basal_ganglia", "Artery_Coronary", "Brain_Cortex",
                        "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Spleen", "Prostate", "Brain_Frontal_Cortex_BA9", "Brain_Anterior_cingulate_cortex_BA24",
                        "Small_Intestine_Terminal_Ileum", "Brain_Hippocampus", "Brain_Putamen_basal_ganglia", "Brain_Hypothalamus", "Ovary", "Cells_EBV-transformed_lymphocytes",
                        "Vagina", "Uterus", "Brain_Amygdala", "Brain_Spinal_cord_cervical_c-1", "Minor_Salivary_Gland", "Brain_Substantia_nigra")

pheno_list = setdiff(dir("/gpfs/group/dxl46/default/private/poom/espresso/sumstats_clean_new"), "check.R")

result = read.table( paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", meth, "_sig_loci.txt", sep = ""), header = TRUE)
#result = result %>% filter(pheno != "gwg") %>% filter(pheno != "gest_dur") 

obs = as.data.frame(matrix(0, length(tissue_list), length(pheno_list)))

for (i in 1:length(pheno_list)) {
        for (j in 1:length(tissue_list)) {
                obs_pheno = result %>% filter(pheno == pheno_list[i])
                idx_row = which(obs_pheno$tissue == tissue_list[j])
                obs[j, i] = obs_pheno$num[idx_row]
        }
}
rownames(obs) = tissue_list
colnames(obs) = pheno_list
colMax <- function(data) sapply(data, max, na.rm = TRUE)
x = colMax(obs)
obs = obs[, which(x > 1)]

##Calcualte Pearson's residuals##
p_residual = chisq.test(obs)$residuals
p_chisquared = p_residual^2
p_chisquared = p_chisquared[, colSums(is.na(p_chisquared)) == 0]
pchisq(sum(p_chisquared), df=(nrow(p_chisquared)-1)*(ncol(p_chisquared)-1), lower.tail=FALSE)

#############################################################################################################################################################
##Calculate adjusted Pearson's residual##

##Calculated expected values for Pearson##
obs_rowsum = rowSums(obs)
obs_colsum = colSums(obs)
exp = as.data.frame(matrix(0, length(tissue_list), ncol(obs)))

for (i in 1:ncol(obs)) {
        for (j in 1:nrow(obs)) {
                exp[j, i] = obs_colsum[i]*obs_rowsum[j]/sum(obs_rowsum)
        }
}
rownames(exp) = rownames(obs)
colnames(exp) = colnames(obs)

##Calculate observed ratios##
obs_p_row = obs_rowsum/sum(obs_rowsum)
obs_p_col = obs_colsum/sum(obs_colsum)

##Calculate residuals##
adjp_residual = as.data.frame(matrix(0, nrow(obs), ncol(obs)))

for (i in 1:ncol(obs)) {
        for (j in 1:nrow(obs)) {
                adjp_residual[j,i] = (obs[j,i] - exp[j,i])/(sqrt(exp[j,i]*(1-obs_p_row[j])*(1-obs_p_col[i])))
        }
}
rownames(adjp_residual) = rownames(obs)
colnames(adjp_residual) = colnames(obs)

##Calculate chi-square statistic##
adjp_chisquared = adjp_residual^2
adjp_chisquared = adjp_chisquared[, colSums(is.na(adjp_chisquared)) == 0]
pchisq(sum(adjp_chisquared), df=(nrow(adjp_chisquared)-1)*(ncol(adjp_chisquared)-1), lower.tail=FALSE)

#############################################################################################################################################################

save(adjp_residual, file = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", meth, "_residual.RDat", sep ="") )
