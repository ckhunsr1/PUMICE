args = commandArgs(trailingOnly=TRUE)

######################################################################################################################################
##Created list of observed vs. predicted expression##

## Function for subsetting gene sets
## args
## N: number of total genes
## total_file_num: number of total files
## values
## idx: matrix of fold by 2 with first col being starting index and second col being ending index
file_helper <- function(N, total_file_num){
        valid_num = ceiling(N/total_file_num)
        idx1 = seq(1,N,valid_num)
        idx2 = c(idx1[-1]-1,N)
        cbind(idx1,idx2)
}


## Function for filling in missing dosage
## args
## geno: genotype dataframe
## values
## geno_new: genotype dataframe with filled in dosage
geno_fill <- function(geno, type){

        ##Filling in missing data##
        geno = data.matrix(geno)
        class(geno) = "numeric"
        k = which(is.na(geno), arr.ind=TRUE)
        geno[k] <- rowMeans(geno, na.rm=TRUE)[k[,1]]
        as.data.frame(geno)
}

## Function for filling in snp coordinates
## args
## geno: genotype dataframe
## values
## geno_new: genotype dataframe with snp coordinates
geno_coord <- function(geno){
        geno$chromosome <- sapply(strsplit(rownames(geno)[1:nrow(geno)], "_"), function(x){as.character(x[1])})
        geno$start <- sapply(strsplit(rownames(geno)[1:nrow(geno)], "_"), function(x){as.numeric(x[2])}) - 1
        start_length = as.data.frame(sapply(sapply(strsplit(rownames(geno)[1:nrow(geno)], "_"), function(x){as.character(x[3])}), nchar))
        end_length = as.data.frame(sapply(sapply(strsplit(rownames(geno)[1:nrow(geno)], "_"), function(x){as.character(x[4])}), nchar))
        df = cbind(start_length, end_length)
        colnames(df)[1] = "start_length"
        colnames(df)[2] = "end_length"
        df$max = apply(df, 1, max)
        rm(start_length, end_length)
        geno$end <- geno$start + df$max
        rm(df)
	geno
}

#########################################################################################################################
                                        ##Import library##
#########################################################################################################################
options(stringsAsFactors=F)
library(tidyr)
library(dplyr)
library(IRanges)
library(GenomicRanges)
library(glmnet)
library(tidyverse)
library(genefilter)
library(caret)

chr_num = args[1]
tissue = args[2]
type = args[3]
ext = args[4] ##one of the following DGN, CMC, GEUVADIS##
gene = args[5]

geno1_path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/input/genotype/GTEx_chr", chr_num, "_processed.traw", sep="")
geno2_path = paste("/gpfs/group/dxl46/default/private/poom/", ext, "/input/genotype/", ext, "_chr", chr_num, "_processed.traw", sep="")

expression1_path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/input/expression_processed/", tissue, "_normalized_expression_final.txt",  sep="")
expression2_path = paste("/gpfs/group/dxl46/default/private/poom/", ext, "/input/expression_processed/", ext, "_normalized_expression_final.txt", sep = "")

#########################################################################################################################
##Process GTEx genotype input##
geno1 = read.table(geno1_path, header = TRUE, na.string = NA, as.is=TRUE,check.names=FALSE) %>% select(-c(1,3,4,5,6))
colnames(geno1)[2:ncol(geno1)] = sapply(strsplit(names(geno1)[2:ncol(geno1)], "_"), function(x){as.character(x[1])})
geno1 = geno1 %>% remove_rownames %>% column_to_rownames(var="SNP")

##Subset samples according to the expression file##
geno1_sample = read.table("/gpfs/group/dxl46/default/private/poom/GTEx/input/GTEx_EUR_sample_exp.txt", header = FALSE, na.string ='.', as.is=TRUE,check.names=FALSE)
geno1 = geno1 %>% select(geno1_sample$V1)

##Create subjectID-sample ID conversion dataframe for cross-validation step##
sample_conversion1 = as.data.frame(cbind( "sample_ID" = seq(1, ncol(geno1)), "subject_ID" = colnames(geno1)))
names(geno1) = sample_conversion1$sample_ID[match(names(geno1), sample_conversion1$subject_ID)]

##Create genotype file for predixcan##
geno1 = geno_fill(geno1)

#########################################################################################################################
##Process external genotype input##
geno2 = read.table(geno2_path, header = TRUE, na.string =NA, as.is=TRUE,check.names=FALSE) %>% select(-c(1,3,4,5,6))

if (ext == "DGN"){
        colnames(geno2)[2:ncol(geno2)] = sapply(strsplit(names(geno2)[2:ncol(geno2)], "_"), function(x){as.character(x[4])})
} else if (ext == "CMC") {
        colnames(geno2)[2:ncol(geno2)] = sapply(strsplit(names(geno2)[2:ncol(geno2)], "_"), function(x){ paste(as.character(x[3]), as.character(x[4]), as.character(x[5]), sep = "_") })
} else if (ext == "GEUVADIS") {
        colnames(geno2)[2:ncol(geno2)] = sapply(strsplit(names(geno2)[2:ncol(geno2)], "_"), function(x){as.character(x[1])})
}

geno2 = geno2 %>% remove_rownames %>% column_to_rownames(var="SNP")

##Create subjectID-sample ID conversion dataframe for cross-validation step##
sample_conversion2 = as.data.frame(cbind( "sample_ID" = seq(1, ncol(geno2)), "subject_ID" = colnames(geno2)))
names(geno2) = sample_conversion2$sample_ID[match(names(geno2), sample_conversion2$subject_ID)]

##Create genotype file for predixcan##
geno2 = geno_fill(geno2)

#########################################################################################################################
                                        ##Matching SNPs from DGN to GTEx##
#########################################################################################################################
##Find overlap SNP##
snp_overlap = intersect(rownames(geno1), rownames(geno2))
geno2_overlap = geno2[snp_overlap,]

##Find SNP needed to be added and fill in SNP with the average dosage of GTEx total samples##
snp_add = setdiff(rownames(geno1), rownames(geno2))
geno1_avg = as.data.frame(apply(geno1[rownames(geno1) %in% snp_add,], 1, mean))
geno2_add = cbind(geno1_avg, replicate(ncol(geno2)-1, geno1_avg[,1]))
colnames(geno2_add) = 1:dim(geno2)[2]

##Create total geno dataframe for DGN cohort##
geno2 = rbind(geno2_overlap, geno2_add)
geno2 = geno2[match( rownames(geno1), rownames(geno2) ), ]

##Center/None/Standardized genotype matrix##
if (type == "utmost") {
        center_x = rowMeans(geno1, na.rm = T)
        geno1 = as.data.frame(t(scale(t(geno1), center = center_x, scale = FALSE)))
        geno2 = as.data.frame(t(scale(t(geno2), center = center_x, scale = FALSE)))
}else if (type == "espresso") {
        geno1 = geno1
        geno2 = geno2
} else {
        center_x = rowMeans(geno1, na.rm = T)
        sd_x = apply(geno1, 1, sd, na.rm = T)
        geno1 = as.data.frame(t(scale(t(geno1), center = center_x, scale = sd_x)))
        geno2 = as.data.frame(t(scale(t(geno2), center = center_x, scale = sd_x)))
}

#########################################################################################################################
                                        ##Processing expression##
#########################################################################################################################
##Process external expression input##
expression2_path = paste("/gpfs/group/dxl46/default/private/poom/", ext, "/input/expression_processed/", ext, "_normalized_expression_final.txt", sep = "")
expression2 = read.table(expression2_path, header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE) %>% filter(chromosome == chr_num)
names(expression2)[5:ncol(expression2)] = sample_conversion2$sample_ID[match(names(expression2)[5:ncol(expression2)], sample_conversion2$subject_ID)]
sample_overlap2 = intersect(colnames(geno2), colnames(expression2))
expression2 = expression2 %>% select(gene_id, sample_overlap2) %>% remove_rownames %>% column_to_rownames(var="gene_id")


#########################################################################################################################
                                        ##Importing weights and do prediction##
#########################################################################################################################

if ( type == "espresso" ) {
	path = paste("/gpfs/group/dxl46/default/private/poom/espresso/espresso_output/result/", tissue, sep = "")
        setwd(path)
        filename = dir(path, pattern = paste("chr", chr_num, "_", sep = ""))

        weight = data.frame()
        for (i in 1:length(filename)) {
                temp = read.table(filename[i], header = TRUE, na.string = NA, as.is=TRUE, check.names=FALSE)
                weight = rbind(weight, temp)
        }

	##Process weight dataframe##
        weight$snp = paste(str_replace(str_replace(weight$snp, ":", "_"), "/", "_"), "b37" ,sep = "_")
        weight$tissue = tissue
        weight$weight = -1*weight$weight ##since espresso effect allele is alternative allele##
        weight = weight %>% select(tissue, gene, snp, weight)
        colnames(weight)[2] = "gene_id"
} else if ( type == "mashr" ) {
	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/output/", tissue, sep = "")
        setwd(path)
        filename = dir(path)
	weight = read.table(filename, header = TRUE, na.string = NA, as.is=TRUE, check.names=FALSE)
} else {
	path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", type, "/output/", tissue, sep = "")
        setwd(path)
        filename = dir(path, pattern = paste("chr", chr_num, "_", sep = ""))

        weight = data.frame()
        for (i in 1:length(filename)) {
                temp = read.table(filename[i], header = TRUE, na.string = NA, as.is=TRUE, check.names=FALSE)
                weight = rbind(weight, temp)
        }
}
#########################################################################################################################
                                        ##Prediction##
#########################################################################################################################

##Subset expression dataframe##
expression2_subset = expression2[gene, ]

##Filter for snp weight##
weight_subset = weight %>% filter(gene_id == gene) %>% select(snp, weight) %>% remove_rownames %>% column_to_rownames(var="snp")
weight_num = nrow(weight_subset)

##Filter genotype matrix##
geno2_subset = as.matrix(geno2[rownames(geno2) %in% rownames(weight_subset), ])
geno2_num = nrow(geno2_subset)

##Make sure to get the same matrix dimension##
snp_overlap = intersect(rownames(weight_subset), rownames(geno2_subset))
weight_subset = t(weight_subset[rownames(weight_subset) %in% snp_overlap,])
colnames(weight_subset) = snp_overlap
geno2_subset = geno2_subset[rownames(geno2_subset) %in% snp_overlap, ]

##Prediction in testing set##
prediction = as.data.frame(weight_subset %*% geno2_subset)

output = cbind( t(expression2_subset), t(prediction) )
colnames(output) = c("observed", "predicted")

write.table(output, paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", gene, "_", ext, "_", type, ".txt", sep = ""), 
		col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

