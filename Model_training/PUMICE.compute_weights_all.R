args = commandArgs(trailingOnly=TRUE)

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

## Function for processing input in Constant method
## args
## geno: genotype dataframe
## expression: expression dataframe
## chr: chromosome number
## window: +/- (window) kb to consider around each gene
## values
## origin_id: original map location
## gene_id: ensembl gene id
## snp_list: list of semicolon-separated snp_id
constant_process <- function(geno, expression, chr, window){
	snp_range = makeGRangesFromDataFrame(geno, seqnames.field=c("chromosome"),start.field="start", end.field=c("end"), ignore.strand=TRUE)
	gene_range = makeGRangesFromDataFrame(expression, seqnames.field=c("chromosome"),start.field="start", end.field=c("end"), ignore.strand=TRUE)
	gene_range = resize(gene_range, width(gene_range) + window*1000, fix = 'start')
	gene_range = resize(gene_range, width(gene_range) + window*1000, fix = 'end')
	range = as.data.frame(ranges(gene_range))
	name = as.data.frame(cbind(range$names, paste(chr_num, "_", range$start, "_", range$end, sep="")))
	colnames(name) = c("gene_id", "origin_id")

	fo = findOverlaps(query=snp_range, subject=gene_range, type = "within")
	df = as.data.frame(cbind(rownames(expression[subjectHits(fo),]), rownames(geno[queryHits(fo),])))
	colnames(df) = c("gene_id", "snp_id")
        df$gene_id = sapply(strsplit(df$gene_id, "\\."), function(x){as.character(x[1])})
	df$snp_id = sapply(strsplit(df$snp_id, "\\."), function(x){as.character(x[1])})
	df = df %>% group_by(gene_id) %>% mutate(snp_list = paste0(snp_id, collapse = ";")) %>% select(gene_id, snp_list) %>% distinct(gene_id, .keep_all =TRUE)
	
	merge(x = name, y = df, by.x = "gene_id", by.y = "gene_id") %>% select(origin_id, gene_id, snp_list)

}

## Function to process method 3D (Domain, Loop, TAD, FIRE)
## args
## geno: genotype dataframe
## expression: expression dataframe
## func: 3D map dataframe
## values
## origin_id: original map location
## gene_id: ensembl gene id
## snp_list: list of semicolon-separated snp_id 
func_process <- function(geno, expression, func){
	snp_range = makeGRangesFromDataFrame(geno, seqnames.field=c("chromosome"),start.field="start", end.field=c("end"), ignore.strand=TRUE)
        gene_range = makeGRangesFromDataFrame(expression, seqnames.field=c("chromosome"),start.field="start", end.field=c("end"), ignore.strand=TRUE)
	func_range = makeGRangesFromDataFrame(func, seqnames.field=c("chromosome"),start.field="start", end.field=c("end"), ignore.strand=TRUE)

	fo_func_snp = findOverlaps(query=snp_range, subject=func_range, type = "within")
	df_func_snp = as.data.frame(cbind( rownames(func[subjectHits(fo_func_snp),]), rownames(geno[queryHits(fo_func_snp),]) ))
	colnames(df_func_snp) = c("func_id", "snp_id")
	df_func_snp$func_id = sapply(strsplit(df_func_snp$func_id, "\\."), function(x){as.character(x[1])})
	df_func_snp$snp_id = sapply(strsplit(df_func_snp$snp_id, "\\."), function(x){as.character(x[1])})
	df_func_snp = df_func_snp %>% group_by(func_id) %>% mutate(snp_list = paste0(snp_id, collapse = ";")) %>% distinct(func_id, .keep_all = TRUE) %>% select(func_id, snp_list)	

	fo_func_gene = findOverlaps(query=gene_range, subject=func_range, type = "within")
	df_func_gene = as.data.frame(cbind( rownames(func[subjectHits(fo_func_gene),]), rownames(expression[queryHits(fo_func_gene),]) ))
	colnames(df_func_gene) = c("func_id", "gene_id")
        df_func_gene$func_id = sapply(strsplit(df_func_gene$func_id, "\\."), function(x){as.character(x[1])})
        df_func_gene$gene_id = sapply(strsplit(df_func_gene$gene_id, "\\."), function(x){as.character(x[1])})
	df_func_gene = df_func_gene %>% group_by(func_id) %>% mutate(gene_list = paste0(gene_id, collapse = ";")) %>% distinct(func_id, .keep_all = TRUE) %>% select(func_id, gene_list)

	df_gene_snp = merge(x = df_func_gene, y = df_func_snp, by.x = "func_id", by.y = "func_id")
	s = strsplit(df_gene_snp$gene_list, split = ";")

	unique(data.frame( origin_id = rep(df_gene_snp$func_id, sapply(s, length)), gene_id = unlist(s), snp_list = rep(df_gene_snp$snp_list, sapply(s, length)) ))

}

## Function to process method 3D (pcHiC)
## args
## geno: genotype dataframe
## expression: expression dataframe
## func: 3D map dataframe
## values
## origin_id: original map location
## gene_id: ensembl gene id
## snp_list: list of semicolon-separated snp_id
pchic_process <- function(geno, expression, func){
	snp_range = makeGRangesFromDataFrame(geno, seqnames.field=c("chromosome"),start.field="start", end.field=c("end"), ignore.strand=TRUE)
	func_range = makeGRangesFromDataFrame(func, seqnames.field=c("chromosome"),start.field="start", end.field=c("end"), ignore.strand=TRUE)
	fo_func_snp = findOverlaps(query=snp_range, subject=func_range, type = "within")
        df_gene_snp = as.data.frame(cbind( func[subjectHits(fo_func_snp),], rownames(geno[queryHits(fo_func_snp),]) )) %>% select(-c("chromosome", "start", "end"))
	colnames(df_gene_snp) = c("origin_id", "gene_id", "snp_id")
	df_gene_snp$snp_id = sapply(strsplit(df_gene_snp$snp_id, "\\."), function(x){as.character(x[1])})
	df_gene_snp$gene_id = sapply(strsplit(df_gene_snp$gene_id, "\\."), function(x){as.character(x[1])})
	df_gene_snp = df_gene_snp %>% group_by(gene_id) %>% mutate(snp_list = paste0(snp_id, collapse = ";")) %>% distinct(gene_id, .keep_all = TRUE) %>% select(origin_id, gene_id, snp_list)
	subset(df_gene_snp, gene_id %in% rownames(expression))
}

##tuning is interchangable with training, and testing is interchangable with validating##
geno_scale <- function(geno_subset, sample_tuning, sample_testing){
	##Normalized genotype according to the tuning set##
        geno_tuning = geno_subset %>% select(sample_tuning)
	geno_tuning = geno_fill(geno_tuning)
	center_x_tune = rowMeans(geno_tuning, na.rm = T)
        sd_x_tune = apply(geno_tuning, 1, sd, na.rm = T)
        geno_tuning = as.data.frame(t(scale(t(geno_tuning), center = center_x_tune, scale = sd_x_tune)))
	geno_tuning = geno_tuning[which(!(rowSums(is.na(geno_tuning)) == ncol(geno_tuning))),]
	
	geno_testing = geno_subset %>% select(sample_testing)
	geno_testing = geno_testing[rownames(geno_testing) %in% rownames(geno_tuning),]
	center_x_tune = center_x_tune[rownames(geno_tuning)]
	sd_x_tune = sd_x_tune[rownames(geno_tuning)]
	geno_testing = as.data.frame(t(scale(t(data.matrix(geno_testing)), center = center_x_tune, scale = sd_x_tune)))
	geno_testing[is.na(geno_testing)] <- 0

	list(geno_tuning, geno_testing)
}

##tuning is interchangable with	training, and testing is interchangable	with validating##
expression_scale <- function(expression_subset, sample_tuning, sample_testing){
	##Normalize expression according to the training set##
	expression_tuning = expression_subset %>% select(sample_tuning)
       	center_y_tune = rowMeans(expression_tuning, na.rm = T)
	expression_tuning = as.data.frame(t(scale(t(expression_tuning), center = center_y_tune, scale = FALSE)))

        expression_testing = expression_subset %>% select(sample_testing)
	expression_testing = as.data.frame(t(scale(t(expression_testing), center = center_y_tune, scale = FALSE)))

	list(expression_tuning, expression_testing)
}

elnet <- function(geno_tuning, expression_tuning, geno_testing, expression_testing, penalty_k, lambda_list) {

        ##Fit tuning set and predict in testing sample
        testing <- tryCatch({
                alpha.fit <- glmnet(t(geno_tuning), t(expression_tuning), alpha = 0.5, family="gaussian", penalty.factor = penalty_k, standardize = F, lambda = lambda_list)
                alpha.predicted_testing <- predict(alpha.fit, newx = t(geno_testing))
		alpha.predicted_testing
		},
                error = function(cond) {
                alpha.predicted_testing <- as.matrix(rep(mean(unlist(expression_tuning)), ncol(geno_testing)))
		alpha.predicted_testing
                })

        rownames(testing) = colnames(geno_testing)
        colnames(testing) = 1:length(lambda_list)
	testing
}

#########################################################################################################################
#################################################COMMAD_LINE_INPUT#######################################################
#########################################################################################################################


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

#########################################################################################################################
                                        ##Processing genotype##
#########################################################################################################################
chr_num = args[1]
tissue = args[2]
fold = as.numeric(args[3])
dir_path = normalizePath(args[4])

geno_path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/input/genotype/GTEx_chr", chr_num, "_processed.traw", sep="")
expression_path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/input/expression_processed/", tissue, "_normalized_expression_final.txt",  sep="")

##Process genotype input##
geno = read.table(geno_path, header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE) %>% select(-c(1,3,4,5,6))
colnames(geno)[2:ncol(geno)] = sapply(strsplit(names(geno)[2:ncol(geno)], "_"), function(x){as.character(x[1])})
geno = geno %>% remove_rownames %>% column_to_rownames(var="SNP")

##Subset samples according to the expression file##
geno_sample = read.table("/gpfs/group/dxl46/default/private/poom/GTEx/input/GTEx_EUR_sample_exp.txt", header = FALSE, na.string ='.', as.is=TRUE,check.names=FALSE)
geno = geno %>% select(geno_sample$V1)

##Create subjectID-sample ID conversion dataframe for cross-validation step##
sample_conversion = as.data.frame(cbind( "sample_ID" = seq(1, ncol(geno)), "subject_ID" = colnames(geno)))
names(geno) = sample_conversion$sample_ID[match(names(geno), sample_conversion$subject_ID)]

#########################################################################################################################
                                        ##Processing expression##
#########################################################################################################################
##Process expression input##
expression = read.table(expression_path, header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE) %>% filter(chromosome == chr_num)
names(expression)[5:ncol(expression)] = sample_conversion$sample_ID[match(names(expression)[5:ncol(expression)], sample_conversion$subject_ID)]

#########################################################################################################################
                                        ##Processing overall##
#########################################################################################################################
##Find sample overlap between genotype and expression data##
sample_overlap = intersect(colnames(geno), colnames(expression))

##Update the sample list in both genotype and expression data##
geno = geno %>% select(sample_overlap)
expression = expression %>% select(gene_id, chromosome, start, end, sample_overlap) %>% remove_rownames %>% column_to_rownames(var="gene_id")

##Filling in SNP coordinates##
geno = geno_coord(geno)

#########################################################################################################################
                                        ##Import result from internal validation##
#########################################################################################################################
##Import result from cv_int_total from SCREEN##
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/cv_int_total/screen_r/output/", tissue, sep = "")
filename = dir(path, pattern = paste("_chr", chr_num, "_", sep = ""))

screen_int_total = data.frame()
for (i in 1:length(filename)){
        print(i)
        try({
	screen_int_temp_total = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE)
        screen_int_total = rbind(screen_int_total, screen_int_temp_total)
        })
}
screen_tissue_int_total = screen_int_total %>% filter(penalty != "NA") %>% group_by(gene_id) %>% dplyr::slice(which.min(cvm_v))
predixcan_tissue_int_total = screen_int_total %>% filter(type == "1000") %>% filter(penalty == "7")

##Write out internal validation result##
screen_iv_output = cbind("tissue" = rep(tissue, nrow(screen_tissue_int_total)), screen_tissue_int_total %>% select(gene_id, type, region, penalty, pear_avg_t, pear_stouffer_pval, spear_avg_t, spear_stouffer_pval))
predixcan_iv_output = cbind("tissue" = rep(tissue, nrow(predixcan_tissue_int_total)), predixcan_tissue_int_total %>% select(gene_id, type, region, penalty, pear_avg_t, pear_stouffer_pval, spear_avg_t, spear_stouffer_pval))

iv_path = paste(dir_path, "/multi-omics/internal_validation/", tissue, "/", tissue, 
		"_chr", chr_num, "_GTEx_multi-omics_internal_validation.txt", sep = "")
write.table(screen_iv_output, iv_path, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

iv_path = paste(dir_path, "/predixcan/internal_validation/", tissue, "/", tissue,
                "_chr", chr_num, "_GTEx_predixcan_internal_validation.txt", sep = "")
write.table(predixcan_iv_output, iv_path, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


#########################################################################################################################
                                        ##Import ENCODE annotation##
#########################################################################################################################
##Import ENCODE file##
penalties <- tryCatch({
        encode_path <- paste("/gpfs/group/dxl46/default/private/poom/GTEx/input/screen/", tissue ,"/GTEx_screen_annotation_chr", chr_num, ".txt", sep="")
        encode <- read.table(encode_path, header = FALSE, na.string ='.', as.is=TRUE,check.names=FALSE)
        colnames(encode) <- c("chromosome", "start", "end", "overlap")
        encode$chromosome <- chr_num
        encode$chromosome <- as.character(encode$chromosome)
        encode$start <- as.character(encode$start)
        encode$end <- as.character(encode$end)
        encode$code <- paste(encode$chromosome, encode$start, encode$end, sep = ":")
        encode <- encode %>% select(code, overlap)

        ##Convert encode bed position to snp name##
        geno_code <- geno %>% select(chromosome, start, end)
        geno_code$chromosome <- as.character(geno_code$chromosome)
        geno_code$start <- as.character(geno_code$start)
        geno_code$end <- as.character(geno_code$end)
        geno_code$code <- paste(geno_code$chromosome, geno_code$start, geno_code$end, sep = ":")
        geno_code$snp <- rownames(geno)

        encode <- distinct(merge(x = encode, y = geno_code, by.x = "code", by.y = "code") %>% select(snp, overlap))
        encode$overlap <- as.numeric(as.character(encode$overlap))
        encode <- encode[match( rownames(geno), encode$snp ),]

        ##Penalty1 is based on cut off of number of annotation and no penalization on established predictors##
        penalty1 <- c(rep(1, nrow(geno)))
        penalty1[which(encode$overlap > 0)] <- 0

        ##Penalty2 is based on cut off of number of annotation and set fixed penalty at 1/6##
        penalty2 <- c(rep(1, nrow(geno)))
        penalty2[which(encode$overlap > 0)] <- 1/6

        ##Penalty3 is based on cut off of number of annotation and set fixed penalty at 2/6##
        penalty3 <- c(rep(1, nrow(geno)))
        penalty3[which(encode$overlap > 0)] <- 2/6

        ##Penalty4 is based on cut off of number of annotation and set fixed penalty at 3/6##
        penalty4 <- c(rep(1, nrow(geno)))
        penalty4[which(encode$overlap > 0)] <- 3/6

        ##Penalty5 is based on cut off of number of annotation and set fixed penalty at 4/6##
        penalty5 <- c(rep(1, nrow(geno)))
        penalty5[which(encode$overlap > 0)] <- 4/6

	##Penalty6 is based on cut off of number of annotation and set fixed penalty at 5/6##
        penalty6 <- c(rep(1, nrow(geno)))
        penalty6[which(encode$overlap > 0)] <- 5/6

        ##Penalty7 is elastic net##
        penalty7 <- c(rep(1, nrow(geno)))
        penalty7[which(encode$overlap > 0)] <- 1

        rbind(penalty1, penalty2, penalty3, penalty4, penalty5, penalty6, penalty7)
        },
	error = function(cond) {
        t(as.matrix(c(rep(1,nrow(geno)))))
        })
colnames(penalties) = rownames(geno)
n_penalties = nrow(penalties)

#########################################################################################################################
                                        ##Create CV folds##
#########################################################################################################################
##Create cv fold
set.seed(123)
cv_fold_total = createFolds(sample_conversion$sample_ID, k = fold, list = TRUE, returnTrain = FALSE)
cv_fold = data.frame()
for (f in 1:fold) {
        temp = cbind(sample = intersect( sample_overlap, sample_conversion$sample_ID[cv_fold_total[[f]]] ), foldid = f)
        cv_fold = rbind(cv_fold, temp)
}

#########################################################################################################################
                                        ##Process Multi-omics##
#########################################################################################################################
coef = data.frame()
summary = data.frame()
for (i in 1:nrow(screen_tissue_int_total)) {
	tryCatch({
	func_type = screen_tissue_int_total$type[i]
	expression_subset = expression[rownames(expression) %in% screen_tissue_int_total$gene_id[i],]	

	if (func_type == "pchic" || func_type == "pcHiC") {
		func_path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/input/3D_genome_processed/", func_type, "_input/", tissue, "_", func_type, ".txt", sep = "")
	        func = read.table(func_path, header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE) %>% filter(chromosome == chr_num)
		func = func %>% select(-c("length", "MinusLog10Pval"))
		map = pchic_process(geno, expression_subset, func)
	} else if (func_type == "Loop" || func_type == "TAD" || func_type == "Domain") {
		func_path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/input/3D_genome_processed/", func_type, "_input/", tissue, "_", func_type, ".txt", sep = "")		
		func = read.table(func_path, header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE) %>% filter(chromosome == chr_num)
		func$func_id = paste(func$chromosome, func$start, func$end, sep="_")
                func = func %>% filter(func_id == screen_tissue_int_total$region[i]) %>% remove_rownames %>% column_to_rownames(var="func_id")
                map = func_process(geno, expression_subset, func)
	} else {
		func_type = as.numeric(func_type)
		map = constant_process(geno, expression_subset, chr_num, func_type)
	}
	
	geno_subset = geno[unlist(strsplit(as.character(map[1,3]), ";")),] %>% select(-c("chromosome", "start", "end"))
	expression_subset = expression_subset %>% select(-c("chromosome", "start", "end"))
	origin_id = as.character(map[1,1])

	##Need to create a list of predefined lambdas##	
	##Standardized genotype matrix##
	geno_subset_all = geno_fill(geno_subset)
	center_x_train_all = rowMeans(geno_subset_all, na.rm = T)
        sd_x_train_all = apply(geno_subset_all, 1, sd, na.rm = T)
        geno_subset_all = as.data.frame(t(scale(t(geno_subset_all), center = center_x_train_all, scale = sd_x_train_all)))
        geno_subset_all = geno_subset_all[which(!(rowSums(is.na(geno_subset_all)) == ncol(geno_subset_all))),]

	##Center expression matrix##
        center_y_train = rowMeans(expression_subset, na.rm = T)
        expression_subset_all = as.data.frame(t(scale(t(expression_subset), center = center_y_train, scale = FALSE)))

	##Retrieve penalty coefficients##
	penalty_k = penalties[screen_tissue_int_total$penalty[i], rownames(geno_subset_all)]

	##Find a list of lambda##
	lambda_list =  (glmnet(t(geno_subset_all), t(expression_subset_all), alpha = 0.5, family="gaussian", penalty.factor = penalty_k, standardize = F))$lambda	
	
	prediction_testing = data.frame()
	for (f in 1:fold) {
		print(f)
                sample_testing = as.list(cv_fold %>% filter(foldid == f) %>% select(sample))$sample
		sample_tuning = as.list(cv_fold %>% filter(foldid != f) %>% select(sample))$sample

		##Scale the predictors accordingly##
		geno_output = geno_scale(geno_subset, sample_tuning, sample_testing)
                geno_tuning = geno_output[[1]]
                geno_testing = geno_output[[2]]

		##Scale the expression accordingly##
		expression_output = expression_scale(expression_subset, sample_tuning, sample_testing)
                expression_tuning = expression_output[[1]]
                expression_testing = expression_output[[2]]

		if (nrow(geno_testing) > 1) {
			alpha.predicted_testing = elnet(geno_tuning, expression_tuning, geno_testing, expression_testing, penalty_k, lambda_list)
                        prediction_testing = rbind(prediction_testing, cbind( t(expression_testing), alpha.predicted_testing ) )
		} 
	}
	prediction_testing = prediction_testing[match( colnames(expression_subset), rownames(prediction_testing) ),]
	
	##Identify best lambda value##
	cvm_testing = sapply(2:ncol(prediction_testing), function(x) { mean( (unlist( prediction_testing[,x] ) - unlist( prediction_testing[,1] ))^2 ) } )
	best_lambda = lambda_list[which.min(cvm_testing)]

	##Create final model##
	alpha.fit = glmnet(t(geno_subset_all), t(expression_subset_all), alpha = 0.5, family="gaussian", penalty.factor = penalty_k, standardize = F, lambda = best_lambda)
	
	if ( all(alpha.fit$beta == 0) ) {
		summary_temp = data.frame( "tissue" = tissue, "gene_id" = rownames(expression_subset_all), 
					   "func_type" = func_type, "region" = origin_id, fit = "elnet", 
					   "alpha" = 0.5, "penalty" = screen_tissue_int_total$penalty[i],
					   "input_snps" = nrow(geno_subset_all), "nonzero_snps" = 0)
		summary = rbind(summary, summary_temp)
	} else {
		coef_temp = as.matrix(alpha.fit$beta[which(alpha.fit$beta[,1] != 0),])
		rownames(coef_temp) = rownames(coef(alpha.fit, s = 'lambda.min'))[coef(alpha.fit, s = 'lambda.min')[,1]!= 0][-1]
		coef_output = data.frame( "tissue" = rep(tissue, nrow(coef_temp)), "gene_id" = rep(rownames(expression_subset_all), nrow(coef_temp)), 
					  "snp" = rownames(coef_temp), "weight" = coef_temp[,1]) %>% remove_rownames
		coef = rbind(coef, coef_output)

		summary_temp = data.frame( "tissue" = tissue, "gene_id" = rownames(expression_subset_all),  
                                           "func_type" = func_type, "region" = origin_id, fit = "elnet",  
                                           "alpha" = 0.5, "penalty" = screen_tissue_int_total$penalty[i], 
                                           "input_snps" = nrow(geno_subset_all), "nonzero_snps" = length(alpha.fit$beta[which(alpha.fit$beta[,1] != 0),]))
		summary = rbind(summary, summary_temp)		
	}
	}, error=function(e){})

	filename_summ_result <- paste(dir_path, "/multi-omics/summary/", tissue, "/", tissue, "_chr", chr_num, "_GTEx_multi-omics_summary.txt", sep="")
	write.table(summary, filename_summ_result, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

	filename_coef_result <- paste(dir_path, "/multi-omics/output/", tissue, "/", tissue, "_chr", chr_num, "_GTEx_multi-omics_coef.txt", sep="")
	write.table(coef, filename_coef_result, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
dummy="I am done"
filename_dummy = paste(dir_path, "/multi-omics/dummy/", tissue, "/dummy_testing_chr", chr_num, ".txt", sep="" )
write.table(dummy, filename_dummy, quote=FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#########################################################################################################################
                                        ##Process PrediXcan##
#########################################################################################################################
coef = data.frame()
summary = data.frame()
for (i in 1:nrow(predixcan_tissue_int_total)) {
	tryCatch({
        func_type = as.numeric(predixcan_tissue_int_total$type[i])
        expression_subset = expression[rownames(expression) %in% predixcan_tissue_int_total$gene_id[i],]
	map = constant_process(geno, expression_subset, chr_num, func_type)

	geno_subset = geno[unlist(strsplit(as.character(map[1,3]), ";")),] %>% select(-c("chromosome", "start", "end"))
        expression_subset = expression_subset %>% select(-c("chromosome", "start", "end"))
        origin_id = as.character(map[1,1])
	
	##Need to create a list of predefined lambdas##
        ##Standardized genotype matrix##
        geno_subset_all = geno_fill(geno_subset)
        center_x_train_all = rowMeans(geno_subset_all, na.rm = T)
        sd_x_train_all = apply(geno_subset_all, 1, sd, na.rm = T)
        geno_subset_all = as.data.frame(t(scale(t(geno_subset_all), center = center_x_train_all, scale = sd_x_train_all)))
        geno_subset_all = geno_subset_all[which(!(rowSums(is.na(geno_subset_all)) == ncol(geno_subset_all))),]

        ##Center expression matrix##
        center_y_train = rowMeans(expression_subset, na.rm = T)
        expression_subset_all = as.data.frame(t(scale(t(expression_subset), center = center_y_train, scale = FALSE)))

        ##Retrieve penalty coefficients##
        penalty_k = penalties[predixcan_tissue_int_total$penalty[i], rownames(geno_subset_all)]

        ##Find a list of lambda##
        lambda_list =  (glmnet(t(geno_subset_all), t(expression_subset_all), alpha = 0.5, family="gaussian", penalty.factor = penalty_k, standardize = F))$lambda

        prediction_testing = data.frame()
	for (f in 1:fold) {
                print(f)
                sample_testing = as.list(cv_fold %>% filter(foldid == f) %>% select(sample))$sample
                sample_tuning = as.list(cv_fold %>% filter(foldid != f) %>% select(sample))$sample

                ##Scale the predictors accordingly##
                geno_output = geno_scale(geno_subset, sample_tuning, sample_testing)
                geno_tuning = geno_output[[1]]
                geno_testing = geno_output[[2]]

                ##Scale the expression accordingly##
                expression_output = expression_scale(expression_subset, sample_tuning, sample_testing)
                expression_tuning = expression_output[[1]]
                expression_testing = expression_output[[2]]

                if (nrow(geno_testing) > 1) {
                        alpha.predicted_testing = elnet(geno_tuning, expression_tuning, geno_testing, expression_testing, penalty_k, lambda_list)
                        prediction_testing = rbind(prediction_testing, cbind( t(expression_testing), alpha.predicted_testing ) )
                }
        }
	prediction_testing = prediction_testing[match( colnames(expression_subset), rownames(prediction_testing) ),]

        ##Identify best lambda value##
        cvm_testing = sapply(2:ncol(prediction_testing), function(x) { mean( (unlist( prediction_testing[,x] ) - unlist( prediction_testing[,1] ))^2 ) } )
        best_lambda = lambda_list[which.min(cvm_testing)]

        ##Create final model##
        alpha.fit = glmnet(t(geno_subset_all), t(expression_subset_all), alpha = 0.5, family="gaussian", penalty.factor = penalty_k, standardize = F, lambda = best_lambda)

        if ( all(alpha.fit$beta == 0) ) {
                summary_temp = data.frame( "tissue" = tissue, "gene_id" = rownames(expression_subset_all),
                                           "func_type" = func_type, "region" = origin_id, fit = "elnet",
                                           "alpha" = 0.5, "penalty" = predixcan_tissue_int_total$penalty[i],
                                           "input_snps" = nrow(geno_subset_all), "nonzero_snps" = 0)
		summary = rbind(summary, summary_temp)
        } else {
                coef_temp = as.matrix(alpha.fit$beta[which(alpha.fit$beta[,1] != 0),])
		rownames(coef_temp) = rownames(coef(alpha.fit, s = 'lambda.min'))[coef(alpha.fit, s = 'lambda.min')[,1]!= 0][-1]
                coef_output = data.frame( "tissue" = rep(tissue, nrow(coef_temp)), "gene_id" = rep(rownames(expression_subset_all), nrow(coef_temp)),
                                          "snp" = rownames(coef_temp), "weight" = coef_temp[,1]) %>% remove_rownames
                coef = rbind(coef, coef_output)

                summary_temp = data.frame( "tissue" = tissue, "gene_id" = rownames(expression_subset_all),
                                           "func_type" = func_type, "region" = origin_id, fit = "elnet",
                                           "alpha" = 0.5, "penalty" = predixcan_tissue_int_total$penalty[i],
                                           "input_snps" = nrow(geno_subset_all), "nonzero_snps" = length(alpha.fit$beta[which(alpha.fit$beta[,1] != 0),]))
                summary = rbind(summary, summary_temp)
        }
	}, error=function(e){})

	filename_summ_result <- paste(dir_path, "/predixcan/summary/", tissue, "/", tissue, "_chr", chr_num, "_GTEx_predixcan_summary.txt", sep="")
        write.table(summary, filename_summ_result, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

        filename_coef_result <- paste(dir_path, "/predixcan/output/", tissue, "/", tissue, "_chr", chr_num, "_GTEx_predixcan_coef.txt", sep="")
        write.table(coef, filename_coef_result, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
dummy="I am done"
filename_dummy = paste(dir_path, "/predixcan/dummy/", tissue, "/dummy_testing_chr", chr_num, ".txt", sep="" )
write.table(dummy, filename_dummy, quote=FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
