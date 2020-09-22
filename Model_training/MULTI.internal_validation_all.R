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

## Function to process method 3D (Subcompartment)
## args
## geno: genotype dataframe
## expression: expression dataframe
## func: 3D map dataframe
## sc: subcompartment selected
## values
## origin_id: original map location
## gene_id: ensembl gene id
## snp_list: list of semicolon-separated snp_id
sc_process <- function(geno, expression, func, sc){
	snp_range = makeGRangesFromDataFrame(geno, seqnames.field=c("chromosome"),start.field="start", end.field=c("end"), ignore.strand=TRUE)
	gene_range = makeGRangesFromDataFrame(expression, seqnames.field=c("chromosome"),start.field="start", end.field=c("end"), ignore.strand=TRUE)
        func_range = makeGRangesFromDataFrame(func, seqnames.field=c("chromosome"),start.field="start", end.field=c("end"), ignore.strand=TRUE)

	fo_func_snp = findOverlaps(query=snp_range, subject=func_range, type = "within")
        df_func_snp = as.data.frame(cbind( func[subjectHits(fo_func_snp),], rownames(geno[queryHits(fo_func_snp),]) ))
        colnames(df_func_snp)[5] = "snp_id"
	df_func_snp$snp_id = sapply(strsplit(df_func_snp$snp_id, "\\."), function(x){as.character(x[1])})
	df_func_snp = df_func_snp %>% select(sc, snp_id) %>% group_by(sc) %>% mutate(snp_list = paste0(snp_id, collapse = ";")) %>% distinct(sc, .keep_all = TRUE) %>% select(sc, snp_list) 

        fo_func_gene = findOverlaps(query=gene_range, subject=func_range, type = "within")
        df_func_gene = as.data.frame(cbind( func[subjectHits(fo_func_gene),], rownames(expression[queryHits(fo_func_gene),]) ))
        colnames(df_func_gene)[5] = "gene_id"
	df_func_gene$gene_id = sapply(strsplit(df_func_gene$gene_id, "\\."), function(x){as.character(x[1])})
	df_func_gene = df_func_gene %>% select(sc, gene_id) %>% group_by(sc) %>% mutate(gene_list = paste0(gene_id, collapse = ";")) %>% distinct(sc, .keep_all = TRUE) %>% select(sc, gene_list)

        df_gene_snp = merge(x = df_func_gene, y = df_func_snp, by.x = "sc", by.y = "sc")
        s = strsplit(df_gene_snp$gene_list, split = ";")
        unique(data.frame( origin_id = rep(df_gene_snp$sc, sapply(s, length)), gene_id = unlist(s), snp_list = rep(df_gene_snp$snp_list, sapply(s, length)) ))
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

elnet <- function(geno_training, expression_training, geno_validating, expression_validating, geno_tuning, expression_tuning, geno_testing, expression_testing, penalty_k) {

	##Fit best lambda in training set and predict in validating sample
        alpha.predicted_validating <- tryCatch({
        	set.seed(123)
		alpha.fit <- cv.glmnet(t(geno_training), t(expression_training), alpha = 0.5, family="gaussian", penalty.factor = penalty_k, standardize = F, nfolds = fold, type.measure='mse')
                predict(alpha.fit, newx = t(geno_validating), s = "lambda.min")
                },
                error = function(cond) {
		as.matrix(rep(mean(unlist(expression_training)), ncol(geno_validating)))
                })
        rownames(alpha.predicted_validating) = colnames(geno_validating)

        ##Fit tuning set and predict in testing sample
        testing <- tryCatch({
		set.seed(123)
                alpha.fit <- cv.glmnet(t(geno_tuning), t(expression_tuning), alpha = 0.5, family="gaussian", penalty.factor = penalty_k, standardize = F, nfolds = fold, type.measure='mse')
                alpha.predicted_testing <- predict(alpha.fit, newx = t(geno_testing), s = "lambda.min")
		alpha_coef = as.matrix(alpha.fit$glmnet.fit$beta[, which.min(alpha.fit$cvm)])                
		list(alpha.predicted_testing, alpha_coef)
		},
                error = function(cond) {
                alpha.predicted_testing <- as.matrix(rep(mean(unlist(expression_tuning)), ncol(geno_testing)))
		alpha_coef = as.matrix(rep(0, nrow(geno_tuning)))
		list(alpha.predicted_testing, alpha_coef)
                })

	alpha.predicted_testing = testing[[1]]
	alpha_coef = testing[[2]]
        rownames(alpha.predicted_testing) = colnames(geno_testing)
        colnames(alpha.predicted_testing) = "test"
	rownames(alpha_coef) = rownames(geno_tuning)
	colnames(alpha_coef) = "weight"

        ##Calculate encode importance##
        enc_imp <- tryCatch({
                coef_nonzero = as.matrix(abs(alpha_coef[which(alpha_coef[,1] != 0),]))
                encode_snp = encode %>% filter(overlap > 0) %>% select(snp)
                coef_nonzero_encode = as.matrix(coef_nonzero[rownames(coef_nonzero) %in% encode_snp$snp,])
                ifelse( !(is.na((sum(coef_nonzero_encode))/(sum(coef_nonzero)))), (sum(coef_nonzero_encode))/(sum(coef_nonzero)), 0)
                },
                error = function(cond) {
                return(0)
                })

	##Calculate Pearson correlation, zscore, and pval of the validating fold##
        pear_folds_v = ifelse(sd(alpha.predicted_validating) != 0, cor( unlist(alpha.predicted_validating), unlist(expression_validating), method = "pearson"), 0)
	spear_folds_v = ifelse(sd(alpha.predicted_validating) != 0, cor( unlist(alpha.predicted_validating), unlist(expression_validating), method = "spearman"), 0)	

        ##Calculate Pearson correlation, zscore, and pval of the test fold##
        pear_folds_t = ifelse(sd(alpha.predicted_testing) != 0, cor( unlist(alpha.predicted_testing), unlist(expression_testing), method = "pearson"), 0)
	spear_folds_t = ifelse(sd(alpha.predicted_testing) != 0, cor( unlist(alpha.predicted_testing), unlist(expression_testing), method = "spearman"), 0)
        pear_zscore_folds <- atanh(pear_folds_t)*sqrt(length(expression_testing) - 3) #Fisher transformation
	spear_zscore_folds <- atanh(spear_folds_t)*sqrt(length(expression_testing) - 3) #Fisher transformation
		
	list(alpha.predicted_validating, pear_folds_v, spear_folds_v, enc_imp, pear_folds_t, spear_folds_t, pear_zscore_folds, spear_zscore_folds)
}

linreg <- function(geno_training, expression_training, geno_validating, expression_validating, geno_tuning, expression_tuning, geno_testing, expression_testing) {

	##Fitting glm to training set and predict in validating set
        x = as.data.frame(t(geno_training))
        colnames(x)[1] = "x"
        y = as.data.frame(t(expression_training))
        colnames(y)[1] = "y"
        df = as.data.frame(cbind(x, y))
        glm.fit = glm(y ~ x, data = df)
        x = as.data.frame(t(geno_validating))
        colnames(x)[1] = "x"
        glm.predicted_validating = as.matrix(predict(glm.fit, newdata = x, interval = "confidence"))

        ##Fitting glm to tuning set and predict in testing set
        x = as.data.frame(t(geno_tuning))
        colnames(x)[1] = "x"
        y = as.data.frame(t(expression_tuning))
        colnames(y)[1] = "y"
        df = as.data.frame(cbind(x, y))
        glm.fit = glm(y ~ x, data = df)
        x = as.data.frame(t(geno_testing))
        colnames(x)[1] = "x"
        glm.predicted_testing = as.matrix(predict(glm.fit, newdata = x, interval = "confidence"))

        ##Calculate encode importance##
        encode_snp = encode %>% filter(overlap > 0) %>% select(snp)
        enc_imp = ifelse( length(intersect(rownames(geno_subset), encode_snp$snp)) == 1 , 1, 0)

	##Calculate Pearson correlation, zscore, and pval of the test fold##
	pear_folds_v = ifelse(sd(glm.predicted_validating) != 0, cor( unlist(glm.predicted_validating), unlist(expression_validating), method = "pearson"), 0)
        spear_folds_v = ifelse(sd(glm.predicted_validating) != 0, cor( unlist(glm.predicted_validating), unlist(expression_validating), method = "spearman"), 0)

        ##Calculate Pearson correlation, zscore, and pval of the test fold##
	pear_folds_t = ifelse(sd(glm.predicted_testing) != 0, cor( unlist(glm.predicted_testing), unlist(expression_testing), method = "pearson"), 0)
        spear_folds_t = ifelse(sd(glm.predicted_testing) != 0, cor( unlist(glm.predicted_testing), unlist(expression_testing), method = "spearman"), 0)

        pear_zscore_folds <- atanh(pear_folds_t)*sqrt(length(expression_testing) - 3) #Fisher transformation
	spear_zscore_folds <- atanh(spear_folds_t)*sqrt(length(expression_testing) - 3) #Fisher transformation

	list(glm.predicted_validating, pear_folds_v, spear_folds_v, enc_imp, pear_folds_t, spear_folds_t, pear_zscore_folds, spear_zscore_folds)
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
total_file_num = as.numeric(args[2])
file_num = as.numeric(args[3])
tissue = args[4]
fold = as.numeric(args[5])
dir_path = normalizePath(args[6])
method = args[7]

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

##Create map file##
if (method == "3D" || method == "3d"){
	func_type = args[8]
	func_path = normalizePath(args[9])
	func = read.table(func_path, header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE) %>% filter(chromosome == chr_num)

	if (func_type == "sc" || func_type == "SC"){
		sc_list = func %>% select(sc) %>% distinct(sc) %>% arrange(sc)
		map = data.frame()
		for (i in sc_list$sc){
			map_temp = sc_process(geno, expression, func %>% filter(sc == i), i)
			map = rbind(map, map_temp)
		}
	} else if (func_type == "pchic" || func_type == "pcHiC"){
		func = func %>% select(-c("length", "MinusLog10Pval"))
		map = pchic_process(geno, expression, func)
	} else {
		func$func_id = paste(func$chromosome, func$start, func$end, sep="_")
        	func = func %>% remove_rownames %>% column_to_rownames(var="func_id")
		map = func_process(geno, expression, func)
	}

} else if (method == "Constant" || method == "constant" || method == "con"){
	func_type = as.numeric(args[8])
	map = constant_process(geno, expression, chr_num, func_type)	
}

##Split files##
map_idx = file_helper(nrow(map), total_file_num)
map = map[map_idx[file_num,1]:map_idx[file_num,2],]

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

	##Penalty3 is based on cut off of number of annotation and set fixed penalty at	2/6##
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

##Update the gene list in expression dataframe according to map file##
geno = geno %>% select(-c("chromosome", "start", "end"))
expression = expression[rownames(expression) %in% map$gene_id,] %>% select(-c("chromosome", "start", "end"))

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
					###Bootstrapping and GLMNET
#########################################################################################################################
summ_result = data.frame()
for (i in 1:nrow(map)) {
	geno_subset = geno[unlist(strsplit(as.character(map[i,3]), ";")),]

	if (nrow(geno_subset) > 1) {
        
	expression_subset = expression[as.character(map[i,2]),]
	origin_id = as.character(map[i,1])

	for (k in 1:nrow(penalties)){	
		print(k)
		prediction_validating = data.frame()
	        prediction_testing = data.frame()      
		pear_folds_v = rep(0,fold)
        	spear_folds_v = rep(0,fold)
		pear_folds_t = rep(0,fold)
		spear_folds_t = rep(0,fold)
                pear_zscore_folds = rep(0,fold)      
		spear_zscore_folds = rep(0,fold)
		enc_imp = rep(0,fold)
		penalty_k = penalties[k, rownames(geno_subset)]			

		for (f in 1:fold) {
			print(f)
			sample_testing = as.list(cv_fold %>% filter(foldid == f) %>% select(sample))$sample
			sample_tuning = as.list(cv_fold %>% filter(foldid != f) %>% select(sample))$sample
			sample_validating = as.list(cv_fold %>% filter(foldid == (f%%fold)+1) %>% select(sample))$sample
			sample_training = setdiff(sample_overlap, union(sample_validating, sample_testing) )

			##Scale the predictors accordingly##				
			geno_output = geno_scale(geno_subset, sample_training, sample_validating)
                        geno_training = geno_output[[1]]
                        geno_validating = geno_output[[2]]

			geno_output = geno_scale(geno_subset, sample_tuning, sample_testing)	
			geno_tuning = geno_output[[1]]
			geno_testing = geno_output[[2]]

			##Scale the expression accordingly##
			expression_output = expression_scale(expression_subset, sample_training, sample_validating)
                        expression_training = expression_output[[1]]
                        expression_validating = expression_output[[2]]
	
			expression_output = expression_scale(expression_subset, sample_tuning, sample_testing)
			expression_tuning = expression_output[[1]]
			expression_testing = expression_output[[2]] 

			result = elnet(geno_training, expression_training, geno_validating, expression_validating, geno_tuning, expression_tuning, geno_testing, expression_testing, penalty_k)
			
			alpha.predicted_validating = result[[1]]
			colnames(alpha.predicted_validating) = "predicted"
			prediction_validating = rbind(prediction_validating, cbind( t(expression_validating), alpha.predicted_validating ) )
			pear_folds_v[f] = result[[2]]
			spear_folds_v[f] = result[[3]]
			enc_imp[f] = result[[4]]
			pear_folds_t[f] = result[[5]]
			spear_folds_t[f] = result[[6]]
			pear_zscore_folds[f] = result[[7]]
			spear_zscore_folds[f] = result[[8]]
		}

		prediction_validating = prediction_validating[match( colnames(expression_subset), rownames(prediction_validating) ),]
	
		##Calculate validating cvm##
	        cvm_validating = min(sapply(2:ncol(prediction_validating), function(x) { mean( (unlist( prediction_validating[,x] ) - unlist( prediction_validating[,1] ))^2 ) } ))		
				
		##Calculate average of correlation across folds##
		pear_avg_v = mean(pear_folds_v)
		spear_avg_v = mean(spear_folds_v)	
	
		pear_avg_t = mean(pear_folds_t)
		spear_avg_t = mean(spear_folds_t)

		##Combine Z-scores via Stouffer's method##
	        pear_zscore_est = sum(pear_zscore_folds) / sqrt(fold)
        	pear_stouffer_pval <- 2*pnorm(abs(pear_zscore_est), lower.tail = FALSE)
			
		spear_zscore_est = sum(spear_zscore_folds) / sqrt(fold)
		spear_stouffer_pval <- 2*pnorm(abs(spear_zscore_est), lower.tail = FALSE)

		##Calculatre average of encode importance##
		enc_avg = mean(enc_imp)

		els_output = data.frame("gene_id" = rownames(expression_subset), "type" = func_type, "region" = origin_id, 
						"penalty" = k, "cvm_v" = cvm_validating, "pear_avg_v" = pear_avg_v, "spear_avg_v" = spear_avg_v,
						"pear_avg_t" = pear_avg_t, "spear_avg_t" = spear_avg_t,
						"pear_stouffer_pval" = pear_stouffer_pval, "spear_stouffer_pval" = spear_stouffer_pval, 
						"enc_avg" = enc_avg) %>% remove_rownames
		summ_result = rbind(summ_result, els_output)	
		filename_summ_result <- paste(dir_path, "/output/", tissue, "/", tissue, "_chr", chr_num, "_GTEx_prediction_", func_type, "_", total_file_num, ".", file_num, ".txt", sep="")
		write.table(summ_result, filename_summ_result, quote=FALSE, sep = "\t", row.names = FALSE)
	}
	
	} else {
                next
        }

}


dummy="I am done"
filename_dummy = paste(dir_path, "/dummy/", tissue, "/dummy_testing_chr", chr_num, "_", func_type, "_", total_file_num, ".", file_num, ".txt", sep="" )
write.table(dummy, filename_dummy, quote=FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
