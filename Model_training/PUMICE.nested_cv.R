## Function for splitting jobs
## args
## N: number of total genes
## total_file_num: number of total files
## values
## idx: matrix of fold by 2 with first col being starting index and second col being ending index
file_helper <- function(x,n){
	split(x, cut(seq_along(x), n, labels = FALSE))
}

## Function for filling in missing genotype
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

## Function for mapping snps to gene input in constant method
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
	name = as.data.frame(cbind(range$names, paste(chr, "_", range$start, "_", range$end, sep="")))
	colnames(name) = c("gene_id", "origin_id")

	fo = findOverlaps(query=snp_range, subject=gene_range, type = "within")
	df = as.data.frame(cbind(rownames(expression[subjectHits(fo),]), rownames(geno[queryHits(fo),])))
	colnames(df) = c("gene_id", "snp_id")
        df$gene_id = sapply(strsplit(df$gene_id, "\\."), function(x){as.character(x[1])})
	df$snp_id = sapply(strsplit(df$snp_id, "\\."), function(x){as.character(x[1])})
	df = df %>% group_by(gene_id) %>% mutate(snp_list = paste0(snp_id, collapse = ";")) %>% select(gene_id, snp_list) %>% distinct(gene_id, .keep_all =TRUE)
	
	merge(x = name, y = df, by.x = "gene_id", by.y = "gene_id") %>% select(origin_id, gene_id, snp_list)

}

## Function for mapping snps to genes in 3D method (e.g. Domain, Loop, TAD)
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

## Function for mapping snps to genes in 3D method pchic
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

## Function for scaling genotype dataframe to mean of zero and sd of one
## tuning is interchangable with training, and testing is interchangable with validating##
geno_scale <- function(geno_subset, sample_tuning, sample_testing){
	##Normalized genotype according to the tuning set##
        geno_tuning = geno_subset %>% select(all_of(sample_tuning))
	geno_tuning = geno_fill(geno_tuning)
	center_x_tune = rowMeans(geno_tuning, na.rm = T)
        sd_x_tune = apply(geno_tuning, 1, sd, na.rm = T)
        geno_tuning = as.data.frame(t(scale(t(geno_tuning), center = center_x_tune, scale = sd_x_tune)))
	geno_tuning = geno_tuning[which(!(rowSums(is.na(geno_tuning)) == ncol(geno_tuning))),]
	
	geno_testing = geno_subset %>% select(all_of(sample_testing))
	geno_testing = geno_testing[rownames(geno_testing) %in% rownames(geno_tuning),]
	center_x_tune = center_x_tune[rownames(geno_tuning)]
	sd_x_tune = sd_x_tune[rownames(geno_tuning)]
	geno_testing = as.data.frame(t(scale(t(data.matrix(geno_testing)), center = center_x_tune, scale = sd_x_tune)))
	geno_testing[is.na(geno_testing)] <- 0

	list(geno_tuning, geno_testing)
}

## Function for scaling expression dataframe to mean of zero
## tuning is interchangable with training, and testing is interchangable with validating##
expression_scale <- function(expression_subset, sample_tuning, sample_testing){
	##Normalize expression according to the training set##
	expression_tuning = expression_subset %>% select(all_of(sample_tuning))
       	center_y_tune = rowMeans(expression_tuning, na.rm = T)
	expression_tuning = as.data.frame(t(scale(t(expression_tuning), center = center_y_tune, scale = FALSE)))

        expression_testing = expression_subset %>% select(all_of(sample_testing))
	expression_testing = as.data.frame(t(scale(t(expression_testing), center = center_y_tune, scale = FALSE)))

	list(expression_tuning, expression_testing)
}

## Function for performing nested cross-validation
elnet <- function(geno_training, expression_training, fold_training, geno_validating, expression_validating, geno_tuning, expression_tuning, fold_tuning, geno_testing, expression_testing, penalty_k) {

	##Fit best lambda in training set and predict in validating sample
        alpha.predicted_validating <- tryCatch({
		alpha.fit <- cv.glmnet(t(geno_training), t(expression_training), alpha = 0.5, family="gaussian", penalty.factor = penalty_k, standardize = F, nfolds = opt$fold, type.measure='mse', foldid = fold_training)
                predict(alpha.fit, newx = t(geno_validating), s = "lambda.min")
                },
                error = function(cond) {
		as.matrix(rep(mean(unlist(expression_training)), ncol(geno_validating)))
                })
        rownames(alpha.predicted_validating) = colnames(geno_validating)

        ##Fit tuning set and predict in testing sample
        testing <- tryCatch({
                alpha.fit <- cv.glmnet(t(geno_tuning), t(expression_tuning), alpha = 0.5, family="gaussian", penalty.factor = penalty_k, standardize = F, nfolds = opt$fold, type.measure='mse', foldid = fold_tuning)
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
                encode_snp = geno_bed %>% filter(overlap > 0) %>% select(snpid)
                coef_nonzero_encode = as.matrix(coef_nonzero[rownames(coef_nonzero) %in% encode_snp$snpid,])
                ifelse( !(is.na((sum(coef_nonzero_encode))/(sum(coef_nonzero)))), (sum(coef_nonzero_encode))/(sum(coef_nonzero)), 0)
                },
                error = function(cond) {
                return(0)
                })

	##Calculate Pearson correlation, zscore, and pval of the validating fold##
	R2_v = calc_R2(unlist(expression_validating), unlist(alpha.predicted_validating))
        pear_folds_v = ifelse(sd(alpha.predicted_validating) != 0, cor( unlist(alpha.predicted_validating), unlist(expression_validating), method = "pearson"), 0)
	spear_folds_v = ifelse(sd(alpha.predicted_validating) != 0, cor( unlist(alpha.predicted_validating), unlist(expression_validating), method = "spearman"), 0)	

        ##Calculate Pearson correlation, zscore, and pval of the test fold##
	R2_t = calc_R2(unlist(expression_testing), unlist(alpha.predicted_testing))
        pear_folds_t = ifelse(sd(alpha.predicted_testing) != 0, cor( unlist(alpha.predicted_testing), unlist(expression_testing), method = "pearson"), 0)
	spear_folds_t = ifelse(sd(alpha.predicted_testing) != 0, cor( unlist(alpha.predicted_testing), unlist(expression_testing), method = "spearman"), 0)
        pear_zscore_folds <- atanh(pear_folds_t)*sqrt(length(expression_testing) - 3) #Fisher transformation
	spear_zscore_folds <- atanh(spear_folds_t)*sqrt(length(expression_testing) - 3) #Fisher transformation
		
	list(alpha.predicted_validating, R2_v, pear_folds_v, spear_folds_v, enc_imp, R2_t, pear_folds_t, spear_folds_t, pear_zscore_folds, spear_zscore_folds)
}

## Function for calculating coefficient of determination
calc_R2 <- function(y, y_pred) {
 	tss <- sum((y-mean(y))**2)
	rss <- sum((y - y_pred)**2)
	1 - rss/tss
}

## Function to clean up intermediate file
cleanup = function() {
        if ( ! opt$noclean ) {
                arg = paste("rm -f ", paste(opt$out, "/genotype_epigenomic_overlap_", opt$type, "_", opt$file_num, ".bed", sep = ""), sep = "")
                system(arg)

                arg = paste("rm -f ", paste(opt$out, "/genotype_temp_", opt$type, "_", opt$file_num, ".bed", sep = ""), sep = "")
                system(arg)
        }
}

#########################################################################################################################
#################################################COMMAD_LINE_INPUT#######################################################
#########################################################################################################################


#########################################################################################################################				
					##Import library##
#########################################################################################################################
options(stringsAsFactors=F)
library(optparse)
library(tidyr)
library(dplyr)
library(IRanges)
library(GenomicRanges)
library(glmnet)
library(tidyverse)
library(genefilter)
library(caret)
library(data.table)

#########################################################################################################################
                                        ##Input parameters##
#########################################################################################################################
option_list = list(
		make_option("--geno", action="store", default=NA, type='character',
				help="Path to genotype data in traw format [required]"),
		make_option("--chr", action="store", default=NA, type='character',
                                help="Chromosome number [required]"),
		make_option("--exp", action="store", default=NA, type='character',
                                help="Path to expression data [required]"),
		make_option("--out", action="store", default=NA, type='character',
				help="Path to output directory where output/temporary files are stored [required]"),
		make_option("--method", action="store", default=NA, type='character',
                                help="Window type to be used for creating models. Options: 3d,constant [required]"),
		make_option("--type", action="store", default=NA, type='character',
                                help="Specific 3D genome windows being used/Specific constant window size being used (in kb). Options: loop, tad, pchic, 250, 1000, etc. [required]"),
		make_option("--window_path", action="store", default=NA, type='character',
				help="Path to 3D genome window file. [only required for method == 3d]"),
		make_option("--bedtools_path", action="store", default=NA, type='character',
                                help="Path to bedtools software. [required]"),
		make_option("--epi_path", action="store", default=NA, type='character',
                                help="Path to epigenomic data in bed file format. [required]"),
		make_option("--fold", action="store", default=5, type='integer',
                                help="Number of folds to be performed for cross-validation [default %default]"),
		make_option("--total_file_num", action="store", default=1, type='integer',
                                help="Number of total jobs to be splitted into. This option is helpful for expression file containing many genes. [default %default]"),
		make_option("--file_num", action="store", default=1, type='integer',
                                help="Job number. [default %default]"),
		make_option("--noclean", action="store_true", default=FALSE,
				help="Do not delete any temporary files (for debugging) [default: %default]")
	      )
opt = parse_args(OptionParser(option_list=option_list))

#opt = list()
#opt$geno = "/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/GEUVADIS_chr22_genotype_subset.traw"
#opt$chr = 22
#opt$exp = "/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/GEUVADIS_chr22_expression_subset.txt"
#opt$out = "/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/output"
#opt$method = "constant"
#opt$type = "250"
#opt$window_path = "/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/LCL_chr22_pchic.txt"
#opt$bedtools_path = "/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/bedtools"
#opt$epi_path = "/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/ENCFF028SGJ_chr22_screen.bed"
#opt$fold = 5
#opt$total_file_num = 10
#opt$file_num = 1
#opt$noclean = FALSE

#########################################################################################################################
                                        ##Processing genotype##
#########################################################################################################################

cat( "UPDATE: Processing genotype data\n" , file=stderr())
geno = as.data.frame(fread(opt$geno, header = TRUE)) %>% select(-c(1,3,4,5,6))
geno = geno %>% remove_rownames %>% column_to_rownames(var="SNP")

##Create subjectID-sample ID conversion dataframe for cross-validation step##
sample_conversion = as.data.frame(cbind( "sample_ID" = seq(1, ncol(geno)), "subject_ID" = colnames(geno)))
names(geno) = sample_conversion$sample_ID[match(names(geno), sample_conversion$subject_ID)]

geno = mutate_all(geno, function(x) as.numeric(as.character(x)))

#########################################################################################################################
                                        ##Processing expression##
#########################################################################################################################

cat( "UPDATE: Processing expression data\n" , file=stderr())
expression = as.data.frame(fread(opt$exp, header = TRUE))
sample_overlap = intersect(sample_conversion$subject_ID, colnames(expression))
expression = expression %>% select(gene_id, chromosome, start, end, all_of(sample_overlap))
names(expression)[5:ncol(expression)] = sample_conversion$sample_ID[match(names(expression)[5:ncol(expression)], sample_conversion$subject_ID)]

#########################################################################################################################
                                        ##Processing epigenomic/3D genomics data##
#########################################################################################################################
##Find sample overlap between genotype and expression data##
sample_overlap = intersect(colnames(geno), colnames(expression))

##Update the sample list in both genotype and expression data##
geno = geno %>% select(all_of(sample_overlap))
expression = expression %>% select(gene_id, chromosome, start, end, all_of(sample_overlap)) %>% remove_rownames %>% column_to_rownames(var="gene_id")

##Filling in SNP coordinates##
geno = geno_coord(geno)

##Filter out indel/cnv##
geno = geno %>% filter(end - start == 1)

##Create map file##
cat( "UPDATE: Processing window", opt$method, "-", opt$type, "\n" , file=stderr())
if (opt$method == "3d"){
	func = as.data.frame(fread(opt$window_path, header = TRUE))

	if (opt$type == "pchic"){
		func = func %>% select(-c("length", "MinusLog10Pval"))
		map = pchic_process(geno, expression, func)
	} else {
		func$func_id = paste(func$chromosome, func$start, func$end, sep="_")
        	func = func %>% remove_rownames %>% column_to_rownames(var="func_id")
		map = func_process(geno, expression, func)
	}
} else if (opt$method == "constant"){
	func_type = as.numeric(as.character(opt$type))
	map = constant_process(geno, expression, opt$chr, func_type)
}

##Split into multiple jobs##
if (opt$total_file_num > 1){
	map_idx = file_helper(1:nrow(map), opt$total_file_num)
	map = map[map_idx[[opt$file_num]], ]
}

##Running bedtool to find overlap between SNP and epigenomics data##
cat( "UPDATE: Processing epigenomic data\n" , file=stderr())
geno_bed = data.frame("chr" = paste("chr", geno$chromosome, sep = ""), "start" = geno$start, "end" = geno$end, "snpid" = rownames(geno))
write.table(geno_bed %>% select(-c("snpid")), paste(opt$out, "/genotype_temp_", opt$type, "_", opt$file_num, ".bed", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
arg = paste( opt$bedtools_path, " intersect -a ", paste(opt$out, "/genotype_temp_", opt$type, "_", opt$file_num, ".bed", sep = ""), " -b ", opt$epi_path,
                " -wa -wb > ", paste(opt$out, "/genotype_epigenomic_overlap_", opt$type, "_", opt$file_num, ".bed", sep = ""), sep = "" )
system(arg)

##Import ENCODE file##
penalties <- tryCatch({
	encode_path <- paste(opt$out, "/genotype_epigenomic_overlap_", opt$type, "_", opt$file_num, ".bed", sep = "")
	encode <- as.data.frame(fread(encode_path, header = FALSE))
	encode$id <- paste(encode$V1, encode$V2, encode$V3, sep = "_")	

	geno_bed$id <- paste(geno_bed$chr, geno_bed$start, geno_bed$end, sep = "_")
	geno_bed$overlap <- 0
	geno_bed[geno_bed$id %in% encode$id, "overlap"] <- 1
	geno_bed <- geno_bed[match( rownames(geno), geno_bed$snpid ),]

	##Penalty1 is based on cut off of number of annotation and no penalization on established predictors##
	penalty1 <- c(rep(1, nrow(geno)))
	penalty1[which(geno_bed$overlap > 0)] <- 0

	##Penalty2 is based on cut off of number of annotation and set fixed penalty at 1/6##
	penalty2 <- c(rep(1, nrow(geno)))
	penalty2[which(geno_bed$overlap > 0)] <- 1/6

	##Penalty3 is based on cut off of number of annotation and set fixed penalty at	2/6##
	penalty3 <- c(rep(1, nrow(geno)))
	penalty3[which(geno_bed$overlap > 0)] <- 2/6

	##Penalty4 is based on cut off of number of annotation and set fixed penalty at 3/6##
        penalty4 <- c(rep(1, nrow(geno)))
        penalty4[which(geno_bed$overlap > 0)] <- 3/6

	##Penalty5 is based on cut off of number of annotation and set fixed penalty at 4/6##
        penalty5 <- c(rep(1, nrow(geno)))
        penalty5[which(geno_bed$overlap > 0)] <- 4/6

	##Penalty6 is based on cut off of number of annotation and set fixed penalty at 5/6##
        penalty6 <- c(rep(1, nrow(geno)))
        penalty6[which(geno_bed$overlap > 0)] <- 5/6

	##Penalty7 is elastic net##
	penalty7 <- c(rep(1, nrow(geno)))
	penalty7[which(geno_bed$overlap > 0)] <- 1

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
                                           ##Create cross-validation folds##
#########################################################################################################################
set.seed(123)
cv_fold_testing = createFolds(sample_conversion$sample_ID, k = opt$fold, list = TRUE)

set.seed(123)
cv_fold_tuning = createFolds(sample_conversion$sample_ID, k = opt$fold, returnTrain = TRUE)

cv_fold_validating = list()
for (f in 1:opt$fold){
	cv_fold_validating[[f]] = cv_fold_testing[[(f%%opt$fold)+1]]
}

cv_fold_training = list()
for (f in 1:opt$fold){
        cv_fold_training[[f]] = setdiff(cv_fold_tuning[[f]], cv_fold_validating[[f]])
}

#########################################################################################################################
					###Perform nested cross-validation##
#########################################################################################################################
cat( "UPDATE: Start running nested cross-validation\n" , file=stderr())

summ_result = data.frame()
for (i in 1:nrow(map)) {
	geno_subset = geno[unlist(strsplit(as.character(map[i,3]), ";")),]

	if (nrow(geno_subset) > 1) {
        
		expression_subset = expression[as.character(map[i,2]),]
		origin_id = as.character(map[i,1])

		cat( "UPDATE: Start running gene", rownames(expression_subset), "\n" , file=stderr())
		for (k in 1:nrow(penalties)){	
			cat( "UPDATE: Start running gene", rownames(expression_subset), rownames(penalties)[k], "\n" , file=stderr())
			prediction_validating = data.frame()
		        prediction_testing = data.frame()
			R2_folds_v = rep(0,opt$fold)
			pear_folds_v = rep(0,opt$fold)
        		spear_folds_v = rep(0,opt$fold)
			R2_folds_t = rep(0,opt$fold)
			pear_folds_t = rep(0,opt$fold)
			spear_folds_t = rep(0,opt$fold)
        	        pear_zscore_folds = rep(0,opt$fold)      
			spear_zscore_folds = rep(0,opt$fold)
			enc_imp = rep(0,opt$fold)
			penalty_k = penalties[k, rownames(geno_subset)]			

			for (f in 1:opt$fold) {

				##Identify individual in each folds##
				sample_testing = intersect(sample_overlap, cv_fold_testing[[f]])
				sample_tuning = intersect(sample_overlap, cv_fold_tuning[[f]])
				sample_validating = intersect(sample_overlap, cv_fold_validating[[f]])
				sample_training = intersect(sample_overlap, cv_fold_training[[f]])

				##Ceate fold-ids for inner cv##
        	                set.seed(123)
        	                fold_tuning = sample(rep(seq(opt$fold), length.out = length(cv_fold_tuning[[f]])))
				fold_tuning = fold_tuning[match(sample_tuning, cv_fold_tuning[[f]])]
        	                set.seed(123)
        	                fold_training = sample(rep(seq(opt$fold), length.out = length(cv_fold_training[[f]])))
				fold_training = fold_training[match(sample_training, cv_fold_training[[f]])]

				##Scale the predictors accordingly##				
				geno_output = geno_scale(geno_subset, as.character(sample_training), as.character(sample_validating))
        	                geno_training = geno_output[[1]]
        	                geno_validating = geno_output[[2]]

				geno_output = geno_scale(geno_subset, as.character(sample_tuning), as.character(sample_testing))	
				geno_tuning = geno_output[[1]]
				geno_testing = geno_output[[2]]

				##Scale the expression accordingly##
				expression_output = expression_scale(expression_subset, as.character(sample_training), as.character(sample_validating))
        	                expression_training = expression_output[[1]]
        	                expression_validating = expression_output[[2]]
	
				expression_output = expression_scale(expression_subset, as.character(sample_tuning), as.character(sample_testing))
				expression_tuning = expression_output[[1]]
				expression_testing = expression_output[[2]] 

				result = elnet(geno_training, expression_training, fold_training, geno_validating, expression_validating, geno_tuning, expression_tuning, fold_tuning, geno_testing, expression_testing, penalty_k)
			
				alpha.predicted_validating = result[[1]]
				colnames(alpha.predicted_validating) = "predicted"
				prediction_validating = rbind(prediction_validating, cbind( t(expression_validating), alpha.predicted_validating ) )
				R2_folds_v[f] = result[[2]]
				pear_folds_v[f] = result[[3]]
				spear_folds_v[f] = result[[4]]
				enc_imp[f] = result[[5]]
				R2_folds_t[f] =	result[[6]]
				pear_folds_t[f] = result[[7]]
				spear_folds_t[f] = result[[8]]
				pear_zscore_folds[f] = result[[9]]
				spear_zscore_folds[f] = result[[10]]
			}

			prediction_validating = prediction_validating[match( colnames(expression_subset), rownames(prediction_validating) ),]
	
			##Calculate validating cvm##
	        	cvm_validating = min(sapply(2:ncol(prediction_validating), function(x) { mean( (unlist( prediction_validating[,x] ) - unlist( prediction_validating[,1] ))^2 ) } ))		
				
			##Calculate average of correlation across folds##
			R2_avg_v = mean(R2_folds_v)
			pear_avg_v = mean(pear_folds_v)
			spear_avg_v = mean(spear_folds_v)	
	
			R2_avg_t = mean(R2_folds_t)
			pear_avg_t = mean(pear_folds_t)
			spear_avg_t = mean(spear_folds_t)

			##Combine Z-scores via Stouffer's method##
	        	pear_zscore_est = sum(pear_zscore_folds) / sqrt(opt$fold)
        		pear_stouffer_pval <- 2*pnorm(abs(pear_zscore_est), lower.tail = FALSE)
			
			spear_zscore_est = sum(spear_zscore_folds) / sqrt(opt$fold)
			spear_stouffer_pval <- 2*pnorm(abs(spear_zscore_est), lower.tail = FALSE)

			##Calculatre average of encode importance##
			enc_avg = mean(enc_imp)

			els_output = data.frame("gene_id" = rownames(expression_subset), "type" = opt$type, "region" = origin_id, 
						"penalty" = k, "cvm_v" = cvm_validating, "R2_avg_v" = R2_avg_v, "pear_avg_v" = pear_avg_v, "spear_avg_v" = spear_avg_v,
						"R2_avg_t" = R2_avg_t, "pear_avg_t" = pear_avg_t, "spear_avg_t" = spear_avg_t,
						"pear_stouffer_pval" = pear_stouffer_pval, "spear_stouffer_pval" = spear_stouffer_pval, 
						"enc_avg" = enc_avg) %>% remove_rownames
			summ_result = rbind(summ_result, els_output)
			filename_summ_result <- paste(opt$out, "/", "result_cv_", opt$type, "_", opt$total_file_num, ".", opt$file_num, ".txt", sep="")
			write.table(summ_result, filename_summ_result, quote=FALSE, sep = "\t", row.names = FALSE)
		}

	} else {
		next
        }
}

dummy="I am done"
filename_dummy = paste(opt$out, "/dummy_cv_", opt$type, "_", opt$total_file_num, ".", opt$file_num, ".txt", sep="" )
write.table(dummy, filename_dummy, quote=FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

##Remove temporary files in the output folder##
cleanup()
