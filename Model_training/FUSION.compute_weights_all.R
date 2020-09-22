args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=F)
library(optparse)
library(plink2R)
library(methods)
library(tidyr)
library(dplyr)
library(IRanges)
library(GenomicRanges)
library(glmnet)
library(tidyverse)
library(genefilter)
library(caret)

##Poom's code##
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

#########################################################################################################################
##Create map file to subset snp##
chr_num = args[1]
tissue = args[2]
total_file_num = as.numeric(args[3])
file_num = as.numeric(args[4])

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

##Process expression input##
expression = read.table(expression_path, header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE) %>% filter(chromosome == chr_num)
names(expression)[5:ncol(expression)] = sample_conversion$sample_ID[match(names(expression)[5:ncol(expression)], sample_conversion$subject_ID)]

##Find sample overlap between genotype and expression data##
sample_overlap = intersect(colnames(geno), colnames(expression))

##Update the sample list in both genotype and expression data##
geno = geno %>% select(sample_overlap)
expression = expression %>% select(gene_id, chromosome, start, end, sample_overlap) %>% remove_rownames %>% column_to_rownames(var="gene_id")

##Filling in SNP coordinates##
geno = geno_coord(geno)

##Create map file##
map = constant_process(geno, expression, chr_num, 1000)

##Create cv fold
fold = 5
set.seed(123)
cv_fold_total = createFolds(sample_conversion$sample_ID, k = fold, list = TRUE, returnTrain = FALSE)
cv_fold = data.frame()
for (f in 1:fold) {
        temp = cbind(sample = intersect( sample_overlap, sample_conversion$sample_ID[cv_fold_total[[f]]] ), foldid = f)
        cv_fold = rbind(cv_fold, temp)
}

rm(geno, expression)
######################################################################################################################################################################################################
opt=list()
opt$models = "blup,bslmm"
opt$PATH_gemma = "/gpfs/group/dxl46/default/private/poom/FUSION/gemma-0.98.1-linux-static"
opt$tmp = "/gpfs/group/dxl46/default/private/poom/FUSION/temp"
opt$noclean = FALSE
opt$bfile = paste("/gpfs/group/dxl46/default/private/poom/GTEx/ref/GTEx_chr", chr_num, "_527_qc2", sep = "")
opt$PATH_plink = "/gpfs/group/dxl46/default/share/Poom/bin/plink"
opt$pheno = paste("/gpfs/group/dxl46/default/private/poom/GTEx/input/expression_processed/", tissue, "_normalized_expression_final.txt", sep = "")
opt$rn = FALSE
opt$covar = FALSE
opt$crossval = 5
opt$verbose = 2
SYS_PRINT = FALSE

models = c("blup", "bslmm", "top1")
M = length(models)


# --- PREDICTION MODELS

# BSLMM
weights.bslmm = function( input , bv_type , snp , out=NA ) {
        if ( is.na(out) ) out = paste("chr", chr_num, ".", total_file_num, ".", file_num, ".cv.BSLMM", sep='')
	
        arg = paste( opt$PATH_gemma , " -miss 1 -maf 0 -r2 1 -rpace 1000 -wpace 1000 -bfile " , input , " -bslmm " , bv_type , 
		     " -outdir ", opt$tmp, " -o " , out , sep='' )
        system( arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
        eff = read.table( paste(opt$tmp, "/", out,".param.txt",sep=''),head=T,as.is=T)
        eff.wgt = rep(NA,length(snp))
        m = match( snp , eff$rs )
        m.keep = !is.na(m)
        m = m[m.keep]
        eff.wgt[m.keep] = (eff$alpha + eff$beta * eff$gamma)[m]
        return( eff.wgt )
}

# Marginal Z-scores (used for top1)
weights.marginal = function( genos , pheno , beta=F ) {
        if ( beta ) eff.wgt = t( genos ) %*% (pheno) / ( nrow(pheno) - 1)
        else eff.wgt = t( genos ) %*% (pheno) / sqrt( nrow(pheno) - 1 )
        return( eff.wgt )
}

# --- CLEANUP
cleanup = function() {
        if ( ! opt$noclean ) {
                arg = paste("rm -f " , paste(opt$tmp, "/chr", chr_num, ".", total_file_num, ".", file_num, ".", sep = "") , "*", sep='')
                system(arg)
        }
}

#########################################################################################################################
fam = read.table(paste(opt$bfile,".fam",sep=''),as.is=T)

# Make/fetch the phenotype file
pheno.file = opt$pheno
pheno_original = read.table(pheno.file, header = TRUE, na.string ='.', as.is=TRUE, check.names=FALSE) %>% filter(chromosome == chr_num)

# Match up data
sample_overlap = intersect(fam[,1], colnames(pheno_original)[5:ncol(pheno_original)])
fam = fam[fam$V1 %in% sample_overlap, ]
pheno_original = pheno_original %>% select(gene_id, sample_overlap)
pheno_original = pheno_original[pheno_original$gene_id %in% map$gene_id, ]
pheno_original = pheno_original %>% remove_rownames %>% column_to_rownames(var="gene_id")

## Splitting jobs
gene_list_idx = file_helper(nrow(pheno_original), total_file_num)
pheno_original = pheno_original[gene_list_idx[file_num,1]:gene_list_idx[file_num,2], ]

iv_blup = data.frame()
iv_bslmm = data.frame()
iv_top = data.frame()
summ_blup = data.frame()
summ_bslmm = data.frame()
summ_top = data.frame()
coef_blup = data.frame()
coef_bslmm = data.frame()
coef_top = data.frame()
for (i in 1:nrow(pheno_original)){
        pheno = as.data.frame(t(pheno_original[i,]))
        pheno$V1 = colnames(pheno_original)
        pheno$V2 = pheno$V1
        colnames(pheno)[1] = "V3"
        pheno = pheno %>% select(V1, V2, V3)
        pheno.file = paste(opt$tmp, "/chr", chr_num, ".", total_file_num, ".", file_num, ".pheno",sep="")
        write.table(pheno,quote=F,row.names=F,col.names=F,file=pheno.file)

        snp_included = map %>% filter(gene_id == rownames(pheno_original)[i])
        snp_included = as.data.frame(unlist(strsplit(snp_included$snp_list, split = ";")))

        snp.file = paste(opt$tmp, "/chr", chr_num, ".", total_file_num, ".", file_num, ".snplist",sep="")
        write.table(snp_included,quote=F,row.names=F,col.names=F,file=snp.file)

        geno.file = paste(opt$tmp, "/chr", chr_num, ".", total_file_num, ".", file_num, sep = "")
        # recode to the intersection of samples and new phenotype
        arg = paste( opt$PATH_plink ," --keep-allele-order --bfile ",opt$bfile," --pheno ",pheno.file," --keep ",pheno.file, " --extract ", snp.file, " --make-bed --out ", geno.file, sep='')
        system(arg , ignore.stdout=FALSE, ignore.stderr=FALSE)

        #########################################################################################################################
        # read in genotypes
        genos = read_plink(geno.file,impute="avg")
        mafs = apply(genos$bed,2,mean)/2
        sds = apply(genos$bed,2,sd)
        # important : genotypes are standardized and scaled here:
        genos$bed = scale(genos$bed)
        pheno = genos$fam[,c(1,2,6)]
        pheno[,3] = scale(pheno[,3])

        # check if any genotypes are NA
        nasnps = apply( is.na(genos$bed) , 2 , sum )
        if ( sum(nasnps) != 0 ) {
                cat( "WARNING :",sum(nasnps != 0),"SNPs could not be scaled and were zeroed out, make sure all SNPs are polymorphic\n" , file=stderr())
                genos$bed[,nasnps != 0] = 0
        }
	N.tot = nrow(genos$bed)
        if ( opt$verbose >= 1 ) cat(nrow(pheno),"phenotyped samples, ",nrow(genos$bed),"genotyped samples, ",ncol(genos$bed)," markers\n")

        # --- CROSSVALIDATION ANALYSES
        set.seed(123)
        cv.performance = matrix(NA,nrow=2,ncol=M)
        rownames(cv.performance) = c("rsq","pval")
        colnames(cv.performance) = models

        if ( opt$crossval <= 1 ) {
                if ( opt$verbose >= 1 ) cat("### Skipping cross-validation\n")
        } else {
                if ( opt$verbose >= 1 ) cat("### Performing",opt$crossval,"fold cross-validation\n")
                cv.all = pheno
                N = nrow(cv.all)

		x = merge(sample_conversion[sample_conversion$subject_ID %in% cv.all$V1, ], cv_fold, by.x = "sample_ID", by.y = "sample")
                x = x[order(x$foldid), ]
                cv.sample = match(x$subject_ID, pheno$V1)
                cv.all = cv.all[ cv.sample , ]
                folds = as.numeric(as.character(x$foldid))
                cv.calls = matrix(NA,nrow=N,ncol=M)

		pear_folds = matrix(0, opt$crossval, M)
	        spear_folds = matrix(0, opt$crossval, M)
		pear_zscore_folds = matrix(0, opt$crossval, M)
        	spear_zscore_folds = matrix(0, opt$crossval, M)
                for ( f in 1:opt$crossval ) {
                        if ( opt$verbose >= 1 ) cat("- Crossval fold",f,"\n")
                        indx = which(folds==f,arr.ind=TRUE)
                        cv.train = cv.all[-indx,]

                        # store intercept
                        intercept = mean( cv.train[,3] )
                        cv.train[,3] = scale(cv.train[,3])

                        # hide current fold
                        cv.file = paste(opt$tmp, "/chr", chr_num, ".", total_file_num, ".", file_num, ".cv",sep='')
                        write.table( cv.train , quote=F , row.names=F , col.names=F , file=paste(cv.file,".keep",sep=''))
                        arg = paste( opt$PATH_plink ," --keep-allele-order --bfile ", geno.file, " --keep ", cv.file, ".keep --out ", cv.file, " --make-bed",sep='')
                        system(arg , ignore.stdout=FALSE,ignore.stderr=FALSE)

			for ( mod in 1:M ) {
                                if ( models[mod] == "blup" ) {
                                        pred.wgt = weights.bslmm( cv.file , bv_type=2 , snp=genos$bim[,2] )
					pred.wgt[ is.na(pred.wgt) ] = 0
					predicted = genos$bed[ cv.sample[ indx ] , ] %*% pred.wgt
					cv.calls[ indx , mod ] = predicted
                                
					##Determine correlation and Z-score##
					pear_folds[f,mod] = ifelse( sd(predicted) != 0, cor(predicted, cv.all[ indx, 3], method = "pearson"), 0)
					spear_folds[f,mod] = ifelse( sd(predicted) != 0, cor(predicted, cv.all[ indx, 3], method = "spearman"), 0) 
					pear_zscore_folds[f,mod] = atanh(pear_folds[f,mod])*sqrt(length(indx) - 3)
					spear_zscore_folds[f,mod] = atanh(spear_folds[f,mod])*sqrt(length(indx) - 3)
				}
                                else if ( models[mod] == "bslmm" ) {
                                        pred.wgt = weights.bslmm( cv.file , bv_type=1 , snp=genos$bim[,2] )
					pred.wgt[ is.na(pred.wgt) ] = 0
                                        predicted = genos$bed[ cv.sample[ indx ] , ] %*% pred.wgt
                                        cv.calls[ indx , mod ] = predicted                                    
	
                                        ##Determine correlation and Z-score##
                                        pear_folds[f,mod] = ifelse( sd(predicted) != 0,	cor(predicted, cv.all[ indx, 3], method = "pearson"), 0)
                                        spear_folds[f,mod] = ifelse( sd(predicted) != 0, cor(predicted, cv.all[ indx, 3], method = "spearman"), 0)
                                        pear_zscore_folds[f,mod] = atanh(pear_folds[f,mod])*sqrt(length(indx) - 3)
                                        spear_zscore_folds[f,mod] = atanh(spear_folds[f,mod])*sqrt(length(indx) - 3)
                                }
                                else if ( models[mod] == "top1" ) {
                                        pred.wgt = weights.marginal( genos$bed[ cv.sample[ -indx ],] , as.matrix(cv.train[,3,drop=F]) , beta=T )
                                        pred.wgt[ - which.max( pred.wgt^2 ) ] = 0
					pred.wgt[ is.na(pred.wgt) ] = 0
                                        predicted = genos$bed[ cv.sample[ indx ] , ] %*% pred.wgt
                                        cv.calls[ indx , mod ] = predicted                                    
	
                                        ##Determine correlation and Z-score##
                                        pear_folds[f,mod] = ifelse( sd(predicted) != 0,	cor(predicted, cv.all[ indx, 3], method = "pearson"), 0)
                                        spear_folds[f,mod] = ifelse( sd(predicted) != 0, cor(predicted, cv.all[ indx, 3], method = "spearman"), 0)
                                        pear_zscore_folds[f,mod] = atanh(pear_folds[f,mod])*sqrt(length(indx) - 3)
                                        spear_zscore_folds[f,mod] = atanh(spear_folds[f,mod])*sqrt(length(indx) - 3)
                                }
                        }
                }

		cvm_blup = mean( (cv.calls[,1] - cv.all[,3])^2 )
		cvm_bslmm = mean( (cv.calls[,2] - cv.all[,3])^2 )
		cvm_top = mean( (cv.calls[,3] - cv.all[,3])^2 )

		##Calculate average of correlation across folds##
		pear_avg_blup = mean( pear_folds[, 1] )
		pear_avg_bslmm = mean( pear_folds[, 2] )
		pear_avg_top = mean( pear_folds[, 3] )

		spear_avg_blup = mean( spear_folds[, 1] )
                spear_avg_bslmm = mean( spear_folds[, 2] )
                spear_avg_top = mean( spear_folds[, 3] )

		##Combine Z-scores via Stouffer's method to get p value##
		pear_zscore_blup = sum( pear_zscore_folds[,1] ) / sqrt(opt$crossval)
                pear_stouffer_blup <- 2*pnorm(abs(pear_zscore_blup), lower.tail = FALSE)	
		pear_zscore_bslmm = sum( pear_zscore_folds[,2] ) / sqrt(opt$crossval)
                pear_stouffer_bslmm <- 2*pnorm(abs(pear_zscore_bslmm), lower.tail = FALSE)
		pear_zscore_top = sum( pear_zscore_folds[,3] ) / sqrt(opt$crossval)
                pear_stouffer_top <- 2*pnorm(abs(pear_zscore_top), lower.tail = FALSE)
		
		spear_zscore_blup = sum( spear_zscore_folds[,1] ) / sqrt(opt$crossval)
                spear_stouffer_blup <- 2*pnorm(abs(spear_zscore_blup), lower.tail = FALSE)
                spear_zscore_bslmm = sum( spear_zscore_folds[,2] ) / sqrt(opt$crossval)
                spear_stouffer_bslmm <- 2*pnorm(abs(spear_zscore_bslmm), lower.tail = FALSE)
                spear_zscore_top = sum( spear_zscore_folds[,3] ) / sqrt(opt$crossval)
                spear_stouffer_top <- 2*pnorm(abs(spear_zscore_top), lower.tail = FALSE)

		blup_output = data.frame("tissue" = tissue, "gene_id" = rownames(pheno_original)[i], "type" = 1000, "region" = map[i, 2],
					"cvm_v" = cvm_blup, "pear_avg_t" = pear_avg_blup, "spear_avg_t" = spear_avg_blup, 
					"pear_stouffer_pval" = pear_stouffer_blup, "spear_stouffer_pval" = spear_stouffer_blup) %>% remove_rownames
		iv_blup = rbind(iv_blup, blup_output)

		bslmm_output = data.frame("tissue" = tissue, "gene_id" = rownames(pheno_original)[i], "type" = 1000, "region" = map[i, 2],
                                        "cvm_v" = cvm_bslmm, "pear_avg_t" = pear_avg_bslmm, "spear_avg_t" = spear_avg_bslmm,
                                        "pear_stouffer_pval" = pear_stouffer_bslmm, "spear_stouffer_pval" = spear_stouffer_bslmm) %>% remove_rownames
		iv_bslmm = rbind(iv_bslmm, bslmm_output)

		top_output = data.frame("tissue" = tissue, "gene_id" = rownames(pheno_original)[i], "type" = 1000, "region" = map[i, 2],
                                        "cvm_v" = cvm_top, "pear_avg_t" = pear_avg_top, "spear_avg_t" = spear_avg_top,
                                        "pear_stouffer_pval" = pear_stouffer_top, "spear_stouffer_pval" = spear_stouffer_top) %>% remove_rownames
		iv_top = rbind(iv_top, top_output)

		# --- FULL ANALYSES
		if ( opt$verbose >= 1 ) cat("Computing full-sample weights\n")

		# call models to get weights
		wgt.matrix = matrix(0,nrow=nrow(genos$bim),ncol=M)
		colnames(wgt.matrix) = models
		rownames(wgt.matrix) = genos$bim[,2]
		for ( mod in 1:M ) {
		        if ( models[mod] == "blup" ) {
		                wgt.matrix[,mod] = weights.bslmm( geno.file , bv_type=2 , snp=genos$bim[,2] )
		        }
			else if ( models[mod] == "bslmm" ) {
		                wgt.matrix[,mod] = weights.bslmm( geno.file , bv_type=1 , snp=genos$bim[,2]  )
		        }
			else if ( models[mod] == "top1" ) {
		                wgt.matrix[,mod] = weights.marginal( genos$bed , as.matrix(pheno[,3]) , beta=F )
		        }
		}

		coef_nz_blup = as.data.frame(wgt.matrix[which(wgt.matrix[,1] != 0), 1])
		summ_blup_temp = data.frame("tissue" = tissue, "gene_id" = rownames(pheno_original)[i], "func_type" = 1000, "region" = map[i, 2],
						"fit" = "blup", "input_snps" = nrow(genos$bim), "nonzero_snps" = nrow(coef_nz_blup))
		summ_blup = rbind(summ_blup, summ_blup_temp)		
		coef_blup_temp = data.frame( "tissue" = rep(tissue, nrow(coef_nz_blup)), "gene_id" = rep(rownames(pheno_original)[i], nrow(coef_nz_blup)),
						"snp" = rownames(coef_nz_blup), "weight" = coef_nz_blup[,1])  %>% remove_rownames
		coef_blup = rbind(coef_blup, coef_blup_temp)

		coef_nz_bslmm = as.data.frame(wgt.matrix[which(wgt.matrix[,2] != 0), 2])
                summ_bslmm_temp = data.frame("tissue" = tissue, "gene_id" = rownames(pheno_original)[i], "func_type" = 1000, "region" = map[i, 2],
                                                "fit" = "bslmm", "input_snps" = nrow(genos$bim), "nonzero_snps" = nrow(coef_nz_bslmm))
                summ_bslmm = rbind(summ_bslmm, summ_bslmm_temp)
                coef_bslmm_temp = data.frame( "tissue" = rep(tissue, nrow(coef_nz_bslmm)), "gene_id" = rep(rownames(pheno_original)[i], nrow(coef_nz_bslmm)),
                                                "snp" = rownames(coef_nz_bslmm), "weight" = coef_nz_bslmm[,1])  %>% remove_rownames
                coef_bslmm = rbind(coef_bslmm, coef_bslmm_temp)

		coef_nz_top = as.data.frame(wgt.matrix[which(wgt.matrix[,3] != 0), 3])
                summ_top_temp = data.frame("tissue" = tissue, "gene_id" = rownames(pheno_original)[i], "func_type" = 1000, "region" = map[i, 2],
                                                "fit" = "top1", "input_snps" = nrow(genos$bim), "nonzero_snps" = nrow(coef_nz_top))
                summ_top = rbind(summ_top, summ_top_temp)
                coef_top_temp = data.frame( "tissue" = rep(tissue, nrow(coef_nz_top)), "gene_id" = rep(rownames(pheno_original)[i], nrow(coef_nz_top)),
                                                "snp" = rownames(coef_nz_top), "weight" = coef_nz_top[,1])  %>% remove_rownames
                coef_top = rbind(coef_top, coef_top_temp)

		##Write output for blup##
		filename = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/blup/internal_validation/", tissue, "/", tissue,  "_chr", chr_num, 
				"_GTEx_blup_internal_validation_", total_file_num, ".", file_num, ".txt", sep = "")
		write.table(iv_blup, filename, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

		filename = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/blup/summary/", tissue, "/", tissue,  "_chr", chr_num, 
				"_GTEx_blup_summary_", total_file_num, ".", file_num, ".txt", sep = "")
	        write.table(summ_blup, filename, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

		filename = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/blup/output/", tissue, "/", tissue,  "_chr", chr_num, 
				"_GTEx_blup_coef", total_file_num, ".", file_num, ".txt", sep = "")
		write.table(coef_blup, filename, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

		##Output for bslmm##
		filename = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/bslmm/internal_validation/", tissue, "/", tissue,  "_chr", chr_num,
	                        "_GTEx_bslmm_internal_validation_", total_file_num, ".", file_num, ".txt", sep = "")
	        write.table(iv_bslmm, filename, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

	        filename = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/bslmm/summary/", tissue, "/", tissue,  "_chr", chr_num,
	                        "_GTEx_bslmm_summary_", total_file_num, ".", file_num, ".txt", sep = "")
	        write.table(summ_bslmm, filename, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

	        filename = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/bslmm/output/", tissue, "/", tissue,  "_chr", chr_num,
	                        "_GTEx_bslmm_coef", total_file_num, ".", file_num, ".txt", sep = "")
	        write.table(coef_bslmm, filename, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

		##Output for top1##
		filename = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/top1/internal_validation/", tissue, "/", tissue,  "_chr", chr_num,
	                        "_GTEx_top1_internal_validation_", total_file_num, ".", file_num, ".txt", sep = "")
	        write.table(iv_top, filename, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

	        filename = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/top1/summary/", tissue, "/", tissue,  "_chr", chr_num,
	                        "_GTEx_top1_summary_", total_file_num, ".", file_num, ".txt", sep = "")
	        write.table(summ_top, filename, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

	        filename = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/top1/output/", tissue, "/", tissue,  "_chr", chr_num,
	                        "_GTEx_top1_coef", total_file_num, ".", file_num, ".txt", sep = "")
	        write.table(coef_top, filename, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

	}
}

dummy="I am done"
filename = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/blup/dummy/", tissue, "/dummy_testing_chr", chr_num,
                        "_", total_file_num, ".", file_num, ".txt", sep = "")
write.table(dummy, filename, quote=FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

cleanup()
