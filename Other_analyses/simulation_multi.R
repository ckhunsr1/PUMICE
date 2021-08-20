read_plink_custom <- function(root, impute = c('none', 'avg', 'random')) {
    if(impute == 'random') {
        stop("The 'impute' random option has not been implemented.", call. = FALSE)
    }

    ## structure from https://github.com/gabraham/plink2R/blob/master/plink2R/R/plink2R.R
    proot <- path.expand(root)

    bedfile <- paste(proot, ".bed", sep="")
    famfile <- paste(proot, ".fam", sep="")
    bimfile <- paste(proot, ".bim", sep="")

    ## Could change this code to use data.table
    bim <- read.table(bimfile, header=FALSE, sep="", stringsAsFactors=FALSE)
    fam <- read.table(famfile, header=FALSE, sep="", stringsAsFactors=FALSE)
    ## Set the dimensions
    geno <- BEDMatrix::BEDMatrix(bedfile, n = nrow(fam), p = nrow(bim))

    ## Convert to a matrix
    geno <- as.matrix(geno)
    if(impute == 'avg') {
        ## Check if any are missing
        geno_na <- is.na(geno)
        if(any(geno_na)) {
            means <- colMeans(geno, na.rm = TRUE)
            geno[geno_na] <- rep(means, colSums(geno_na))
        }
    }
    colnames(geno) <- bim[,2]
    rownames(geno) <- paste(fam[,1], fam[, 2], sep=":")

    list(bed=geno, fam=fam, bim=bim)
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

##Function to simulate multi-tissue gene expression##
exp_sim	<- function(geno_ss_epi, geno_ss_noepi, result_share, cor_num, n_causal, n_causal_epi, h2e, h2e_epi, h2e_noepi, rho, seed_w, seed_e){
	set.seed(seed_w)

	##Pick causal variants and their associated weights in causal tissues##
        snps_causal_epi = colnames(geno_ss_epi)[ sample( 1:ncol(geno_ss_epi), n_causal_epi ) ]
	snps_causal_noepi = colnames(geno_ss_noepi)[ sample( 1:ncol(geno_ss_noepi), n_causal_noepi ) ]
        w_causal_epi = rnorm(length(snps_causal_epi), mean=0, sd=1)*sqrt((h2e_epi)/length(snps_causal_epi))
	w_causal_noepi = rnorm(length(snps_causal_noepi), mean=0, sd=1)*sqrt((h2e_noepi)/length(snps_causal_noepi))

	##Predefined list of correlated tissue##
	k_list = sort(sample( setdiff(1:length(tissue_list), which(tissue_list == tissue)), cor_num))
	k_tissue_list = tissue_list[k_list]

	################################################################################################################################################################
	##Pick causal variants and their associated weights in correlated tissues by taking into account of percent shared from real data##
	n_causal_epi_shared = round(n_causal_epi * result_share[, k_tissue_list, drop = FALSE])

	##Create list of rho across different tissues for epi snps##
	rho_causal_epi_list = matrix(rho, n_causal_epi, length(tissue_list) - 1)
	rownames(rho_causal_epi_list) = snps_causal_epi
	colnames(rho_causal_epi_list) = setdiff(tissue_list, tissue)
	
	##For non-correlated tissues##
	tissue_leave = setdiff( setdiff(tissue_list, tissue), k_tissue_list)
	if (length(tissue_leave) > 0){
		for (i in 1:length(tissue_leave)){
			tis = tissue_leave[i]
			rho_causal_epi_list[, which(colnames(rho_causal_epi_list) == tis) ] = 0
		}
	}

	##For correlated tissues##
	if (length(k_tissue_list) > 0 ){
		for (i in 1:length(k_tissue_list)){
			tis = k_tissue_list[i]
			n_notshared = n_causal_epi - n_causal_epi_shared[, tis]
			rho_causal_epi_list[ sample(1:length(snps_causal_epi), n_notshared), which(colnames(rho_causal_epi_list) == tis) ] = 0
		}	
	}

        ##Create effect size for correlated/uncorrelated tissues##
        w_cor_epi = matrix(NA, length(snps_causal_epi), length(tissue_list) - 1)
        rownames(w_cor_epi) = snps_causal_epi
        colnames(w_cor_epi) = setdiff(tissue_list, tissue)
        for (j in 1:length(snps_causal_epi)){
                beta_causal_epi = w_causal_epi[j]

		mean = rep(0, length(tissue_list) - 1 )		
		cov = diag( ((h2e_epi)/length(snps_causal_epi)), nrow = length(tissue_list) - 1, ncol = length(tissue_list) - 1)

		rho_temp = rho_causal_epi_list[snps_causal_epi[j], , drop = FALSE]
		rho_temp = colnames(rho_temp)[which(rho_temp > 0)]

		if (length(rho_temp) > 0){
			for (r in 1:length(rho_temp)){
				mean[which(setdiff(tissue_list, tissue) == rho_temp[r])] = rho*beta_causal_epi
				diag(cov)[which(setdiff(tissue_list, tissue) == rho_temp[r])] = ((h2e_epi)/length(snps_causal_epi))*(1-rho^2)
			}
		}
		beta_cor_epi = mvrnorm(1, mean, cov)
                w_cor_epi[j, ] = beta_cor_epi
		rm(j, mean, cov)
        }

	rm(tissue_leave, tis)
	################################################################################################################################################################
	if (length(snps_causal_noepi) >	0){
		##Pick causal variants and their associated weights in correlated tissues by taking into account of percent shared from real data##
		n_causal_noepi_shared = round( (n_causal - n_causal_epi) * result_share[, k_tissue_list, drop = FALSE] )

		##Create list of rho across different tissues for epi snps##
	        rho_causal_noepi_list = matrix(rho, n_causal_noepi, length(tissue_list) - 1)
	        rownames(rho_causal_noepi_list) = snps_causal_noepi
	        colnames(rho_causal_noepi_list) = setdiff(tissue_list, tissue)

	        ##For non-correlated tissues##
	        tissue_leave = setdiff( setdiff(tissue_list, tissue), k_tissue_list)
	        if (length(tissue_leave) > 0){
			for (i in 1:length(tissue_leave)){
		                tis = tissue_leave[i]
	        	        rho_causal_noepi_list[, which(colnames(rho_causal_noepi_list) == tis) ] = 0
	        	}
		}

		##For correlated tissues##
	        if (length(k_tissue_list) > 0 ){
			for (i in 1:length(k_tissue_list)){
		                tis = k_tissue_list[i]
	        	        n_notshared = (n_causal - n_causal_epi) - n_causal_noepi_shared[, tis]
	                	rho_causal_noepi_list[ sample(1:length(snps_causal_noepi), n_notshared), which(colnames(rho_causal_noepi_list) == tis) ] = 0
	        	}
		}

                ##Create effect size for correlated tissues##
                w_cor_noepi = matrix(NA, length(snps_causal_noepi), length(tissue_list) - 1)
                rownames(w_cor_noepi) = snps_causal_noepi
                colnames(w_cor_noepi) = setdiff(tissue_list, tissue)
                for (j in 1:length(snps_causal_noepi)){
                        beta_causal_noepi = w_causal_noepi[j]
			
			mean = rep(0, length(tissue_list) - 1 )
	                cov = diag( ((h2e_noepi)/length(snps_causal_noepi)), nrow = length(tissue_list) - 1, ncol = length(tissue_list) - 1)

        	        rho_temp = rho_causal_noepi_list[snps_causal_noepi[j], , drop = FALSE]
        	        rho_temp = colnames(rho_temp)[which(rho_temp > 0)]
	
			if (length(rho_temp) > 0){
		                for (r in 1:length(rho_temp)){
		                        mean[which(setdiff(tissue_list, tissue) == rho_temp[r])] = rho*beta_causal_noepi
	        	                diag(cov)[which(setdiff(tissue_list, tissue) == rho_temp[r])] = ((h2e_noepi)/length(snps_causal_noepi))*(1-rho^2)
	                	}
			}

			beta_cor_noepi = mvrnorm(1, mean, cov)
                        w_cor_noepi[j, ] = beta_cor_noepi
			rm(j, mean, cov)
                }
		rm(tissue_leave, tis)
        }
	
	################################################################################################################################################################
	set.seed(seed_e)
        ##Create residuals##
	cor_tissue = diag(x = 1, nrow = length(tissue_list), ncol = length(tissue_list)) ##make independent assumption##
        e = mvrnorm(length(gtex_id), rep(0, length(tissue_list)), diag(x = sqrt(1-h2e), nrow = length(tissue_list), ncol = length(tissue_list))%*% cor_tissue %*%
                                                diag(x = sqrt(1-h2e), nrow = length(tissue_list), ncol = length(tissue_list)))

        ##Simulate expression values##
        geno_sim = cbind( geno_ss_epi[, snps_causal_epi, drop = FALSE], geno_ss_noepi[, snps_causal_noepi, drop = FALSE])
        geno_sim = scale(geno_sim) ##scale genotype data##
	if ( length(snps_causal_noepi) > 0 ) {
                w_sim = cbind( c(w_causal_epi, w_causal_noepi), rbind(w_cor_epi, w_cor_noepi))
        } else {
                w_sim = cbind(w_causal_epi, w_cor_epi)
        }
        colnames(w_sim)[1] = tissue
        exp_sim = (geno_sim %*% w_sim) + e
	rownames(exp_sim) = rownames(geno_ss_epi)

        return(exp_sim)
}

##Function to simulate gene expression##
exp_sim_single <- function(genos_ss_epi, genos_ss_noepi, n_causal, n_causal_epi, h2e, h2e_epi, h2e_noepi, seed_w, seed_e){
        set.seed(seed_w)
        snps_causal_epi = colnames(genos_ss_epi)[ sample( 1:ncol(genos_ss_epi), n_causal_epi ) ]
        snps_causal_noepi = colnames(genos_ss_noepi)[ sample( 1:ncol(genos_ss_noepi), n_causal - n_causal_epi ) ]

        w_causal_epi = rnorm(length(snps_causal_epi), mean=0, sd=1)*sqrt((h2e_epi)/length(snps_causal_epi))
        w_causal_noepi = rnorm(length(snps_causal_noepi), mean=0, sd=1)*sqrt((h2e_noepi)/length(snps_causal_noepi))

        set.seed(seed_e)
        ##Create residuals##
        e = rnorm(nrow(genos_ss_epi), mean=0, sd=sqrt(1-h2e))

        ##Simulate expression values##
        genos_sim = cbind( genos_ss_epi[, snps_causal_epi], genos_ss_noepi[, snps_causal_noepi])
        genos_sim = scale(genos_sim) ##scale genotype data##
        w_sim = matrix( c(w_causal_epi, w_causal_noepi) )
        exp_sim = (genos_sim %*% w_sim) + matrix(e)

        return(exp_sim)
}

phe_sim <- function(exp_sim_ss, h2p, seed_e){
        ##Create beta values##
        beta = sqrt(h2p)

	set.seed(seed_e)
        ##Create residuals##
        e = rnorm(length(exp_sim_ss), mean=0, sd=sqrt(1-h2p))

        ##Simulate phenotype values##
        phe_sim = beta*exp_sim_ss + e

        return(phe_sim)
}

###############################################################################################################################################
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=F)
library(data.table)
library(BEDMatrix)
library(tidyr)
library(MASS)
library(dplyr)
library(IRanges)
library(GenomicRanges)
library(glmnet)
library(tidyverse)
library(genefilter)
library(caret)
library(data.table)

ext = "DGN"
tissue = "Whole_Blood"
chr = 22
type = args[1] ##window, 3d##
grid_num = as.numeric(args[2])
gene = args[3] ##ENSG00000133422##
rho = as.numeric(args[4]) ##correlation coefficient between the causal and other tissues##
cor_num = as.numeric(args[5]) ##number of correlated tissues##

##Import genotype data from GTEx datasets##
geno.file = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/simulation/sim_multi/geno/GTEx_chr", chr, "_final", sep = "")
genos = read_plink_custom(geno.file, impute = 'avg')
genos$bed = scale(genos$bed) ##scale genotype##

##Import expression data##
exp.file = "/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/simulation/sim_multi/GTEx_tissue_count.txt"
exp = read.table(exp.file, header = TRUE) %>% filter( chromosome == chr )

##Calcualte correlation among gene expression levels across tissues##
exp_ss = exp %>% filter(gene_id == gene)
tissue_list = unlist(strsplit(exp_ss$tissue_list, ";"))

exp_total = matrix( NA, nrow(genos$bed), length(tissue_list) )
rownames(exp_total) = rownames(genos$bed)
colnames(exp_total) = tissue_list
gtex_id = c()
for (i in 1:length(tissue_list)){
	print(i)
	tis = tissue_list[i]
	expression_path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/input/expression_processed/", tis, "_normalized_expression_final.txt",  sep="")
        temp = as.data.frame( fread( expression_path, header = TRUE) ) %>% filter(gene_id == gene) %>% select(-c("gene_id", "chromosome", "start", "end"))
	colnames(temp) = paste0(colnames(temp), ":", colnames(temp))
	gtex_id = c(gtex_id, colnames(temp))
	temp_add_name = setdiff(rownames(genos$bed), colnames(temp))
	mean_exp = rowMeans(temp)
	temp_add = as.data.frame(t(rep(mean_exp, length(temp_add_name))))
	colnames(temp_add) = temp_add_name 
	temp = cbind( temp, temp_add )
	temp = temp[, match(rownames(exp_total), colnames(temp))]
	exp_total[,i] = as.numeric(temp)
}
#cor_tissue = cor(exp_total)
gtex_id = unique(gtex_id)

###############################################################################################################################################
##Determine snps-gene pair per 3D genome window##
rownames(exp) = exp$gene_id

snps_total = data.frame( "snp" = colnames(genos$bed) )
rownames(snps_total) = snps_total$snp
snps_total = geno_coord(snps_total)

func_list = c("pcHiC", "TAD", "Domain", "Loop")
map_func = data.frame()
for (i in 1:length(func_list)){
	func_type = func_list[i]
	func_path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/input/3D_genome_processed/", func_type, "_input/", tissue, "_", func_type, ".txt", sep = "")
	func = read.table(func_path, header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE) %>% filter(chromosome == chr)
	if ( func_type == "pcHiC" ){
		func = func %>% select(-c("length", "MinusLog10Pval"))
                map = pchic_process(snps_total, exp, func)
		map = data.frame("func" = func_type, map)
	} else {
		func$func_id = paste(func$chromosome, func$start, func$end, sep="_")
       	        func = func %>% remove_rownames %>% column_to_rownames(var="func_id")
       	        map = func_process(snps_total, exp, func)
		map = data.frame("func" = func_type, map)
	}
	map_func = rbind(map_func, map)
	rm(map)
}
map_func = map_func[map_func$gene_id %in% gene, ]
if (nrow(map_func) > 0){
	map_func$snpcount = sapply(strsplit(as.character(map_func$snp_list), ";"), length)
	map_func$windowsize = sapply(strsplit(map_func$origin_id, "_"), function(x){as.numeric(x[3])}) - sapply(strsplit(map_func$origin_id, "_"), function(x){as.numeric(x[2])})
}

##For constant window##
func_list = c(250, 1000)
map_w = data.frame()
for (i in 1:length(func_list)){
        func_type = func_list[i]
        map = constant_process(snps_total, exp, chr, func_type)
        map = data.frame("func" = func_type, map)
        map_w = rbind(map_w, map)
}
map_w = map_w[map_w$gene_id %in% gene, ]
if (nrow(map_w) > 0){
	map_w$snpcount = sapply(strsplit(as.character(map_w$snp_list), ";"), length)
	map_w$windowsize = sapply(strsplit(map_w$origin_id, "_"), function(x){as.numeric(x[3])}) - sapply(strsplit(map_w$origin_id, "_"), function(x){as.numeric(x[2])})
}

map_total = rbind(map_func, map_w)

genos = as.data.frame(genos$bed)
if ( strsplit(type, "-")[[1]][1] == "3d"){
	snp_inc = unique(strsplit( as.character( (map_func %>% filter(func == strsplit(type, "-")[[1]][2]) %>% filter(origin_id == strsplit(type, "-")[[1]][3]))$snp_list), ";")[[1]])
	genos_ss = genos[, colnames(genos) %in% snp_inc]
} else if (strsplit(type, "-")[[1]][1] == "window"){
	snp_inc = unique(strsplit( as.character( (map_w %>% filter(func == strsplit(type, "-")[[1]][2]) %>% filter(origin_id == strsplit(type, "-")[[1]][3]))$snp_list), ";")[[1]])
        genos_ss = genos[, colnames(genos) %in% snp_inc]
} 

###############################################################################################################################################
##Import SCREEN annotation##
encode_path <- paste("/gpfs/group/dxl46/default/private/poom/GTEx/input/screen/", tissue ,"/GTEx_screen_annotation_chr", chr, ".txt", sep="")
encode <- read.table(encode_path, header = FALSE, na.string ='.', as.is=TRUE,check.names=FALSE)
colnames(encode) <- c("chromosome", "start", "end", "overlap")
encode$chromosome <- chr
encode$chromosome <- as.character(encode$chromosome)
encode$start <- as.character(encode$start)
encode$end <- as.character(encode$end)
encode$code <- paste(encode$chromosome, encode$start, encode$end, sep = ":")
encode <- encode %>% select(code, overlap)
encode = encode %>% filter(overlap > 0)

snps_total = data.frame( "snp" = colnames(genos_ss) )
rownames(snps_total) = colnames(genos_ss)
snps_total = geno_coord(snps_total)

##Convert encode bed position to snp name##
geno_code <- snps_total %>% select(chromosome, start, end)
geno_code$chromosome <- as.character(geno_code$chromosome)
geno_code$start <- as.character(geno_code$start)
geno_code$end <- as.character(geno_code$end)
geno_code$code <- paste(geno_code$chromosome, geno_code$start, geno_code$end, sep = ":")
geno_code$snp <- snps_total$snp

encode <- distinct(merge(x = encode, y = geno_code, by.x = "code", by.y = "code") %>% select(snp, overlap))

###############################################################################################################################################
##Simulation settings##
sim_num = 40 ##number of simulations##

if ( strsplit(type, "-")[[1]][2] %in% c("pcHiC", "TAD", "Domain", "Loop") ){
	sim_grid = read.table( paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/simulation/sim_multi/sim_grid_3d.txt", sep = ""), header = TRUE)
} else {
	sim_grid = read.table( paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/simulation/sim_multi/sim_grid_", strsplit(type, "-")[[1]][2], ".txt", sep = ""), header = TRUE)
}
sim_grid = sim_grid[grid_num, ]

n_causal = sim_grid$n_snps ##number of causal snps out of all cis-snps in the region##
n_causal_epi = sim_grid$n_episnps ##number of causal snps that lie in epigenomic regions##
n_causal_noepi = n_causal - n_causal_epi
epi_factor = sim_grid$epi_factor
h2e = sim_grid$h2e ##proportion of expression variance explained by total causal snps##
h2e_epi = (n_causal_epi * epi_factor * h2e) / ( (n_causal_epi*epi_factor) + (n_causal_noepi) ) ##proportion of expression variance explained by causal snps in epigenomics region##
h2e_noepi = h2e - h2e_epi ##proportion of expression variance explained by causal snps NOT in epigenomics region##
h2p = sim_grid$h2p ##proportion of phenotypic variance explained by simulated gene expression levels##

##Load in % significant shared eqtl across GTEx tissue##
load("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/simulation/sim_multi/gtex_eqtl/per_eqtl_shared.RData") ##result##
result_share = result[tissue, , drop = FALSE]
rm(result)

##Do expression simulations in ext, GTEx, and UKB##
path_list = c( paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/simulation/sim_multi/geno/GTEx_chr", chr, "_final", sep = ""),
	       paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/simulation/sim_multi/geno/", ext, "_chr", chr, "_final", sep = ""),
	       paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/simulation/sim_multi/geno/UKB_chr", chr, "_final_ss", sep = ""))


multiexp_sim_list = list()
exp_sim_list = list()
phe_sim_list = list()
set.seed(123)
seeds1<-sample(10**4,length(path_list))
seeds2<-lapply(seeds1, function(jj){set.seed(jj); sample(10**6,sim_num)})
seeds3<-lapply(seeds1 + 1, function(jj){set.seed(jj); sample(10**6,sim_num)})

#####################################################################################################################################################################################
##Simulation for GTEx##
i = 1
multiexp_sim_list[[i]] = list()
set.seed(seeds1[i])	

geno = read_plink_custom(path_list[i], impute = 'avg')
geno_ss = geno$bed[ , colnames(genos_ss)]
geno_ss = genos_ss[rownames(geno_ss) %in% gtex_id, ]

geno_ss_epi = geno_ss[, colnames(geno_ss) %in% encode$snp, drop = FALSE]	
geno_ss_noepi = geno_ss[, !(colnames(geno_ss) %in% encode$snp), drop = FALSE]

##Simulate gene expression levels##
exp_sim_total = as.data.frame( matrix(NA, nrow(geno_ss), sim_num) )
for (j in 1:sim_num){
	temp = exp_sim(geno_ss_epi, geno_ss_noepi, result_share, cor_num, n_causal, n_causal_epi, h2e, h2e_epi, h2e_noepi, rho, j, seeds2[[i]][j])
	multiexp_sim_list[[i]][[j]] = temp
        exp_sim_total[, j] = temp[,1]
        rm(temp)
}
exp_sim_list[[i]] = exp_sim_total
rownames(exp_sim_list[[i]]) = rownames( geno_ss )

##Simulate phenotypes##
phe_sim_total = as.data.frame( matrix(NA, nrow(geno_ss), sim_num) )
for (k in 1:sim_num){
	temp = phe_sim(exp_sim_total[, sim_num], h2p, seeds3[[i]][k])
        phe_sim_total[, k] = temp
        rm(temp)
}
phe_sim_list[[i]] = phe_sim_total
rownames(phe_sim_list[[i]]) = rownames( geno_ss )

#####################################################################################################################################################################################
##Simulation for ext and UKB##
for (i in 2:length(path_list)){
        set.seed(seeds1[i])

        geno = read_plink_custom(path_list[i], impute = 'avg')
        geno_ss = geno$bed[ , colnames(genos_ss)]

	geno_ss_epi = geno_ss[, colnames(geno_ss) %in% encode$snp, drop = FALSE]
        geno_ss_noepi = geno_ss[, !(colnames(geno_ss) %in% encode$snp), drop = FALSE]

        ##Simulate gene expression levels##
        exp_sim_total = as.data.frame( matrix(NA, nrow(geno_ss), sim_num) )
        for (j in 1:sim_num){
                temp = exp_sim_single(geno_ss_epi, geno_ss_noepi, n_causal, n_causal_epi, h2e, h2e_epi,h2e_noepi, j, seeds2[[i]][j])
                exp_sim_total[, j] = temp[,1]
                rm(temp)
        }
	exp_sim_list[[i]] = exp_sim_total
        rownames(exp_sim_list[[i]]) = rownames( geno_ss )

	##Simulate phenotypes##
        phe_sim_total = as.data.frame( matrix(NA, nrow(geno_ss), sim_num) )
        for (k in 1:sim_num){
                temp = phe_sim(exp_sim_total[, sim_num], h2p, seeds3[[i]][k])
                phe_sim_total[, k] = temp
                rm(temp)
        }
	phe_sim_list[[i]] = phe_sim_total
        rownames(phe_sim_list[[i]]) = rownames( geno_ss )
}

save(map_total, multiexp_sim_list, exp_sim_list, phe_sim_list, 
	file = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/simulation/sim_multi/input/", rho, "/", cor_num, "/", type, "_gridnum", grid_num, "_", gene, "_sim.RData", sep = "") )
