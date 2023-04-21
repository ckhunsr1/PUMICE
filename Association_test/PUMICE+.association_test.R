args = commandArgs(trailingOnly=TRUE)

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

convertToComplement <- function(x){
	bases=c("A","C","G","T")
	xx<-unlist(strsplit(toupper(x),NULL))
	paste(unlist(lapply(xx,function(bbb){
	if(bbb=="A") compString<-"T"
	if(bbb=="C") compString<-"G"
	if(bbb=="G") compString<-"C"
	if(bbb=="T") compString<-"A"
	if(!bbb %in% bases) compString<-"N"
	return(compString)
	})),collapse="")
}

pumice_plus_twas <- function(pumice.weight, utmost.weight, gwas.result) {
    gene.vec <- unique( c(pumice.weight$gene, utmost.weight$gene) );
    twas.z.p <- 0; twas.z.u <- 0; twas.z.minp <-0; twas.z.cauchy <- 0;
    pval.p <- 0; pval.u <- 0; pval.minp <- 0; pval.cauchy <- 0;

    for(ii in 1:length(gene.vec)) {
        cat('Analyzing', gene.vec[ii],'\n');
        a <- Sys.time();

        ##Get weight matrix##
	wgt.pumice = pumice.weight %>% filter(gene == gene.vec[ii])
        wgt.utmost = utmost.weight %>% filter(gene == gene.vec[ii])

        ###############################################################################################################################################################
        ##UTMOST imputation and calculation of statistics##

        ##Initialization##
        stat.u <- NA;
        twas.z.u.temp <- NA;

        ##Impute sumstats##
        # Match up the reference data to predicive SNPs (get rid of snps that are not in LD reference panel)
        m = match( wgt.utmost$snp_id , genos$bim[,2] )
        m.keep = !is.na(m)
        wgt.utmost = wgt.utmost[m.keep, ]
        cur.genos = scale(genos$bed[,m[m.keep]])
        cur.bim = genos$bim[m[m.keep],]

        # Match up the SNPs and the summary stats
        m = match(cur.bim[,2] , gwas.result$snp_id)
        cur.Z = gwas.result$beta[m]/gwas.result$se[m]

	if (nrow(wgt.utmost) > 0 & sum(is.na(cur.Z)) != length(cur.Z)) {

                cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)
		cur.LD <- regMat(cur.LD,0.1)
                cur.miss = is.na(cur.Z)

		if (sum(!cur.miss) > opt$max_nsnp) {
                        set.seed(123)
                        index = sample(which(!cur.miss), opt$max_nsnp)
                        index = setdiff(which(!cur.miss), index)
                        cur.not.miss.downsample = !cur.miss
                        cur.not.miss.downsample[index] = FALSE
                } else {
                        cur.not.miss.downsample = !cur.miss
                }

                if ( mean(cur.miss) > opt$max_impute ) {
                        stat.u = NA
                } else {
			cur.wgt =  cur.LD[cur.miss, cur.not.miss.downsample] %*% solve( cur.LD[cur.not.miss.downsample, cur.not.miss.downsample] + 0.1 * diag(sum(cur.not.miss.downsample)) )
                        cur.impz = cur.wgt %*% cur.Z[cur.not.miss.downsample]
                        cur.r2pred = diag( cur.wgt %*% cur.LD[cur.not.miss.downsample, cur.not.miss.downsample] %*% t(cur.wgt) )
                        cur.Z[cur.miss] = cur.impz / sqrt(cur.r2pred)

                        all.r2pred = rep(1,length(cur.Z))
                        all.r2pred[ cur.miss ] = cur.r2pred

                        if ( sum(is.na(all.r2pred)) != 0 ) {
                                stat.u = NA
                        } else if ( mean( all.r2pred[ wgt.utmost$w != 0 ] ) < opt$min_r2pred ) {
                                stat.u = NA
			} else if (sum(is.na(all.r2pred)) == 0 & mean( all.r2pred[ wgt.utmost$w != 0 ] ) >= opt$min_r2pred ) {
                                ##Correlation matrix##
                                cor.u <- cur.LD

                                ##Add the imputed Z-score into the dataframe##
                                z.u = cur.Z

                                ##Weight vector##
                                w.u <- wgt.utmost$w

                                ##Calculate statistics##
                                twas.z.u.temp = sum(z.u*w.u, na.rm=TRUE)/sqrt(as.numeric(t(w.u)%*%cor.u%*%w.u));
                                stat.u <- twas.z.u.temp^2
                        }
                }

        }
	###############################################################################################################################################################
        ##PUMICE imputation and calculation of statistics##

        ##Initialization##
        stat.p <- NA;
        twas.z.p.temp <- NA;

        ##Impute sumstats##
        # Match up the reference data to predicive SNPs (get rid of snps that are not in LD reference panel)
        m = match( wgt.pumice$snp_id, genos$bim[,2] )
        m.keep = !is.na(m)
        wgt.pumice = wgt.pumice[m.keep, ]
        cur.genos = scale(genos$bed[,m[m.keep]])
        cur.bim = genos$bim[m[m.keep],]

        # Match up the SNPs and the summary stats
        m = match(cur.bim[,2] , gwas.result$snp_id )
        cur.Z = gwas.result$beta[m]/gwas.result$se[m]

	if (nrow(wgt.pumice) > 0 & sum(is.na(cur.Z)) != length(cur.Z)) {

                cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)
		cur.LD <- regMat(cur.LD,0.1);
                cur.miss = is.na(cur.Z)

		if (sum(!cur.miss) > opt$max_nsnp) {
                        set.seed(123)
                        index = sample(which(!cur.miss), opt$max_nsnp)
                        index = setdiff(which(!cur.miss), index)
                        cur.not.miss.downsample = !cur.miss
                        cur.not.miss.downsample[index] = FALSE
                } else {
                        cur.not.miss.downsample = !cur.miss
                }

                if ( mean(cur.miss) > opt$max_impute ) {
                        stat.p = NA
                } else {
			cur.wgt =  cur.LD[cur.miss, cur.not.miss.downsample] %*% solve( cur.LD[cur.not.miss.downsample, cur.not.miss.downsample] + 0.1 * diag(sum(cur.not.miss.downsample)) )
                        cur.impz = cur.wgt %*% cur.Z[cur.not.miss.downsample]
                        cur.r2pred = diag( cur.wgt %*% cur.LD[cur.not.miss.downsample, cur.not.miss.downsample] %*% t(cur.wgt) )
                        cur.Z[cur.miss] = cur.impz / sqrt(cur.r2pred)

			all.r2pred = rep(1,length(cur.Z))
                        all.r2pred[ cur.miss ] = cur.r2pred

                        if ( sum(is.na(all.r2pred)) != 0 ) {
                                stat.p = NA
                        } else if ( mean( all.r2pred[ wgt.pumice$w != 0 ] ) < opt$min_r2pred ) {
                                stat.p = NA
                        } else if (sum(is.na(all.r2pred)) == 0 & mean( all.r2pred[ wgt.pumice$w != 0 ] ) >= opt$min_r2pred) {
                                ##Correlation matrix##
                                cor.p <- cur.LD

                                ##Add the imputed Z-score into the dataframe##
                                z.p = cur.Z

                                ##Weight vector##
                                w.p <- wgt.pumice$w

                                ##Calculate statistics##
				twas.z.p.temp = sum(z.p*w.p, na.rm=TRUE)/sqrt(as.numeric(t(w.p)%*%cor.p%*%w.p));
	                        stat.p <- twas.z.p.temp^2
                        }
                }
        }
	###############################################################################################################################################################
        ##Combining p-values##
        twas.z.u[ii] = twas.z.u.temp;
        twas.z.p[ii] = twas.z.p.temp;
        twas.z.minp[ii] = NA;
        twas.z.cauchy[ii] = NA;
        pval.u[ii] <- pchisq(stat.u,df=1,lower.tail=FALSE);
        pval.p[ii] <- pchisq(stat.p,df=1,lower.tail=FALSE);
        pval.minp[ii] <- NA;
        pval.cauchy[ii] <- NA;
        if(!is.na(stat.u) & !is.na(stat.p)) {

            ##Correlation matrices between SNPs from UTMOST and those from PUMICE##
            ##Impute sumstats##
            m = match( union(wgt.utmost$snp_id, wgt.pumice$snp_id) , genos$bim[,2] )
            m.keep = !is.na(m)
            cur.genos = scale(genos$bed[,m[m.keep]])
            cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)
            cur.LD <- regMat(cur.LD,0.1);
            rownames(cur.LD) = genos$bim[m,2]
            colnames(cur.LD) = genos$bim[m,2]
            cor.u.p <- cur.LD[wgt.utmost$snp_id , wgt.pumice$snp_id]

            stat.max <- max(c(stat.u,stat.p),na.rm=TRUE);
            ##pvt <- function(statistic,mu,sigma,alternative=c('greater','less','two.sided'))
            cov.stat <- diag(2);
            cov.stat[1,1] <-t(w.u)%*%cor.u%*%w.u;
            cov.stat[2,2] <-t(w.p)%*%cor.p%*%w.p;
            cov.stat[1,2] <- t(w.u)%*%cor.u.p%*%(w.p);
            cov.stat[2,1] <- cov.stat[1,2];
            cor.stat <- cov2cor(cov.stat);

            pval.minp[ii] <- pvt(stat.max,rep(0,2),cor.stat,alternative='two.sided');
            pval.cauchy[ii] <- cauchy.p(c(pval.u[ii],pval.p[ii]));
	    
	    if ( abs(twas.z.u[ii]) > abs(twas.z.p[ii]) ) {
		     twas.z.minp[ii] <- sign(twas.z.u[ii])*sqrt(qchisq(pval.minp[ii],df=1,lower.tail=F))
		     twas.z.cauchy[ii] <- sign(twas.z.u[ii])*sqrt(qchisq(pval.cauchy[ii],df=1,lower.tail=F))
	    } else {
		     twas.z.minp[ii] <-	sign(twas.z.p[ii])*sqrt(qchisq(pval.minp[ii],df=1,lower.tail=F))
		     twas.z.cauchy[ii] <- sign(twas.z.p[ii])*sqrt(qchisq(pval.cauchy[ii],df=1,lower.tail=F))
	    }

	}
	if(is.na(stat.u) & !is.na(stat.p)) {
            twas.z.minp[ii] <- twas.z.p[ii]
            twas.z.cauchy[ii] <- twas.z.p[ii]
            pval.minp[ii] <- pval.p[ii];
            pval.cauchy[ii] <- pval.p[ii];
        }
	if(!is.na(stat.u) & is.na(stat.p)) {
            twas.z.minp[ii] <- twas.z.u[ii]
            twas.z.cauchy[ii] <- twas.z.u[ii]
            pval.minp[ii] <- pval.u[ii];
            pval.cauchy[ii] <- pval.u[ii];
        }
	b <- Sys.time();
        cat('Analysis of',gene.vec[ii],'  finishes in',format(b-a),'\n');
        print(c(gene.vec[ii],pval.u[ii],pval.p[ii],pval.cauchy[ii]));
    }
    res.out <- cbind(gene.vec,twas.z.u,pval.u,twas.z.p,pval.p,twas.z.cauchy,pval.cauchy);
    return(list(res.out=res.out));
}


#' combine p-values from Cauchy distribution;
#'
#' @param pval a vector of p-values
#' @param weight a vector of weights, the default is equal weight, and the weights add up to 1
#' @return combined p-values;
#' @export
cauchy.p <- function(pval,weight=NULL) {
    if(is.null(weight)) weight <- rep(1/length(pval),length(pval));
    cauchy.stat <- qcauchy(pval,lower.tail=FALSE);
    cauchy.combined <- sum(cauchy.stat*weight);
    p.value <- pcauchy(cauchy.combined,lower.tail=FALSE);
    return(p.value);
}

##########################################################################################################################################################
library(rareGWAMA)
library(BEDMatrix)
library("optparse")
library(dplyr)
library(RSQLite)
options(stringsAsFactors=F)

##########################################################################################################################################################
                                        ##Input parameters##
##########################################################################################################################################################
option_list = list(
                make_option("--geno", action="store", default=NA, type='character',
                                help="Path to genotype data in PLINK format [required]"),
                make_option("--chr", action="store", default=NA, type='character',
                                help="Chromosome number [required]"),
                make_option("--gwas", action="store", default=NA, type='character',
                                help="Path to GWAS summary statistics [required]"),
                make_option("--pumice_weight", action="store", default=NA, type='character',
                                help="Path to PUMICE db file [required]"),
                make_option("--utmost_weight", action="store", default=NA, type='character',
                                help="Path to UTMOST db file. [required]"),
		make_option("--out", action="store", default=NA, type='character',
                                help="Path to output directory and filename [required]"),
		make_option("--max_impute", action="store", default=0.5 , type='double',
		                help="Maximum fraction of SNPs allowed to be missing per gene (will be imputed using LD). [default: %default]"),			  
		make_option("--max_nsnp", action="store", default=100 , type='double',
                                help="Maximum number of SNPs allowed to use for LD-based imputation per gene. [default: %default]"),
  		make_option("--min_r2pred", action="store", default=0.7 , type='double',
              			help="Minimum average LD-based imputation accuracy allowed for expression weight SNP Z-scores. [default: %default]")
              )
opt = parse_args(OptionParser(option_list=option_list))

#opt = list()
#opt$geno = "/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/sample_data/EUR.chr22.final"
#opt$chr = "22"
#opt$gwas = "/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/sample_data/gwas_ss_chr22.txt"
#opt$pumice_weight = "/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/figure_data/db_create/db/PUMICE_Whole_Blood.db"
#opt$utmost_weight = "/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/figure_data/db_create/db/UTMOST_Whole_Blood.db"
#opt$out = "/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/output/pumice+_output_sle_chr22.txt"
#opt$min_r2pred = 0.7
#opt$max_impute = 0.5
#opt$max_nsnp = 100
##########################################################################################################################################################

refpanel_filename = opt$geno
chr_num = as.numeric(as.character(opt$chr))
gwas_filename = opt$gwas
twas_pumice_filename = opt$pumice_weight
twas_utmost_filename = opt$utmost_weight

##########################################################################################################################################################

##Load in reference panel##
genos = read_plink_custom(refpanel_filename, impute = 'avg')
genos$bim = genos$bim[!duplicated(genos$bim[,4]),] ##Only keep one row for row with duplicated position##
genos$bed = genos$bed[,genos$bim[,2]] ##Update bed file according to new bim file##

##########################################################################################################################################################

##Load in GWAS summary statistics##
gwas.result = as.data.frame( fread(gwas_filename, header = TRUE) ) ##effect allele is ref allele##
colnames(gwas.result) = c("snp_id", "chr", "pos", "ea", "oa", "beta", "se")
gwas.result = gwas.result[!(is.na(gwas.result$beta)), ]

##Subset GWAS sumstat for specificied chromosome##
gwas.result = gwas.result[gwas.result$chr == chr_num, ]

##########################################################################################################################################################

##Load in PUMICE TWAS models from db file##
sqlite.driver <- dbDriver("SQLite")
db_p <- dbConnect(sqlite.driver, dbname = twas_pumice_filename)
dbListTables(db_p)
wgtlist_p <- dbReadTable(db_p, "weight")
genename_p <- dbReadTable(db_p, "molecularfeature")

##Subset wgtlist for specified chromosome##
wgtlist_p <- wgtlist_p[wgtlist_p$chrom == chr_num, ]
wgtlist_p <- merge(wgtlist_p, genename_p %>% select(id, ens_gene_id), by.x = "model_id", by.y = "id" )
wgtlist_p <- wgtlist_p %>% select(ens_gene_id, snp, chrom, pos, effect_allele, alt_allele, effect) ##effect allele is ref allele##
colnames(wgtlist_p) = c("gene", "snp_id", "chr", "pos", "ea", "oa", "w")

##########################################################################################################################################################

##Load in UTMOST TWAS models from db file##
sqlite.driver <- dbDriver("SQLite")
db_u <- dbConnect(sqlite.driver, dbname = twas_utmost_filename)
dbListTables(db_u)
wgtlist_u <- dbReadTable(db_u, "weight")
genename_u <- dbReadTable(db_u, "molecularfeature")

##Subset wgtlist for specified chromosome##
wgtlist_u <- wgtlist_u[wgtlist_u$chrom == chr_num, ]
wgtlist_u <- merge(wgtlist_u, genename_u %>% select(id, ens_gene_id), by.x = "model_id", by.y = "id" )
wgtlist_u <- wgtlist_u %>% select(ens_gene_id, snp, chrom, pos, effect_allele, alt_allele, effect) ##effect allele is ref allele##
colnames(wgtlist_u) = c("gene", "snp_id", "chr", "pos", "ea", "oa", "w")

##########################################################################################################################################################

##Run PUMICE+##
x = pumice_plus_twas( wgtlist_p, wgtlist_u, gwas.result)	
result = as.data.frame(x$res.out)
write.table(result, opt$out, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

##########################################################################################################################################################
