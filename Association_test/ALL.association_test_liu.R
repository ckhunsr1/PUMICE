args = commandArgs(trailingOnly=TRUE)

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

twas <- function(gene_overlap, wgtlist, genos, sumstat_new){
FAIL.ctr = 0
## For each wgt file:
for ( w in 1:length(gene_overlap) ) {
        # Load weights
        wgt.matrix = wgtlist %>% filter(gene_id == gene_overlap[w])

        # Match up the reference data to predicive SNPs (get rid of snps that are not in LD reference panel)
        m = match( wgt.matrix[,3] , genos$bim[,2] )
        m.keep = !is.na(m)
        wgt.matrix = wgt.matrix[m.keep, ]

        cur.genos = scale(genos$bed[,m[m.keep]])
        cur.bim = genos$bim[m[m.keep],]

        cur.FAIL = FALSE

        # Match up the SNPs and the summary stats
        m = match(cur.bim[,2] , sumstat_new[,3])
        cur.Z = sumstat_new$Z[m]

        if ( sum(wgt.matrix[, 4] != 0) == 0 ) {
                cat( "WARNING : " , gene_overlap[w] , "had", length(cur.Z) , "overlapping SNPs, but none with non-zero expression weights, try more SNPS or a different model\n")
                cur.FAIL = TRUE
        }

	# Compute LD matrix
        if ( length(cur.Z) == 0 ) {
                cat( "WARNING : " , gene_overlap[w] , " had no overlapping SNPs\n")
                cur.FAIL = TRUE
                out.tbl$NSNP[w] = NA
        } else if ( !cur.FAIL ) {
		cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)
		cur.LD <- regMat(cur.LD,0.1); ##Added by Dr. Liu##
                cur.miss = is.na(cur.Z)
                # Impute missing Z-scores
                if ( sum(cur.miss) != 0 ) {
                        if ( sum(!cur.miss) == 0 ) {
                                cat( "WARNING : " , gene_overlap[w] , "had no overlapping GWAS Z-scores, skipping this gene\n")
                                cur.FAIL = TRUE
                        } else if ( mean(cur.miss) > opt$max_impute ) {
                                cat( "WARNING : " , gene_overlap[w], "had" , sum(cur.miss) , "/" , length(cur.miss) , "non-overlapping GWAS Z-scores, skipping this gene.\n")
                                cur.FAIL = TRUE
                        } else {
                                cur.wgt =  cur.LD[cur.miss,!cur.miss] %*% solve( cur.LD[!cur.miss,!cur.miss] + 0.1 * diag(sum(!cur.miss)) )
                                cur.impz = cur.wgt %*% cur.Z[!cur.miss]
                                cur.r2pred = diag( cur.wgt %*% cur.LD[!cur.miss,!cur.miss] %*% t(cur.wgt) )
                                cur.Z[cur.miss] = cur.impz / sqrt(cur.r2pred)

                                all.r2pred = rep(1,length(cur.Z))
                                all.r2pred[ cur.miss ] = cur.r2pred
                                if ( sum(is.na(all.r2pred)) != 0 ) {
                                        cat( "WARNING : " , gene_overlap[w] , "had missing GWAS Z-scores that could not be imputed, skipping this gene.\n" )
                                        cur.FAIL = TRUE
                                } else if ( mean( all.r2pred[ wgt.matrix[,4] != 0 ] ) < opt$min_r2pred ) {
                                        cat( "WARNING : " , gene_overlap[w] , "had mean GWAS Z-score imputation r2 of" , mean( all.r2pred[ wgt.matrix[,4] != 0 ] ) , "at expression weight SNPs, skipping this gene.\n")
                                        cur.FAIL = TRUE
                                }
                        }
                }
		if ( !cur.FAIL ) {
                        # Compute TWAS Z-score
                        cur.twasz = wgt.matrix[,4] %*% cur.Z
                        cur.twasr2pred = wgt.matrix[,4] %*% cur.LD %*% wgt.matrix[,4]

                        if ( cur.twasr2pred > 0 ) {
                                cur.twas = cur.twasz / sqrt(cur.twasr2pred)
                        } else {
                                cur.FAIL=T
                                cat( "WARNING : " , unlist(wgtlist[w,]) , " had zero predictive accuracy, try a different model.\n")
                        }
                }
        }

	# populate the output
        out.tbl$CHR[w] = chr_num
        out.tbl$ID[w] = gene_overlap[w]
        out.tbl$MODEL[w] = method
        out.tbl$NSNP[w] = nrow(wgtlist %>% filter(gene_id == gene_overlap[w]))

        topgwas = which.max( cur.Z^2 )
        if ( !cur.FAIL && length(topgwas) != 0 && !is.na(topgwas) ) {
                out.tbl$BEST.GWAS.ID[w] = wgt.matrix[ topgwas, 3 ]
                out.tbl$BEST.GWAS.Z[w] = cur.Z[ topgwas ]
        } else {
                out.tbl$BEST.GWAS.ID[w] = NA
                out.tbl$BEST.GWAS.Z[w] = NA
        }

	if ( !cur.FAIL ) {
                out.tbl$NWGT[w] = nrow(wgt.matrix)
                out.tbl$TWAS.Z[w] = cur.twas
                out.tbl$TWAS.P[w] = 2*(pnorm( abs(out.tbl$TWAS.Z[w]) , lower.tail=F))
        } else {
                out.tbl$TWAS.Z[w] = NA
                out.tbl$TWAS.P[w] = NA
        }

	if ( cur.FAIL ) FAIL.ctr = FAIL.ctr + 1

        filename = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/twas_assoc_new/", method, "/", pheno, "/", tissue, "/output/chr", chr_num, "_result.txt", sep = "")
        write.table(out.tbl, filename, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
dummy="I am done"
filename = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/twas_assoc_new/", method, "/", pheno, "/", tissue, "/dummy/chr", chr_num, "_dummy.txt", sep = "")
write.table(dummy, filename, quote=FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

espresso.twas.core.poom <- function(utmost.weight,espresso.weight,gwas.result,vcf.ref.file) {
    gene.vec <- unique(c(utmost.weight$gene,espresso.weight$gene));
    twas.z.minp <-0;
    twas.z.u <- 0;twas.z.e <- 0;twas.z.cauchy <- 0;
    pval.minp <- 0;
    pval.u <- 0;pval.e <- 0;pval.cauchy <- 0;
    for(ii in 1:length(gene.vec)) {
        cat('Analyzing', gene.vec[ii],'\n');
        a <- Sys.time();

        ##Get weight matrix##
        wgt.utmost = utmost.weight %>% filter(gene == gene.vec[ii])
        wgt.espresso = espresso.weight %>% filter(gene == gene.vec[ii])

        ###############################################################################################################################################################
        ##UTMOST imputation and calculation of statistics##

        ##Initialization##
        stat.u <- NA;
        twas.z.u.temp <- NA;

        ##Impute sumstats##
        # Match up the reference data to predicive SNPs (get rid of snps that are not in LD reference panel)
        m = match( wgt.utmost[,8] , genos$bim[,2] )
        m.keep = !is.na(m)
        wgt.utmost = wgt.utmost[m.keep, ]
        cur.genos = scale(genos$bed[,m[m.keep]])
        cur.bim = genos$bim[m[m.keep],]

        # Match up the SNPs and the summary stats
        m = match(cur.bim[,2] , gwas.result[,3])
        cur.Z = gwas.result$beta[m]/gwas.result$se[m]

	if (nrow(wgt.utmost) > 0 & sum(is.na(cur.Z)) != length(cur.Z)) {

                cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)
                cur.LD <- regMat(cur.LD,0.1); ##Added by Dr. Liu##
                cur.miss = is.na(cur.Z)
                if ( mean(cur.miss) > opt$max_impute ) {
                        stat.e = NA
                } else {
			cur.wgt =  cur.LD[cur.miss,!cur.miss] %*% solve( cur.LD[!cur.miss,!cur.miss] + 0.1 * diag(sum(!cur.miss)) )
                        cur.impz = cur.wgt %*% cur.Z[!cur.miss]
                        cur.r2pred = diag( cur.wgt %*% cur.LD[!cur.miss,!cur.miss] %*% t(cur.wgt) )
                        cur.Z[cur.miss] = cur.impz / sqrt(cur.r2pred)
                        all.r2pred = rep(1,length(cur.Z))
                        all.r2pred[ cur.miss ] = cur.r2pred

                        if ( sum(is.na(all.r2pred)) != 0 ) {
                                stat.e = NA
                        } else if ( mean( all.r2pred[ wgt.utmost[, 6] != 0 ] ) < opt$min_r2pred ) {
                                stat.e = NA
			} else if (sum(is.na(all.r2pred)) == 0 & mean( all.r2pred[ wgt.utmost[, 6] != 0 ] ) >= opt$min_r2pred ) {
                                ##Correlation matrix##
                                cor.u <- cur.LD

                                ##Add the imputed Z-score into the dataframe##
                                z.u = cur.Z

                                ##Weight vector##
                                w.u <- wgt.utmost[, 6]

                                ##Calculate statistics##
                                twas.z.u.temp = sum(z.u*w.u, na.rm=TRUE)/sqrt(as.numeric(t(w.u)%*%cor.u%*%w.u));
                                stat.u <- twas.z.u.temp^2
                        }
                }

        }
	###############################################################################################################################################################
        ##Espresso imputation and calculation of statistics##

        ##Initialization##
        stat.e <- NA;
        twas.z.e.temp <- NA;

        ##Impute sumstats##
        # Match up the reference data to predicive SNPs (get rid of snps that are not in LD reference panel)
        m = match( wgt.espresso[,8] , genos$bim[,2] )
        m.keep = !is.na(m)
        wgt.espresso = wgt.espresso[m.keep, ]
        cur.genos = scale(genos$bed[,m[m.keep]])
        cur.bim = genos$bim[m[m.keep],]

        # Match up the SNPs and the summary stats
        m = match(cur.bim[,2] , gwas.result[,3] )
        cur.Z = gwas.result$beta[m]/gwas.result$se[m]

	if (nrow(wgt.espresso) > 0 & sum(is.na(cur.Z)) != length(cur.Z)) {

                cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)
                cur.LD <- regMat(cur.LD,0.1); ##Added by Dr. Liu##
                cur.miss = is.na(cur.Z)

                if ( mean(cur.miss) > opt$max_impute ) {
                        stat.e = NA
                } else {
                        cur.wgt =  cur.LD[cur.miss,!cur.miss] %*% solve( cur.LD[!cur.miss,!cur.miss] + 0.1 * diag(sum(!cur.miss)) )
                        cur.impz = cur.wgt %*% cur.Z[!cur.miss]
                        cur.r2pred = diag( cur.wgt %*% cur.LD[!cur.miss,!cur.miss] %*% t(cur.wgt) )
                        cur.Z[cur.miss] = cur.impz / sqrt(cur.r2pred)
                        all.r2pred = rep(1,length(cur.Z))
                        all.r2pred[ cur.miss ] = cur.r2pred

                        if ( sum(is.na(all.r2pred)) != 0 ) {
                                stat.e = NA
                        } else if ( mean( all.r2pred[ wgt.espresso[, 6] != 0 ] ) < opt$min_r2pred ) {
                                stat.e = NA
                        } else if (sum(is.na(all.r2pred)) == 0 & mean( all.r2pred[ wgt.espresso[, 6] != 0 ] ) >= opt$min_r2pred) {
                                ##Correlation matrix##
                                cor.e <- cur.LD

                                ##Add the imputed Z-score into the dataframe##
                                z.e = cur.Z

                                ##Weight vector##
                                w.e <- wgt.espresso[, 6]

                                ##Calculate statistics##
                                twas.z.e.temp = sum(z.e*w.e, na.rm=TRUE)/sqrt(as.numeric(t(w.e)%*%cor.e%*%w.e));
                                stat.e <- twas.z.e.temp^2
                        }
                }
        }

	###############################################################################################################################################################
        ##Combining p-values##
        twas.z.u[ii] = twas.z.u.temp;
        twas.z.e[ii] = twas.z.e.temp;
        twas.z.minp[ii] = NA;
        twas.z.cauchy[ii] = NA;
        pval.u[ii] <- pchisq(stat.u,df=1,lower.tail=FALSE);
        pval.e[ii] <- pchisq(stat.e,df=1,lower.tail=FALSE);
        pval.minp[ii] <- NA;
        pval.cauchy[ii] <- NA;
        if(!is.na(stat.u) & !is.na(stat.e)) {

            ##Correlation matrices between SNPs from UTMOST and those from Espresso##
            ##Impute sumstats##
            m = match( union(wgt.utmost[,8], wgt.espresso[,8]) , genos$bim[,2] )
            m.keep = !is.na(m)
            cur.genos = scale(genos$bed[,m[m.keep]])
            cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)
            cur.LD <- regMat(cur.LD,0.1); ##Added by Dr. Liu##
            rownames(cur.LD) = genos$bim[m,2]
            colnames(cur.LD) = genos$bim[m,2]
            cor.u.e <- cur.LD[wgt.utmost$snp , wgt.espresso$snp]

            stat.max <- max(c(stat.u,stat.e),na.rm=TRUE);
            ##pvt <- function(statistic,mu,sigma,alternative=c('greater','less','two.sided'))
            cov.stat <- diag(2);
            cov.stat[1,1] <-t(w.u)%*%cor.u%*%w.u;
            cov.stat[2,2] <-t(w.e)%*%cor.e%*%w.e;
            cov.stat[1,2] <- t(w.u)%*%cor.u.e%*%(w.e);
            cov.stat[2,1] <- cov.stat[1,2];
            cor.stat <- cov2cor(cov.stat);

            pval.minp[ii] <- pvt(stat.max,rep(0,2),cor.stat,alternative='two.sided');
            pval.cauchy[ii] <- cauchy.p(c(pval.u[ii],pval.e[ii]));

            #if (sign(twas.z.u[ii])*sign(twas.z.e[ii]) == 1){
            #        twas.z.minp[ii] <- sign(twas.z.u[ii])*sqrt(qchisq(pval.minp[ii],df=1,lower.tail=F))
            #        twas.z.cauchy[ii] <- sign(twas.z.u[ii])*sqrt(qchisq(pval.cauchy[ii],df=1,lower.tail=F))
            #} else {
            #        twas.z.minp[ii] <- NA
            #        twas.z.cauchy[ii] <- NA
            #}
	    
	    if ( abs(twas.z.u[ii]) > abs(twas.z.e[ii]) ) {
		     twas.z.minp[ii] <- sign(twas.z.u[ii])*sqrt(qchisq(pval.minp[ii],df=1,lower.tail=F))
		     twas.z.cauchy[ii] <- sign(twas.z.u[ii])*sqrt(qchisq(pval.cauchy[ii],df=1,lower.tail=F))
	    } else {
		     twas.z.minp[ii] <-	sign(twas.z.e[ii])*sqrt(qchisq(pval.minp[ii],df=1,lower.tail=F))
		     twas.z.cauchy[ii] <- sign(twas.z.e[ii])*sqrt(qchisq(pval.cauchy[ii],df=1,lower.tail=F))
	    }

	}
	if(is.na(stat.u) & !is.na(stat.e)) {
            twas.z.minp[ii] <- twas.z.e[ii]
            twas.z.cauchy[ii] <- twas.z.e[ii]
            pval.minp[ii] <- pval.e[ii];
            pval.cauchy[ii] <- pval.e[ii];
        }
	if(!is.na(stat.u) & is.na(stat.e)) {
            twas.z.minp[ii] <- twas.z.u[ii]
            twas.z.cauchy[ii] <- twas.z.u[ii]
            pval.minp[ii] <- pval.u[ii];
            pval.cauchy[ii] <- pval.u[ii];
        }
	b <- Sys.time();
        cat('Analysis of',gene.vec[ii],'  finishes in',format(b-a),'\n');
        print(c(gene.vec[ii],pval.u[ii],pval.e[ii],pval.minp[ii],pval.cauchy[ii]));
    }
    res.out <- cbind(gene.vec,twas.z.u,pval.u,twas.z.e,pval.e,twas.z.minp,pval.minp,twas.z.cauchy,pval.cauchy);
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
library(fast.lasso)
library(espresso)
library(plink2R)
library("optparse")
library(dplyr)
options(stringsAsFactors=F)

chr_num = args[1]
tissue = args[2]
pheno = args[3]
opt = list()
opt$min_r2pred = 0.7
opt$max_impute = 0.5

##########################################################################################################################################################
# Load in summary stats
sumstat_new = read.table(paste("/gpfs/group/dxl46/default/private/poom/espresso/sumstats_clean_new/", pheno, "/ss_clean_chr", chr_num, ".txt",sep = ""),head=T,as.is=T)
sumstat_new$Z = sumstat_new$beta/sumstat_new$se
sumstat_new = sumstat_new %>% select(-c("p", "n", "chrpos"))

# Load in reference data
genos = read_plink(paste("/gpfs/group/dxl46/default/private/poom/FUSION/input/1000G_hg19_0.05/EUR.chr", chr_num, ".final", sep=""), impute="avg")
genos$bim = genos$bim[!duplicated(genos$bim[,4]),] ##Only keep one row for row with duplicated position##
genos$bed = genos$bed[,genos$bim[,2]] ##Update bed file according to new bim file##

##########################################################################################################################################################
method = "multi-omics"
##Load in internal validation for multi-omics##
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/internal_validation/", tissue, sep = "")
filename = dir(path, pattern = paste("_chr", chr_num, "_", sep = ""))
int = read.table(paste(path, "/", filename, sep =""), header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
multi_sig = (int %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05))$gene_id

# Load in list of weights
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/output/", tissue, sep = "")
filename = dir(path, pattern = paste("_chr", chr_num, "_", sep = ""))
wgtlist = data.frame()
for (i in 1:length(filename)) {
        try({
	temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
        wgtlist = rbind(wgtlist, temp)
	})
}

##Identify gene to be included##
gene_overlap = intersect( unique(wgtlist$gene_id), multi_sig)
wgtlist = wgtlist[wgtlist$gene_id %in% gene_overlap, ]

##Initialize output##
N = length(gene_overlap)
out.tbl = data.frame( "ID" = character(N) , "CHR" = numeric(N) , "BEST.GWAS.ID" = character(N) , "BEST.GWAS.Z" = numeric(N) ,
                      "NSNP" = numeric(N) , "NWGT" = numeric(N) , "MODEL" = character(N) , "TWAS.Z" = numeric(N) , "TWAS.P" = numeric(N) , stringsAsFactors=FALSE )

##Perform TWAS analysis##
twas(gene_overlap, wgtlist, genos, sumstat_new)

##########################################################################################################################################################
method = "predixcan"
##Load in internal validation for predixcan##
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/internal_validation/", tissue, sep = "")
filename = dir(path, pattern = paste("_chr", chr_num, "_", sep = ""))
int = read.table(paste(path, "/", filename, sep =""), header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
predixcan_sig = (int %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05))$gene_id

# Load in list of weights
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/output/", tissue, sep = "")
filename = dir(path, pattern = paste("_chr", chr_num, "_", sep = ""))
wgtlist = data.frame()
for (i in 1:length(filename)) {
        try({
	temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
        wgtlist = rbind(wgtlist, temp)
	})
}

##Identify gene to be included##
gene_overlap = intersect( unique(wgtlist$gene_id), predixcan_sig)
wgtlist = wgtlist[wgtlist$gene_id %in% gene_overlap, ]

##Initialize output##
N = length(gene_overlap)
out.tbl = data.frame( "ID" = character(N) , "CHR" = numeric(N) , "BEST.GWAS.ID" = character(N) , "BEST.GWAS.Z" = numeric(N) ,
                      "NSNP" = numeric(N) , "NWGT" = numeric(N) , "MODEL" = character(N) , "TWAS.Z" = numeric(N) , "TWAS.P" = numeric(N) , stringsAsFactors=FALSE )

##Perform TWAS analysis##
twas(gene_overlap, wgtlist, genos, sumstat_new)

##########################################################################################################################################################
method = "utmost"

# Load in list of weights
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/output/", tissue, sep = "")
filename = dir(path, pattern = paste("_chr", chr_num, "_", sep = ""))
wgtlist = data.frame()
for (i in 1:length(filename)) {
        try({
	temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
        wgtlist = rbind(wgtlist, temp)
	})
}

##Load in internal validation for utmost##
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/internal_validation/", tissue, sep = "")
filename = dir(path, pattern = paste("_chr", chr_num, "_", sep = ""))
int = data.frame()
for (i in 1:length(filename)) {
	try({
        temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
        int = rbind(int, temp)
	})
}
utmost_sig = (int %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05))$gene_id

##Identify gene to be included##
gene_overlap = intersect( unique(wgtlist$gene_id), utmost_sig)
wgtlist = wgtlist[wgtlist$gene_id %in% gene_overlap, ]

##Initialize output##
N = length(gene_overlap)
out.tbl = data.frame( "ID" = character(N) , "CHR" = numeric(N) , "BEST.GWAS.ID" = character(N) , "BEST.GWAS.Z" = numeric(N) ,
                      "NSNP" = numeric(N) , "NWGT" = numeric(N) , "MODEL" = character(N) , "TWAS.Z" = numeric(N) , "TWAS.P" = numeric(N) , stringsAsFactors=FALSE )

##Perform TWAS analysis##
twas(gene_overlap, wgtlist, genos, sumstat_new)

##########################################################################################################################################################

method = "dapg"
##Load in internal validation for dapg##
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/internal_validation/", tissue, sep = "")
filename = dir(path, pattern = paste("_chr", chr_num, "_", sep = ""))
int = read.table(paste(path, "/", filename, sep =""), header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
dapg_sig = (int %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05))$gene_id

##Only for dapg, we also need to filter out genes with no variants with PIP > 0.01##
load(paste("/gpfs/group/dxl46/default/private/poom/dapg/prior/", tissue, "/", "prior_chr", chr_num, ".RData", sep = ""))
final_df = final_df %>% filter(pip > 0.01)
dapg_pip = unique(final_df$gene)
dapg_sig = intersect(dapg_sig, dapg_pip)

# Load in list of weights
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/output/", tissue, sep = "")
filename = dir(path, pattern = paste("_chr", chr_num, "_", sep = ""))
wgtlist = data.frame()
for (i in 1:length(filename)) {
        try({
        temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
        wgtlist = rbind(wgtlist, temp)
        })
}

##Identify gene to be included##
gene_overlap = intersect( unique(wgtlist$gene_id), dapg_sig)
wgtlist = wgtlist[wgtlist$gene_id %in% gene_overlap, ]

##Initialize output##
N = length(gene_overlap)
out.tbl = data.frame( "ID" = character(N) , "CHR" = numeric(N) , "BEST.GWAS.ID" = character(N) , "BEST.GWAS.Z" = numeric(N) ,
                      "NSNP" = numeric(N) , "NWGT" = numeric(N) , "MODEL" = character(N) , "TWAS.Z" = numeric(N) , "TWAS.P" = numeric(N) , stringsAsFactors=FALSE )

##Perform TWAS analysis##
twas(gene_overlap, wgtlist, genos, sumstat_new)

##########################################################################################################################################################
method = "mashr"

##Load in significant model for mashr##
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/sig/", tissue, "/", tissue, "_allchr_GTEx_mashr_sig.txt", sep = "")
mashr_sig = read.table(path, header = TRUE)$sig_gene

# Load in list of weights
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/output/", tissue, sep = "")
filename = dir(path)
wgtlist = read.table(paste(path, "/", filename, sep = ""), header = TRUE)
wgtlist$chr = sapply(strsplit(wgtlist$snp, "_"), function(x){as.character(x[1])})
wgtlist = wgtlist %>% filter(chr == chr_num) %>% select(-c("chr"))
wgtlist$weight = -1*wgtlist$weight

##Identify gene to be included##
gene_overlap = intersect( unique(wgtlist$gene_id), mashr_sig)
wgtlist = wgtlist[wgtlist$gene_id %in% gene_overlap, ]

##Initialize output##
N = length(gene_overlap)
out.tbl = data.frame( "ID" = character(N) , "CHR" = numeric(N) , "BEST.GWAS.ID" = character(N) , "BEST.GWAS.Z" = numeric(N) ,
                      "NSNP" = numeric(N) , "NWGT" = numeric(N) , "MODEL" = character(N) , "TWAS.Z" = numeric(N) , "TWAS.P" = numeric(N) , stringsAsFactors=FALSE )

##Perform TWAS analysis##
twas(gene_overlap, wgtlist, genos, sumstat_new)

##########################################################################################################################################################
##Perform combined method##
##Import reference filename##
vcf.ref.file <- paste("/gpfs/group/dxl46/default/private/dajiang/projects/twas-pred/espresso-test/EUR.chr", 1:22, ".final.vcf.gz", sep ="")

# Load in reference data
genos = read_plink(paste("/gpfs/group/dxl46/default/private/poom/FUSION/input/1000G_hg19_0.05/EUR.chr", chr_num, ".final", sep=""), impute="avg")
genos$bim = genos$bim[!duplicated(genos$bim[,4]),] ##Only keep one row for row with duplicated position##
genos$bed = genos$bed[,genos$bim[,2]] ##Update bed file according to new bim file##

sumstatFile <- paste("/gpfs/group/dxl46/default/private/poom/espresso/sumstats_clean_new/", pheno, "/ss_clean_chr", chr_num, ".txt",sep = "")
gwas.result = read.table(sumstatFile, header = TRUE) %>% select(-c("p", "n"))
gwas.result$beta = -1*gwas.result$beta ##Change the effect allele to alternative allele##
##########################################################################################################################################################
method = "multi-omics"
##Load in internal validation for multi-omics##
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/internal_validation/", tissue, sep = "")
filename = dir(path, pattern = paste("_chr", chr_num, "_", sep = ""))
int = read.table(paste(path, "/", filename, sep =""), header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
multi_sig = (int %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05))$gene_id

# Load in list of weights
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/output/", tissue, sep = "")
filename = dir(path, pattern = paste("_chr", chr_num, "_", sep = ""))
wgtlist = data.frame()
for (i in 1:length(filename)) {
        try({
	temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
        wgtlist = rbind(wgtlist, temp)
        })
}

##Identify gene to be included##
gene_overlap = intersect( unique(wgtlist$gene_id), multi_sig)
wgtlist = wgtlist[wgtlist$gene_id %in% gene_overlap, ]
multi.weight = wgtlist
multi.weight$chr = sapply(strsplit(multi.weight$snp, "_"), function(x){as.character(x[1])})
multi.weight$pos = sapply(strsplit(multi.weight$snp, "_"), function(x){as.character(x[2])})
multi.weight$ref = sapply(strsplit(multi.weight$snp, "_"), function(x){as.character(x[3])})
multi.weight$alt = sapply(strsplit(multi.weight$snp, "_"), function(x){as.character(x[4])})
multi.weight$w = -1*multi.weight$weight
multi.weight$gene = multi.weight$gene_id
multi.weight$chrpos = paste(multi.weight$chr, multi.weight$pos, sep = ":")
multi.weight = multi.weight %>% select(gene, chr, pos, ref, alt, w, chrpos, snp)
espresso.weight = multi.weight
rm(wgtlist, multi.weight)

##########################################################################################################################################################
method = "utmost"
##Load in internal validation for utmost##
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/internal_validation/", tissue, sep = "")
filename = dir(path, pattern = paste("_chr", chr_num, "_", sep = ""))
int = data.frame()
for (i in 1:length(filename)) {
        try({
	temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
        int = rbind(int, temp)
        })
}
utmost_sig = (int %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05))$gene_id

# Load in list of weights
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/output/", tissue, sep = "")
filename = dir(path, pattern = paste("_chr", chr_num, "_", sep = ""))
wgtlist = data.frame()
for (i in 1:length(filename)) {
        try({
        temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
        wgtlist = rbind(wgtlist, temp)
        })
}

##Identify gene to be included##
gene_overlap = intersect( unique(wgtlist$gene_id), multi_sig)
wgtlist = wgtlist[wgtlist$gene_id %in% gene_overlap, ]
utmost.weight = wgtlist
utmost.weight$chr = sapply(strsplit(utmost.weight$snp, "_"), function(x){as.character(x[1])})
utmost.weight$pos = sapply(strsplit(utmost.weight$snp, "_"), function(x){as.character(x[2])})
utmost.weight$ref = sapply(strsplit(utmost.weight$snp, "_"), function(x){as.character(x[3])})
utmost.weight$alt = sapply(strsplit(utmost.weight$snp, "_"), function(x){as.character(x[4])})
utmost.weight$w = -1*utmost.weight$weight
utmost.weight$gene = utmost.weight$gene_id
utmost.weight$chrpos = paste(utmost.weight$chr, utmost.weight$pos, sep = ":")
utmost.weight = utmost.weight %>% select(gene, chr, pos, ref, alt, w, chrpos, snp)

#################################################################################################################################
x = espresso.twas.core.poom(utmost.weight,espresso.weight,gwas.result,vcf.ref.file)

##Write output##
method = "combined"
filename = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/twas_assoc_new/", method, "/", pheno, "/", tissue, "/output/chr", chr_num, "_result.txt", sep = "")
write.table(x$res.out, filename, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

dummy="I am done"
filename = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/twas_assoc_new/", method, "/", pheno, "/", tissue, "/dummy/chr", chr_num, "_dummy.txt", sep = "")
write.table(dummy, filename, quote=FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



