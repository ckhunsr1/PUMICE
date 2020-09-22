args = commandArgs(trailingOnly=TRUE)

#' calculate correlation between twas results and faciliate fine mapping;
espressso.twas.cor <- function(utmost.weight,espresso.weight,gwas.result,vcf.ref.file,gene.locale) {

    #gene.vec <- intersect(unique(c(utmost.weight$gene,espresso.weight$gene)),gene.locale);
    gene.vec <- intersect(unique(intersect(utmost.weight$gene,espresso.weight$gene)),gene.locale);
    cor.cauchy.mat <- matrix(nrow=length(gene.vec),ncol=length(gene.vec));
    cor.utmost.mat <- cor.cauchy.mat;
    cor.espresso.mat <- cor.cauchy.mat;
    twas.z.minp <-0;
    twas.z.u <- 0;twas.z.e <- 0;twas.z.cauchy <- 0;
    pval.minp <- 0;
    pval.u <- 0;pval.e <- 0;pval.cauchy <- 0;
    pval.cu <- 0;pval.ce <- 0;
    dat.list <- list();
    
    for(ii in 1:length(gene.vec)) {
        cat('Analyzing', gene.vec[ii],'\n');
        a <- Sys.time();

        ##Get weight matrix##
        wgt.utmost = utmost.weight %>% filter(gene == gene.vec[ii])
        wgt.espresso = espresso.weight %>% filter(gene == gene.vec[ii])

        ##UTMOST imputation and calculation of statistics##
        ##Initialization##
        stat.u <- NA;
        twas.z.u.temp <- NA;

        ##Impute sumstats##
        ## Match up the reference data to predicive SNPs (get rid of snps that are not in LD reference panel)
        m = match( wgt.utmost[,8] , genos$bim[,2] )
        m.keep = !is.na(m)
        wgt.utmost = wgt.utmost[m.keep, ]
        pos.u <- genos$bim[m[m.keep],2];
        cur.genos = scale(genos$bed[,m[m.keep]])
        cur.bim = genos$bim[m[m.keep],]

	# Match up the SNPs and the summary stats
        m = match(cur.bim[,2] , gwas.result[,3])
        cur.Z = gwas.result$beta[m]/gwas.result$se[m]
	cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)
	if ( nrow(cur.LD) > 0 ){ 
		cur.LD <- regMat(cur.LD,0.1); ##Added by Dr. Liu##
	}
	cur.miss = is.na(cur.Z)
        ##print(cur.miss);
        if(sum(cur.miss)<length(cur.miss)) {
            cur.wgt =  cur.LD[cur.miss,!cur.miss] %*% solve( cur.LD[!cur.miss,!cur.miss] + 0.1 * diag(sum(!cur.miss)) )
            cur.impz = cur.wgt %*% cur.Z[!cur.miss]
            cur.r2pred = diag( cur.wgt %*% cur.LD[!cur.miss,!cur.miss] %*% t(cur.wgt) )
            cur.Z[cur.miss] = cur.impz / sqrt(cur.r2pred)
            all.r2pred = rep(1,length(cur.Z))
            all.r2pred[ cur.miss ] = cur.r2pred
            z.u = cur.Z
            w.u <- wgt.utmost[, 6]
            rm(cur.Z)
            twas.z.u.temp = sum(z.u*w.u, na.rm=TRUE)/sqrt(as.numeric(t(w.u)%*%cur.LD%*%w.u));
            stat.u <- twas.z.u.temp^2
        }
        if(sum(cur.miss)==length(cur.miss)) {w.u <- NULL;z.u <- NULL}

	##Espresso imputation and calculation of statistics##
        ##Initialization##
        stat.e <- NA;
        twas.z.e.temp <- NA;

        ##Impute sumstats##
        ## Match up the reference data to predicive SNPs (get rid of snps that are not in LD reference panel)
        m = match( wgt.espresso[,8] , genos$bim[,2] )
        m.keep = !is.na(m)
        wgt.espresso = wgt.espresso[m.keep, ]
        pos.e <- genos$bim[m[m.keep],2];
        cur.genos = scale(genos$bed[,m[m.keep]])
        cur.bim = genos$bim[m[m.keep],]

        # Match up the SNPs and the summary stats
        m = match(cur.bim[,2] , gwas.result[,3] )
        cur.Z = gwas.result$beta[m]/gwas.result$se[m]

	##Impute z-score in the missing snps##
	cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)
	if ( nrow(cur.LD) > 0 ){ 
		cur.LD <- regMat(cur.LD,0.1); ##Added by Dr. Liu## 
	}
        cur.miss = is.na(cur.Z)
        #print(cur.miss);
        if(sum(cur.miss)!=length(cur.miss)) { 
            cur.wgt =  cur.LD[cur.miss,!cur.miss] %*% solve( cur.LD[!cur.miss,!cur.miss] + 0.1 * diag(sum(!cur.miss)) )
            cur.impz = cur.wgt %*% cur.Z[!cur.miss]
            cur.r2pred = diag( cur.wgt %*% cur.LD[!cur.miss,!cur.miss] %*% t(cur.wgt) )
            cur.Z[cur.miss] = cur.impz / sqrt(cur.r2pred)
            all.r2pred = rep(1,length(cur.Z))
            all.r2pred[ cur.miss ] = cur.r2pred
            z.e = cur.Z
            w.e = wgt.espresso[, 6]
            rm(cur.Z)
            twas.z.e.temp = sum(z.e*w.e, na.rm=TRUE)/sqrt(as.numeric(t(w.e)%*%cur.LD%*%w.e));
            stat.e <- twas.z.e.temp^2
        }
        if(sum(cur.miss)==length(cur.miss)) {w.e <- NULL;z.e <- NULL};
        pval.e[ii] <- pchisq(stat.e,df=1,lower.tail=FALSE);
        pval.u[ii] <- pchisq(stat.u,df=1,lower.tail=FALSE);
        pval.ce[ii] <- pval.e[ii]
	pval.cu[ii] <- pval.u[ii]

	if(is.na(pval.u[ii]) & is.na(pval.e[ii])){
		 pval.cauchy[ii] <- NA
	} else {
		if(is.na(pval.e[ii])) pval.ce[ii] <- pval.u[ii];
        	if(is.na(pval.u[ii])) pval.cu[ii] <- pval.e[ii];
        	if(!is.na(pval.cu[ii]) & !is.na(pval.ce[ii])) pval.cauchy[ii] <- cauchy.p(c(pval.ce[ii],pval.cu[ii]));
	}

	dat.list[[ii]] <- list(z.e=z.e,
                               w.e=w.e,
                               z.u=z.u,
                               w.u=w.u,
                               chisq.e=stat.e,
                               twas.z.e=twas.z.e.temp,
                               chisq.u=stat.u,
                               twas.z.u=twas.z.u.temp,
                               pos.e=pos.e,
                               pos.u=pos.u);


    }
    a <- Sys.time();
    cat('calculating covariance matrix\n');
     for(ii in 1:length(gene.vec)) {
        for(jj in 1:length(gene.vec)) {
            if(ii>=jj) {
                ##calculate the correlation of twas.e.ii, twas.e.jj, twas.u.ii, twas.u.jj;
                geno.e.ii <- Matrix(genos$bed[,match(dat.list[[ii]]$pos.e,genos$bim[,2])])
                geno.u.ii <- Matrix(genos$bed[,match(dat.list[[ii]]$pos.u,genos$bim[,2])]);
                geno.e.jj <- Matrix(genos$bed[,match(dat.list[[jj]]$pos.e,genos$bim[,2])]);
                geno.u.jj <- Matrix(genos$bed[,match(dat.list[[jj]]$pos.u,genos$bim[,2])]);
                geno.list <- list(geno.u.ii,geno.e.ii,geno.u.jj,geno.e.jj);
                w.list <- list(dat.list[[ii]]$w.u,dat.list[[ii]]$w.e,
                               dat.list[[jj]]$w.u,dat.list[[jj]]$w.e);
                cov.twas.mat <- matrix(nrow=4,ncol=4);
                for(kk1 in 1:length(geno.list)) {
                    for(kk2 in 1:length(geno.list)) {
                        if(length(w.list[[kk1]])>0 & length(w.list[[kk2]])>0) 
                            cov.twas.mat[kk1,kk2] <- t(w.list[[kk1]])%*%corSparse(geno.list[[kk1]],geno.list[[kk2]])%*%w.list[[kk2]];
                        
                    }
                }
                ##if some entries are missing due to the missing models, we force the statistic to be the other model;
                if(is.na(cov.twas.mat[1,1])) {
                    cov.twas.mat[2,1] <- cov.twas.mat[2,2];
                    cov.twas.mat[1,] <- cov.twas.mat[2,];
                    cov.twas.mat[,1] <- cov.twas.mat[,2];
                }

                if(is.na(cov.twas.mat[2,2])) {
                    cov.twas.mat[1,2] <- cov.twas.mat[1,1];
                    cov.twas.mat[2,] <- cov.twas.mat[1,];
                    cov.twas.mat[,2] <- cov.twas.mat[,1];
                }


                if(is.na(cov.twas.mat[3,3])) {
                    cov.twas.mat[4,3] <- cov.twas.mat[4,4];
                    cov.twas.mat[3,] <- cov.twas.mat[4,];
                    cov.twas.mat[,3] <- cov.twas.mat[,4];
                }

                if(is.na(cov.twas.mat[4,4])) {
                    cov.twas.mat[3,4] <- cov.twas.mat[3,3];
                    cov.twas.mat[4,] <- cov.twas.mat[3,];
                    cov.twas.mat[,4] <- cov.twas.mat[,3];
                }

                
                ##print(cov.twas.mat);
                cov.twas.mat <- rm.na(cov.twas.mat);
                twas.stat.simu <- rmvnorm(1000,sigma=cov.twas.mat);
                pval.twas.mat <- pchisq(twas.stat.simu^2,df=1,lower.tail=F);
                pval.cauchy.mat <- matrix(nrow=nrow(pval.twas.mat),ncol=2);
                ##pval.utmost.mat <- matrix(nrow=nrow(pval.twas.mat),ncol=2);
                pval.utmost.mat <- pchisq((twas.stat.simu[,c(1,3)])^2,df=1,lower.tail=FALSE);
                pval.espresso.mat <- pchisq((twas.stat.simu[,c(2,4)])^2,df=1,lower.tail=FALSE);          
                
                for(ll in 1:nrow(pval.cauchy.mat)) {
                    pval.cauchy.mat[ll,1] <- cauchy.p(pval.twas.mat[ll,1],pval.twas.mat[ll,2]);
                    pval.cauchy.mat[ll,2] <- cauchy.p(pval.twas.mat[ll,3],pval.twas.mat[ll,4]);
                    

                }
                z.cauchy.mat <- qnorm(pval.cauchy.mat,lower.tail=FALSE);
                z.utmost.mat <- qnorm(pval.utmost.mat,lower.tail=FALSE);
                z.espresso.mat <- qnorm(pval.espresso.mat,lower.tail=FALSE);
                
                cor.cauchy.mat[ii,jj] <- cor(z.cauchy.mat)[1,2];
                cor.cauchy.mat[jj,ii] <- cor.cauchy.mat[ii,jj];
                cor.utmost.mat[ii,jj] <- cor(z.utmost.mat)[1,2];
                cor.espresso.mat[ii,jj] <- cor(z.espresso.mat)[1,2];
                cor.utmost.mat[jj,ii] <- cor.utmost.mat[ii,jj];
                cor.espresso.mat[jj,ii] <- cor.espresso.mat[ii,jj];
                
                ## cor.u.ii.u.jj <- corSparse(geno.u.ii,geno.u.jj);
                ## cor.u.ii.e.jj <- corSparse(geno.u.ii,geno.e.jj);
                ## cor.e.ii.u.jj <- corSparse(geno.e.ii,geno.u.jj);
                ## cor.e.ii.e.jj <- corSparse(geno.e.ii,geno.e.jj);
                ## cov.twas.u.ii.u.jj <- t(dat.list[[ii]]$w.u)%*%(cor.u.ii.u.jj)%*%(dat.list[[jj]]$w.u);
                ## cov.twas.u.ii.e.jj <- t(dat.list[[ii]]$w.u)%*%(cor.u.ii.e.jj)%*%(dat.list[[jj]]$w.e);
                ## cov.twas.e.ii.u.jj <- t(dat.list[[ii]]$w.e)%*%(cor.e.ii.u.jj)%*%(dat.list[[jj]]$w.u);
                ## cov.twas.e.ii.e.jj <- t(dat.list[[ii]]$w.u)%*%(cor.e.ii.e.jj)%*%(dat.list[[jj]]$w.e);
            }
	}
    }
    b <- Sys.time();
    cat('finish covariance matrix calcualtion in ',format(b-a),'\n');
    cor.cauchy.mat.raw <- cor.cauchy.mat;
    cor.cauchy.mat <- rm.na(cor.cauchy.mat);
    diag(cor.cauchy.mat) <- 1;
    cor.utmost.mat <- rm.na(cor.utmost.mat)
    cor.espresso.mat <- rm.na(cor.espresso.mat);
    diag(cor.utmost.mat) <- 1;
    diag(cor.espresso.mat) <- 1;
    return(list(cor.cauchy.mat=cor.cauchy.mat,
                pval.e=pval.e,
                pval.u=pval.u,
                gene.vec=gene.vec,
                pval.cauchy=pval.cauchy,
                pval.utmost=pval.u,
                pval.espresso=pval.e,
                cor.utmost.mat=cor.utmost.mat,
                cor.espresso.mat=cor.espresso.mat,
                cor.cauchy.mat.raw=cor.cauchy.mat.raw));
}

#' compute ABF from summary stat;
#' 
#' @param beta genetic effect estimates;
#' @param se the SE of the genetic effect estimates;
#' @return
#'  @export
zscore.abf <- function(z,tau=10000) {
    beta <- z;se <- 1;
    lambda <- sqrt((se)^2/((se)^2+tau^2)) * exp(tau^2*(beta)^2/(2*(se)^2*((se)^2+tau^2)))
    lambda[is.infinite(lambda)] = .Machine$double.xmax
    Lambda <- sum(lambda)
    ppi <- lambda/Lambda;
    return(ppi);
}

zscore2pip <- function(z,cov.z,gene.vec,cs.level=.90,tau=10000) {
    ##z <- qnorm(pval,lower.tail=FALSE);
    ppi <- zscore.abf(z,tau);
    ix.order <- order(ppi,decreasing=TRUE);
    ix.cs95 <- min(which(cumsum(ppi[ix.order])>cs.level));
    ##print(ix.cs95);
    ix.gene <- ix.order[1:ix.cs95];
    gene.cs95 <- gene.vec[ix.gene];
    ix.top <- ix.gene[1];
    ##rownames(ppi) <- gene.vec;
    ppi <- as.data.frame(cbind(gene.vec,ppi));
    return(list(pip=ppi,
                gene.cs=ppi[which(ppi[,1]%in%gene.cs95),],
                gene.top=gene.vec[ix.top]));
}

zscore.finemap <- function(pval,cov.z,gene.vec,max.iter=3,cs.level=.90,alpha=5e-8,tau=10000) {
   z <- qnorm(pval,lower.tail=FALSE);
   no.iter <- 1;
   res <- list();
   res[[1]] <- zscore2pip(z,cov.z,gene.vec,tau=tau);
   rr <- 2;
   z.in <- z;
   p.z.in <- pnorm(z.in,lower.tail=FALSE);
   p.z.in[p.z.in == "0"] = .Machine$double.xmin
   cov.z.in <- cov.z;
   gene.vec.in <- gene.vec;
   while(rr<=max.iter & min(p.z.in)<alpha) {
       ##conditional on top;
       ix.top <- which.max(abs(z.in));
       z.in <- z.in[-ix.top]-matrix(cov.z.in[-ix.top,ix.top],ncol=1)%*%ginv(cov.z.in[ix.top,ix.top])%*%z.in[ix.top];
       cov.z.in <- as.matrix(cov.z.in[-ix.top,-ix.top])-matrix(cov.z.in[-ix.top,ix.top],ncol=length(ix.top))%*%ginv(cov.z.in[ix.top,ix.top])%*%matrix(cov.z.in[ix.top,-ix.top],nrow=length(ix.top));
       p.z.in <- pnorm(z.in/sqrt(abs(diag(cov.z.in))),lower.tail=FALSE);
       ##z.in <- z.rr;
       gene.vec.in <- gene.vec.in[-ix.top];

       if( min(p.z.in)<alpha && !(all(zscore.abf(z.in,tau) == 0)) ) {
           res[[rr]] <- zscore2pip(z.in,cov.z.in,gene.vec.in,cs.level);
       }
       rr <- rr+1;
   }
   return(res);
}

cauchy.p <- function(pval,weight=NULL) {
    if(is.null(weight)) weight <- rep(1/length(pval),length(pval));
    cauchy.stat <- qcauchy(pval,lower.tail=FALSE);
    cauchy.combined <- sum(cauchy.stat*weight);
    p.value <- pcauchy(cauchy.combined,lower.tail=FALSE);
    return(p.value);
}

espresso.output <- function(x, tau_value) {

	pval <- x$pval.cauchy;
        pval[which(is.na(pval))] <- 0.99;
	pval[pval == "1"] <- 0.99;
	pval[pval < .Machine$double.xmin] = .Machine$double.xmin
        cov.z <- x$cor.cauchy.mat;
        gene.vec <- x$gene.vec;
        ##zscore.finemap <- function(pval,cov.z,gene.vec,max.iter=3,cs.level=.90)
        res.finemap.cauchy <- zscore.finemap(pval,cov.z,gene.vec, tau = tau_value);
        ##Output cauchy result##
        cauchy_temp = res.finemap.cauchy[[1]]$pip
        cauchy_temp$in_cred_set = 0
        cauchy_temp[cauchy_temp$gene.vec %in% (res.finemap.cauchy[[1]]$gene.cs)$gene.vec, "in_cred_set"] = 1
        cauchy_temp = cauchy_temp[order(as.numeric(as.character(cauchy_temp$ppi)), decreasing = TRUE),]
        cauchy_temp$region = unique(df_reg$V1)[idx]
        colnames(cauchy_temp) = c("gene", "pip", "in_cred_set", "region")
	rm(pval, cov.z, gene.vec)

        pval <- x$pval.espresso;
        pval[which(is.na(pval))] <- 0.99;
	pval[pval == "1"] <- 0.99;
	pval[pval < .Machine$double.xmin] = .Machine$double.xmin
        cov.z <- x$cor.espresso.mat;
        gene.vec <- x$gene.vec;
        ##zscore.finemap <- function(pval,cov.z,gene.vec,max.iter=3,cs.level=.90)
        res.finemap.espresso <- zscore.finemap(pval,cov.z,gene.vec, tau = tau_value);
        ##Output espresso result##
        espresso_temp = res.finemap.espresso[[1]]$pip
        espresso_temp$in_cred_set = 0
        espresso_temp[espresso_temp$gene.vec %in% (res.finemap.espresso[[1]]$gene.cs)$gene.vec, "in_cred_set"] = 1
        espresso_temp = espresso_temp[order(as.numeric(as.character(espresso_temp$ppi)), decreasing = TRUE),]
        espresso_temp$region = unique(df_reg$V1)[idx]
        colnames(espresso_temp) = c("gene", "pip", "in_cred_set", "region")
	rm(pval, cov.z,	gene.vec)

        pval <- x$pval.utmost;
        pval[which(is.na(pval))] <- 0.99;
	pval[pval == "1"] <- 0.99;
	pval[pval < .Machine$double.xmin] = .Machine$double.xmin
        cov.z <- x$cor.utmost.mat;
        gene.vec <- x$gene.vec;
        ##zscore.finemap <- function(pval,cov.z,gene.vec,max.iter=3,cs.level=.90)
        res.finemap.utmost <- zscore.finemap(pval,cov.z,gene.vec,tau = tau_value);
        ##Output utmost result##
        utmost_temp = res.finemap.utmost[[1]]$pip
        utmost_temp$in_cred_set = 0
        utmost_temp[utmost_temp$gene.vec %in% (res.finemap.utmost[[1]]$gene.cs)$gene.vec, "in_cred_set"] = 1
        utmost_temp  = utmost_temp[order(as.numeric(as.character(utmost_temp$ppi)), decreasing = TRUE),]
        utmost_temp$region = unique(df_reg$V1)[idx]
        colnames(utmost_temp) = c("gene", "pip", "in_cred_set", "region")
	rm(pval, cov.z,	gene.vec)

	temp = list()
	temp[[1]] = cauchy_temp
	temp[[2]] = espresso_temp
	temp[[3]] = utmost_temp
	
	return(temp)
}

############################################################################################################################################################################################
library(rareGWAMA)
library(fast.lasso)
library(espresso)
library(dplyr)
library(plink2R)
library(GenomicRanges)
library(IRanges)
options(stringsAsFactors=F)

chr_num = as.numeric(args[1])
tissue = args[2]
pheno = args[3]
tau_value = as.numeric(args[4])
#chr_num = 22
#tissue = "Adipose_Subcutaneous"
#pheno = "ast"
#tau_value = 2

opt = list()
opt$min_r2pred = 0.7
opt$max_impute = 0.5

##Import reference filename##
vcf.ref.file <- paste("/gpfs/group/dxl46/default/private/dajiang/projects/twas-pred/espresso-test/EUR.chr", 1:22, ".final.vcf.gz", sep ="")

# Load in reference data
genos = read_plink(paste("/gpfs/group/dxl46/default/private/poom/FUSION/input/1000G_hg19_0.05/EUR.chr", chr_num, ".final", sep=""), impute="avg")
genos$bim = genos$bim[!duplicated(genos$bim[,4]),] ##Only keep one row for row with duplicated position##
genos$bed = genos$bed[,genos$bim[,2]] ##Update bed file according to new bim file##

sumstatFile <- paste("/gpfs/group/dxl46/default/private/poom/espresso/sumstats_clean_new/", pheno, "/ss_clean_chr", chr_num, ".txt",sep = "")
gwas.result = read.table(sumstatFile, header = TRUE) %>% select(-c("n"))
gwas.result$beta = -1*gwas.result$beta ##Change the effect allele to alternative allele##

##Load in model counts##
count = read.table("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/iv_all_count.txt", header = TRUE)
count = as.data.frame(cbind("tissue" = count[,1], "count" = apply(count[,2:4], 1, max)))
count$count = as.numeric(as.character(count$count))
#################################################################################################################################
##Load in internal validation for utmost##
method = "utmost"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/internal_validation/", tissue, sep = "")
filename = dir(path, pattern = paste("_chr", chr_num, "_", sep = ""))
int = data.frame()
for (i in 1:length(filename)) {
        temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
        int = rbind(int, temp)
}
utmost_sig = (int %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05))$gene_id

##Load in TWAS result##
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/twas_assoc_new/", method, "/", pheno, "/", tissue, "/output/chr", chr_num, "_result.txt", sep = "")
twas_utmost = read.table(path, header = TRUE) %>% filter(TWAS.P < 0.05/count[which(count$tissue == tissue), 2])
utmost_reg = intersect(utmost_sig, twas_utmost$ID)

##Import UTMOST weight##
utmostFile <- dir(paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/output/", tissue, sep = ""), 
		pattern = paste("_chr", chr_num, "_", sep = ""))
utmost.weight = data.frame()
for (i in 1:length(utmostFile)) {
	temp = read.table(paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/output/", 
				tissue, "/", utmostFile[i], sep = ""), header = TRUE)
	utmost.weight = rbind(utmost.weight, temp)
}
utmost.weight = utmost.weight[utmost.weight$gene_id %in% utmost_sig, ]
utmost.weight$chr = sapply(strsplit(utmost.weight$snp, "_"), function(x){as.character(x[1])})
utmost.weight$pos = sapply(strsplit(utmost.weight$snp, "_"), function(x){as.character(x[2])})
utmost.weight$ref = sapply(strsplit(utmost.weight$snp, "_"), function(x){as.character(x[3])})
utmost.weight$alt = sapply(strsplit(utmost.weight$snp, "_"), function(x){as.character(x[4])})
utmost.weight$w = -1*utmost.weight$weight
utmost.weight$gene = utmost.weight$gene_id
utmost.weight$chrpos = paste(utmost.weight$chr, utmost.weight$pos, sep = ":")
utmost.weight = utmost.weight %>% select(gene, chr, pos, ref, alt, w, chrpos, snp)


#################################################################################################################################
##Load in internal validation for multi-omics##
method = "multi-omics"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/internal_validation/", tissue, sep = "")
filename = dir(path, pattern = paste("_chr", chr_num, "_", sep = ""))
int = data.frame()
for (i in 1:length(filename)) {
        temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
        int = rbind(int, temp)
}
multi_sig = (int %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05))$gene_id

##Load in TWAS result##
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/twas_assoc_new/", method, "/", pheno, "/", tissue, "/output/chr", chr_num, "_result.txt", sep = "")
twas_multi = read.table(path, header = TRUE) %>% filter(TWAS.P < 0.05/count[which(count$tissue == tissue), 2])
multi_reg = intersect(multi_sig, twas_multi$ID)

##Import Multi-omics weight#
multiFile <- dir(paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/output/", tissue, sep = ""),
                pattern = paste("_chr", chr_num, "_", sep = ""))
multi.weight = read.table(paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/output/", 
				tissue, "/", multiFile, sep = ""), header = TRUE)
multi.weight = multi.weight[multi.weight$gene_id %in% multi_sig, ]
multi.weight$chr = sapply(strsplit(multi.weight$snp, "_"), function(x){as.character(x[1])})
multi.weight$pos = sapply(strsplit(multi.weight$snp, "_"), function(x){as.character(x[2])})
multi.weight$ref = sapply(strsplit(multi.weight$snp, "_"), function(x){as.character(x[3])})
multi.weight$alt = sapply(strsplit(multi.weight$snp, "_"), function(x){as.character(x[4])})
multi.weight$w = -1*multi.weight$weight
multi.weight$gene = multi.weight$gene_id
multi.weight$chrpos = paste(multi.weight$chr, multi.weight$pos, sep = ":")
multi.weight = multi.weight %>% select(gene, chr, pos, ref, alt, w, chrpos, snp)

espresso.weight = multi.weight

#################################################################################################################################
##Load in internal validation for predixcan##
method = "predixcan"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/internal_validation/", tissue, sep = "")
filename = dir(path, pattern = paste("_chr", chr_num, "_", sep = ""))
int = data.frame()
for (i in 1:length(filename)) {
        temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
        int = rbind(int, temp)
}
predixcan_sig = (int %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05))$gene_id

##Load in TWAS result##
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/twas_assoc_new/", method, "/", pheno, "/", tissue, "/output/chr", chr_num, "_result.txt", sep = "")
twas_predixcan = read.table(path, header = TRUE) %>% filter(TWAS.P < 0.05/count[which(count$tissue == tissue), 2])
predixcan_reg = intersect(predixcan_sig, twas_predixcan$ID)

##Import PrediXcan weight#
predixcanFile <- dir(paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/output/", tissue, sep = ""),
                pattern = paste("_chr", chr_num, "_", sep = ""))
predixcan.weight = read.table(paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/output/",
                                tissue, "/", predixcanFile, sep = ""), header = TRUE)
predixcan.weight = predixcan.weight[predixcan.weight$gene_id %in% predixcan_sig, ]
predixcan.weight$chr = sapply(strsplit(predixcan.weight$snp, "_"), function(x){as.character(x[1])})
predixcan.weight$pos = sapply(strsplit(predixcan.weight$snp, "_"), function(x){as.character(x[2])})
predixcan.weight$ref = sapply(strsplit(predixcan.weight$snp, "_"), function(x){as.character(x[3])})
predixcan.weight$alt = sapply(strsplit(predixcan.weight$snp, "_"), function(x){as.character(x[4])})
predixcan.weight$w = -1*predixcan.weight$weight
predixcan.weight$gene = predixcan.weight$gene_id
predixcan.weight$chrpos = paste(predixcan.weight$chr, predixcan.weight$pos, sep = ":")
predixcan.weight = predixcan.weight %>% select(gene, chr, pos, ref, alt, w, chrpos, snp)

#################################################################################################################################
##Load in internal validation for dapg##
method = "dapg"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/internal_validation/", tissue, sep = "")
filename = dir(path, pattern = paste("_chr", chr_num, "_", sep = ""))
int = read.table(paste(path, "/", filename, sep =""), header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
dapg_sig = (int %>% filter(pear_avg_t > 0.1) %>% filter(pear_stouffer_pval < 0.05))$gene_id

##Only for dapg, we also need to filter out genes with no variants with PIP > 0.01##
load(paste("/gpfs/group/dxl46/default/private/poom/dapg/prior/", tissue, "/", "prior_chr", chr_num, ".RData", sep = ""))
final_df = final_df %>% filter(pip > 0.01)
dapg_pip = unique(final_df$gene)
dapg_sig = intersect(dapg_sig, dapg_pip)

##Load in TWAS result##
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/twas_assoc_new/", method, "/", pheno, "/", tissue, "/output/chr", chr_num, "_result.txt", sep = "")
twas_dapg = read.table(path, header = TRUE) %>% filter(TWAS.P < 0.05/count[which(count$tissue == tissue), 2])
dapg_reg = intersect(dapg_sig, twas_dapg$ID)

##Import DAPG weight#
dapgFile <- dir(paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/output/", tissue, sep = ""),
                pattern = paste("_chr", chr_num, "_", sep = ""))
dapg.weight = read.table(paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/output/",
                                tissue, "/", dapgFile, sep = ""), header = TRUE)
dapg.weight = dapg.weight[dapg.weight$gene_id %in% dapg_sig, ]
dapg.weight$chr = sapply(strsplit(dapg.weight$snp, "_"), function(x){as.character(x[1])})
dapg.weight$pos = sapply(strsplit(dapg.weight$snp, "_"), function(x){as.character(x[2])})
dapg.weight$ref = sapply(strsplit(dapg.weight$snp, "_"), function(x){as.character(x[3])})
dapg.weight$alt = sapply(strsplit(dapg.weight$snp, "_"), function(x){as.character(x[4])})
dapg.weight$w = -1*dapg.weight$weight
dapg.weight$gene = dapg.weight$gene_id
dapg.weight$chrpos = paste(dapg.weight$chr, dapg.weight$pos, sep = ":")
dapg.weight = dapg.weight %>% select(gene, chr, pos, ref, alt, w, chrpos, snp)

#################################################################################################################################
##Load in significant model for mashr##
method = "mashr"
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/sig/", tissue, "/", tissue, "_allchr_GTEx_mashr_sig.txt", sep = "")
mashr_sig = read.table(path, header = TRUE)$sig_gene

##Load in TWAS result##
path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/twas_assoc_new/", method, "/", pheno, "/", tissue, "/output/chr", chr_num, "_result.txt", sep = "")
twas_mashr = read.table(path, header = TRUE) %>% filter(TWAS.P < 0.05/count[which(count$tissue == tissue), 2])
mashr_reg = intersect(mashr_sig, twas_mashr$ID)

##Import mashr weight#
mashrFile <- dir(paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/output/", tissue, sep = ""))
mashr.weight = read.table(paste("/gpfs/group/dxl46/default/private/poom/GTEx/coef_total/", method, "/output/",
                                tissue, "/", mashrFile, sep = ""), header = TRUE)
mashr.weight$chr = sapply(strsplit(mashr.weight$snp, "_"), function(x){as.character(x[1])})
mashr.weight = mashr.weight %>% filter(chr == chr_num)
mashr.weight = mashr.weight[mashr.weight$gene_id %in% mashr_sig, ]
mashr.weight$pos = sapply(strsplit(mashr.weight$snp, "_"), function(x){as.character(x[2])})
mashr.weight$ref = sapply(strsplit(mashr.weight$snp, "_"), function(x){as.character(x[3])})
mashr.weight$alt = sapply(strsplit(mashr.weight$snp, "_"), function(x){as.character(x[4])})
mashr.weight$w = 1*mashr.weight$weight
mashr.weight$gene = mashr.weight$gene_id
mashr.weight$chrpos = paste(mashr.weight$chr, mashr.weight$pos, sep = ":")
mashr.weight = mashr.weight %>% select(gene, chr, pos, ref, alt, w, chrpos, snp)

#################################################################################################################################
##Identify genes common to all 5 methods##
genes_common = intersect(unique(mashr.weight$gene), intersect(unique(dapg.weight$gene), intersect(unique(predixcan.weight$gene), intersect(unique(espresso.weight$gene), unique(utmost.weight$gene)))))
predixcan.weight = predixcan.weight[predixcan.weight$gene %in% genes_common, ]
espresso.weight = espresso.weight[espresso.weight$gene %in% genes_common, ]
utmost.weight = utmost.weight[utmost.weight$gene %in% genes_common, ]
dapg.weight = dapg.weight[dapg.weight$gene %in% genes_common, ]
mashr.weight = mashr.weight[mashr.weight$gene %in% genes_common, ]

##Identify regions to be included##
genes_reg = union(unique(mashr_reg), union(unique(dapg_reg) ,union(unique(multi_reg), union(unique(predixcan_reg), unique(utmost_reg)))))

#################################################################################################################################
##Determine local genes within each LD block##
##ld = read.table("/storage/home/cxk502/focus/pyfocus/data/ld_blocks/grch37.eur.loci.bed", header = TRUE) %>% filter(chrom == chr_num)
ld = read.table("/gpfs/group/dxl46/default/private/poom/espresso/grch37.eur.loci.bed", header = TRUE) %>% filter(chrom == chr_num)
ld$range = paste(ld$chrom, ":", ld$start, "-", ld$stop, sep = "")
genes = read.table("/gpfs/group/dxl46/default/private/poom/FOCUS/ensembl_ref.txt", header = TRUE) %>% filter(chrom == chr_num)

##Map genes into the ld block##
ld_range = makeGRangesFromDataFrame(ld, seqnames.field=c("chrom"),start.field="start", end.field=c("stop"), ignore.strand=TRUE)
gene_range = makeGRangesFromDataFrame(genes, seqnames.field=c("chrom"),start.field="tx_start", end.field=c("tx_stop"), ignore.strand=TRUE)

fo = findOverlaps(query=gene_range, subject=ld_range, type = "within")
df = as.data.frame(cbind(ld[subjectHits(fo), "range"], genes[queryHits(fo), "ens_gene_id"]))
df_reg = df[df$V2 %in% genes_reg, ]

##Identify region with significant GWAS hit##
#gwas_sig = gwas.result %>% filter(p < 5*10^(-8))
gwas_sig = gwas.result %>% filter(p <= 1)

if (nrow(gwas_sig) >= 1 & nrow(df_reg) >= 1){
	snp_range = makeGRangesFromDataFrame(gwas_sig, seqnames.field=c("chr"),start.field="pos", end.field=c("pos"), ignore.strand=TRUE)
	fo = findOverlaps(query=snp_range, subject=ld_range, type = "within")
	sig_loci = as.data.frame(cbind(ld[subjectHits(fo), "range"], gwas_sig[queryHits(fo), "id"]))

	##Update ld block##
	df = df[df$V1 %in% unique(sig_loci$V1), ]
	ld = ld[ld$range %in% df$V1, ]

	gwas.result = gwas.result %>% select(-c("p"))
#################################################################################################################################

	##Perform fine-mapping##
	cauchymu_result = data.frame()
	cauchypu_result = data.frame()
	cauchydu_result = data.frame()
	cauchymau_result = data.frame()
	cauchymam_result = data.frame()
	cauchymap_result = data.frame()
	cauchymad_result = data.frame()

	multi_result = data.frame()
	predixcan_result = data.frame()
	utmost_result = data.frame()
	dapg_result = data.frame()
	mashr_result = data.frame()
	for (idx in 1:length(unique(df_reg$V1))){
    		gene.locale = (df %>% filter(V1 == unique(df_reg$V1)[idx]))$V2

    		if (length(intersect(unique(intersect(utmost.weight$gene,espresso.weight$gene)),gene.locale)) > 0){
    			x = espressso.twas.cor(utmost.weight,espresso.weight,gwas.result,vcf.ref.file,gene.locale)
			temp = espresso.output(x, tau_value)
			cauchymu_result = rbind(cauchymu_result, temp[[1]])
			multi_result = rbind(multi_result, temp[[2]])
			utmost_result = rbind(utmost_result, temp[[3]])
			rm(x, temp)

			x = espressso.twas.cor(utmost.weight,predixcan.weight,gwas.result,vcf.ref.file,gene.locale)
			temp = espresso.output(x, tau_value)
                        cauchypu_result = rbind(cauchypu_result, temp[[1]])
                        predixcan_result = rbind(predixcan_result, temp[[2]])
			rm(x, temp)

			x = espressso.twas.cor(utmost.weight,dapg.weight,gwas.result,vcf.ref.file,gene.locale)
			temp = espresso.output(x, tau_value)
                        cauchydu_result = rbind(cauchydu_result, temp[[1]])
                        dapg_result = rbind(dapg_result, temp[[2]])
                        rm(x, temp)

			x = espressso.twas.cor(utmost.weight,mashr.weight,gwas.result,vcf.ref.file,gene.locale)
			temp = espresso.output(x, tau_value)
                        cauchymau_result = rbind(cauchymau_result, temp[[1]])
                        mashr_result = rbind(mashr_result, temp[[2]])
                        rm(x, temp)

			x = espressso.twas.cor(mashr.weight,espresso.weight,gwas.result,vcf.ref.file,gene.locale)
                        temp = espresso.output(x, tau_value)
                        cauchymam_result = rbind(cauchymam_result, temp[[1]])
                        rm(x, temp)

			x = espressso.twas.cor(mashr.weight,predixcan.weight,gwas.result,vcf.ref.file,gene.locale)
                        temp = espresso.output(x, tau_value)
                        cauchymap_result = rbind(cauchymap_result, temp[[1]])
                        rm(x, temp)		
		
			x = espressso.twas.cor(mashr.weight,dapg.weight,gwas.result,vcf.ref.file,gene.locale)
                        temp = espresso.output(x, tau_value)
                        cauchymad_result = rbind(cauchymad_result, temp[[1]])
                        rm(x, temp)	

			write.table(cauchymu_result, paste("/gpfs/group/dxl46/default/private/poom/espresso/espresso_finemap_new1/", pheno, "/", tissue, "/cauchymu_chr", chr_num, "_", tau_value, "_result.txt", sep = ""),  
			col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
		
			write.table(cauchypu_result, paste("/gpfs/group/dxl46/default/private/poom/espresso/espresso_finemap_new1/", pheno, "/", tissue, "/cauchypu_chr", chr_num, "_", tau_value, "_result.txt", sep = ""),
                        col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
			
			write.table(cauchydu_result, paste("/gpfs/group/dxl46/default/private/poom/espresso/espresso_finemap_new1/", pheno, "/", tissue, "/cauchydu_chr", chr_num, "_", tau_value, "_result.txt", sep = ""),
                        col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)		

			write.table(cauchymau_result, paste("/gpfs/group/dxl46/default/private/poom/espresso/espresso_finemap_new1/", pheno, "/", tissue, "/cauchymau_chr", chr_num, "_", tau_value, "_result.txt", sep = ""),
                        col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

			write.table(cauchymam_result, paste("/gpfs/group/dxl46/default/private/poom/espresso/espresso_finemap_new1/", pheno, "/", tissue, "/cauchymam_chr", chr_num, "_", tau_value, "_result.txt", sep = ""),
                        col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

			write.table(cauchymap_result, paste("/gpfs/group/dxl46/default/private/poom/espresso/espresso_finemap_new1/", pheno, "/", tissue, "/cauchymap_chr", chr_num, "_", tau_value, "_result.txt", sep = ""),
                        col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

			write.table(cauchymad_result, paste("/gpfs/group/dxl46/default/private/poom/espresso/espresso_finemap_new1/", pheno, "/", tissue, "/cauchymad_chr", chr_num, "_", tau_value, "_result.txt", sep = ""),
                        col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

			write.table(multi_result, paste("/gpfs/group/dxl46/default/private/poom/espresso/espresso_finemap_new1/", pheno, "/", tissue, "/multi-omics_chr", chr_num, "_", tau_value, "_result.txt", sep = ""),
        		col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

			write.table(utmost_result, paste("/gpfs/group/dxl46/default/private/poom/espresso/espresso_finemap_new1/", pheno, "/", tissue, "/utmost_chr", chr_num, "_", tau_value, "_result.txt", sep = ""),
        		col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

			write.table(predixcan_result, paste("/gpfs/group/dxl46/default/private/poom/espresso/espresso_finemap_new1/", pheno, "/", tissue, "/predixcan_chr", chr_num, "_", tau_value, "_result.txt", sep = ""),
        		col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

			write.table(dapg_result, paste("/gpfs/group/dxl46/default/private/poom/espresso/espresso_finemap_new1/", pheno, "/", tissue, "/dapg_chr", chr_num, "_", tau_value, "_result.txt", sep = ""),
                        col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

			write.table(mashr_result, paste("/gpfs/group/dxl46/default/private/poom/espresso/espresso_finemap_new1/", pheno, "/", tissue, "/mashr_chr", chr_num, "_", tau_value, "_result.txt", sep = ""),
                        col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

		} else {
			next
		}
	}
}

dummy="I am done"
filename = paste("/gpfs/group/dxl46/default/private/poom/espresso/dummy_finemap_new1/", pheno, "/", tissue, "/dummy_chr", chr_num, "_", tau_value, ".txt", sep = "")
write.table(dummy, filename, quote=FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
