args = commandArgs(trailingOnly=TRUE)

if_verbose = FALSE

###################################################################################################################################################################################
##Copy from UTMOST code##
### optimization part ###
## pre-calculate some metrics for gradient
## args
## X: a list of covariate matrices corresponding to each response
## Y: a list of response vectors
## value
## XY: a list of matrices X^TY for each response
grad_prep <- function(X, Y){
        ll = length(Y)
        P = ncol(X[[1]])
        XY = matrix(0,P,ll)
        for(i in 1:ll){
                XY[,i] = t(X[[i]])%*%Y[[i]]/nrow(X[[i]])
        }
	XY
}

## get the minimum and maximum of lambda searched in cross-validation of an elastic net model
## args
## lst: an object returned by glmnet
## value
## min_lam: smallest lambda searched in glmnet cross-validation
## max_lam: largest lambda searched in glmnet cross-validation
minmax_lambda <- function(lst){
        max_lam = max(unlist(lapply(lst, function(x){max(x$lambda)})))
        min_lam = min(unlist(lapply(lst, function(x){min(x$lambda)})))
        c(min_lam, max_lam)
}

## evaluate the performance of elastic net on each response
## args
## lst: a list of glmnet object (fitted elastic net model for each response)
## X_tune: a list of covariate matrices corresponding for each response (for tuning lambda)
## Y_tune: a list of response vectors (for tuning lambda)
## X_test: a list of covariate matrices corresponding for each response (for testing performance)
## Y_test: a list of response vectors (for testing performance)
## value
## lam: best performing lambda (on (X_tune,Y_tune)) for each response
## mse: list of matrices with each element being a matrix of predicted vs observed response
## est: estimated effect sizes for each response (B matrix)
elastic_net_mse <- function(lst, X_tune, Y_tune, X_test, Y_test){
        P = length(lst)
        M = ncol(X_tune[[1]])
        lam_V = rep(0, P)
        test_res = list()
        test_beta = matrix(0, M, P)
        for(t in 1:P){
                ncv = length(lst[[t]]$lambda)
                tmp_mse = rep(0, ncv)
                for(k in 1:ncv){
                        tmp_mse[k] = mean((Y_tune[[t]] - X_tune[[t]]%*%lst[[t]]$glmnet.fit$beta[,k])^2)
                }
                ss = which.min(tmp_mse)
                test_beta[,t] = lst[[t]]$glmnet.fit$beta[,ss]
                lam_V[t] = lst[[t]]$lambda[ss]
                predicted = X_test[[t]]%*%lst[[t]]$glmnet.fit$beta[,ss]
                test_res[[t]] = cbind(Y_test[[t]], predicted)
        }
	list(lam = lam_V, mse = test_res, est = test_beta)
}

## create a list of Y_true and Y_predicted for analysis
multi_mse <- function(theta_est, X_test, Y_test){
  answer = list()
        P = ncol(theta_est)
        for(t in 1:P){
                predicted = X_test[[t]]%*%theta_est[,t]
                answer[[t]] = cbind(Y_test[[t]], predicted)
        }
	answer
}

## average prediction metrics (rsq, mse, adjusted mse on validation sets)
avg_perm <- function(mse_lst){
        fd = length(mse_lst)
        P = length(mse_lst[[1]])
        rsq = mse = adj_mse = matrix(0, fd, P)
        for(f in 1:fd){
                for(t in 1:P){
                        rsq[f,t] = (cor(mse_lst[[f]][[t]])[1,2])^2
                        mse[f,t] = mean((mse_lst[[f]][[t]][,1]-mse_lst[[f]][[t]][,2])^2)
                        adj_mse[f,t] = mse[f,t]/var(mse_lst[[f]][[t]][,1])
                }
        }
	cbind(apply(rsq, 2, mean), apply(mse, 2, mean), apply(adj_mse, 2, mean))
}

## test Y_true ~ Y_predicted
pred_test <- function(Y){
        if(sum(Y[,2]==0)==nrow(Y)|var(Y[,2])==0){
                return(2)
        }else{
              	summary(lm(Y[,1]~Y[,2]))$coefficients[2,4]
        }
}

## group lasso on data with missing covariates, use validation data for stopping criterion
glasso <- function(X, Y, X1, Y1, XX, XY, Xnorm, lambda1, lambda2, theta, stepsize = 1e-4, maxiter = 50, eps = 1e-3, verbose = FALSE){
        Bgt = Sys.time()
        M = nrow(XY)
        P = length(X)
        NN = unlist(lapply(X, nrow))
        old_objV1 = rep(0,P)
        for(t in 1:P){
                old_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2)
        }
  if (verbose) {
    cat("Training error: ", old_objV1, '\n')
  }
   	old_objV2 = rep(0,P)
        for(t in 1:P){
                old_objV2[t] = 1/2*mean((Y1[[t]]-X1[[t]]%*%theta[,t])^2)
        }
	if (verbose) {
    cat("Testing error: ", old_objV2, '\n')
  }
   	beta_j_lasso = rep(0, P)
        tmp_XYj = 0
        if(!is.loaded("wrapper")){
                dyn.load("/gpfs/group/dxl46/default/private/poom/GTEx/UTMOST/CTIMP/optim.so") # change this to the abs path to optim.so
        }
	for(i in 1:maxiter){
                bgt = Sys.time()
                res = .Call("wrapper", XX, XY, theta, M, P, beta_j_lasso, lambda1, lambda2, Xnorm)
                edt = Sys.time()

                new_objV1 = new_objV2 = rep(0,P)
                for(t in 1:P){
                        new_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2)
                }
                if (verbose) {cat("Training error: ", new_objV1, '\n')}
                for(t in 1:P){
                        new_objV2[t] = 1/2*mean((Y1[[t]]-X1[[t]]%*%theta[,t])^2)
                }
                if (verbose) {cat("Testing error: ", new_objV2, '\n')}
                if(mean(new_objV2) > mean(old_objV2)|mean(new_objV1) > mean(old_objV1)){
                        break
                }else{
                      	old_objV2 = new_objV2
                }
                if(max(abs(new_objV1-old_objV1)) < eps){
                        break
                }else{
                      	old_objV1 = new_objV1
                }
        }
	Edt = Sys.time()
        cat("total training time: ", Edt-Bgt, "\n")
        list(est = theta, avg_tune_err = mean(new_objV2), tune_err=new_objV2)
}

## simpler version of glasso, train model until converges
glasso_no_early_stopping <- function(X, Y, XX, XY, Xnorm, lambda1, lambda2, theta, stepsize = 1e-4, maxiter = 50, eps = 1e-3, verbose = FALSE){
  cat("running glasso_no_early_stopping\n")
        M = nrow(XY)
        P = length(X)
        NN = unlist(lapply(X, nrow))
        old_objV1 = rep(0,P)
        for(t in 1:P){
                old_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2)
        }
	if (verbose) {cat("Training error: ", mean(old_objV1), '\n')}
        beta_j_lasso = rep(0, P)
        tmp_XYj = 0
        if(!is.loaded("wrapper")){
                dyn.load("optim.so")
        }
	for(i in 1:maxiter){
                res = .Call("wrapper", XX, XY, theta, M, P, beta_j_lasso, lambda1, lambda2, Xnorm)
                new_objV1 = rep(0,P)
                for(t in 1:P){
                        new_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2)
                }
                if (verbose) {cat("Training error: ", mean(new_objV1), '\n')}
                if(max(abs(new_objV1-old_objV1)) < eps|mean(new_objV1) > mean(old_objV1)){
                        break
                }else{
                      	old_objV1 = new_objV1
                }
        }
	list(est = theta, avg_train_err = mean(new_objV1), train_err = new_objV1)
}

###################################################################################################################################################################################
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


## Function for filling in missing dosage
## args
## geno: genotype dataframe
## values
## geno_new: genotype dataframe with filled in dosage
geno_fill <- function(geno, type){

        ##Filling in missing data##
        geno = as.matrix(geno)
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
library(foreach)

#########################################################################################################################
                                        ##Processing genotype##
#########################################################################################################################
chr_num = args[1]
total_file_num = as.numeric(args[2])
file_num = as.numeric(args[3])
fold = as.numeric(args[4])
dir_path = normalizePath(args[5])

geno_path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/input/genotype/GTEx_chr", chr_num, "_processed.traw", sep="")

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

##Filling in missing data according to the training and validating set##
geno = geno_fill(geno)

##Standardize genotype matrix (only center)##
center_x_train = rowMeans(geno, na.rm = T)
geno = as.data.frame(t(scale(t(geno), center = center_x_train, scale = FALSE)))

#########################################################################################################################
                                        ##Processing expression##
#########################################################################################################################
tissue_list = c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary",
                "Artery_Tibial", "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia",
                "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9",
                "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia",
                "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes",
                "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa",
                "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal",
                "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg",
                "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")

expression = list()
gene_list = data.frame()
for ( i in 1:length(tissue_list) ) {
	print(i)
	##Process expression input##
	tissue = tissue_list[i]
	expression_path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/input/expression_processed/", tissue, "_normalized_expression_final.txt",  sep="")
	expression[[i]] = read.table(expression_path, header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE) %>% filter(chromosome == chr_num)
	names(expression[[i]])[5:ncol(expression[[i]])] = sample_conversion$sample_ID[match(names(expression[[i]])[5:ncol(expression[[i]])], sample_conversion$subject_ID)]

	##Create gene list for mapping##
	gene_list_temp  = expression[[i]] %>% select(gene_id, chromosome, start, end) %>% remove_rownames %>% column_to_rownames(var="gene_id")
	gene_list = rbind(gene_list, gene_list_temp)
	
	##Center expression##
	expression[[i]] = expression[[i]] %>% select(-c("chromosome", "start", "end")) %>% remove_rownames %>% column_to_rownames(var="gene_id")
	center_y_train = rowMeans(expression[[i]], na.rm = T)
	expression[[i]] = as.data.frame(t(scale(t(expression[[i]]), center = center_y_train, scale = FALSE)))
}
gene_list = unique(gene_list)
set.seed(123)
rows = sample(nrow(gene_list))
gene_list = gene_list[rows, ]
gene_list_idx = file_helper(nrow(gene_list), total_file_num)
gene_list = gene_list[gene_list_idx[file_num,1]:gene_list_idx[file_num,2],]

##For long job##
#gene_list_idx = file_helper(nrow(gene_list), 10)
#gene_list = gene_list[gene_list_idx[file_num_new,1]:gene_list_idx[file_num_new,2],]


#########################################################################################################################
                                        ##Processing overall##
#########################################################################################################################
##Filling in SNP coordinates##
geno = geno_coord(geno)

func_type = 1000
map = constant_process(geno, gene_list, chr_num, func_type)	

##Update the gene list in expression dataframe according to map file##
geno = geno %>% select(-c("chromosome", "start", "end"))

##Update expression list according to map##
for (i in seq(1, length(tissue_list))) {
	expression[[i]] = expression[[i]][rownames(expression[[i]]) %in% map$gene_id,]
}

#########################################################################################################################
                                        ##Create CV folds##
#########################################################################################################################
##Create cv fold
set.seed(123)
cv_fold_total = createFolds(sample_conversion$sample_ID, k = fold, list = TRUE, returnTrain = FALSE)
cv_fold = data.frame()
for (f in 1:fold) {
        temp = cbind(sample = sample_conversion$sample_ID[cv_fold_total[[f]]], foldid = f)
        cv_fold = rbind(cv_fold, temp)
}
#########################################################################################################################
					###Perform UTMOST##
#########################################################################################################################
internal_validate = vector("list", length(tissue_list))
coef = vector("list", length(tissue_list))
summary = vector("list", length(tissue_list))

for (i in 1:nrow(map)) {
	origin_id = as.character(map[i,1])
	
	##Identify which tissue has a specific gene expression##
        tissue_gene = c()
        for ( index in 1:length(tissue_list) ) {
                if ( nrow(expression[[index]][rownames(expression[[index]]) %in% map$gene_id[i],]) == 0) {
                        tissue_gene = tissue_gene
                } else {
                        tissue_gene = c(tissue_gene, index)
                }
        }

	##Create expression matrix##
	Y = list()
	for ( t in 1:length(tissue_gene) ) {

		Y[[t]] = t(expression[[tissue_gene[t]]][rownames(expression[[tissue_gene[t]]]) %in% map$gene_id[i],])
	}
	ssize = unlist(lapply(Y, nrow))
        T_num = length(Y)
	
	dose = t(geno[unlist(strsplit(as.character(map[i,3]), ";")),])
	XX = t(dose)%*%as.matrix(dose)/(nrow(dose))
	Xnorm = diag(XX)
	remove(XX)

	sub_id = as.numeric(as.character(rownames(dose)))
	sub_id_map = list()

	for (t in 1:length(tissue_gene)){
		tmp = rep(0, nrow(Y[[t]]))
		for(j in 1:length(tmp)){
                        tmp[j] = which(sub_id == rownames(Y[[t]])[j])
                }
                sub_id_map[[t]] = tmp
        }
		
	single_res_test = list()
        single_lam = matrix(0,fold,length(tissue_gene))
        single_theta_est = list()

        multi_res_tune = list()

        multi_res_test = list()
        multi_lam = matrix(0,fold,2)
        multi_theta_est = list()

        multi_res_test2 = list()
        multi_lam2 = array(0, dim=c(fold, length(tissue_gene), 2))
        multi_theta_est2 = list()

	#-------------------------- dividing cross-validation sets ---------------------------#
	res_tune = list()
	rec_lamv = matrix(0, fold, fold)

	for(f in 1:fold){
		bgt = Sys.time()
		test_index = as.numeric(as.character(as.list(cv_fold %>% filter(foldid == f) %>% select(sample))$sample))
		test_id = sub_id[test_index]
		tuning_index = as.numeric(as.character(as.list(cv_fold %>% filter(foldid == (f%%fold)+1) %>% select(sample))$sample))
		tuning_id = sub_id[tuning_index]

		X_test = list()
                Y_test = list()
                X_tune = list()
                Y_tune = list()
                X_train = list()
                Y_train = list()
		for(t in 1:length(tissue_gene)){
			X_train_tmp = sub_id_map[[t]][!(sub_id_map[[t]]%in%c(test_index,tuning_index))]
			Y_train_tmp = !(sub_id_map[[t]]%in%c(test_index,tuning_index))
			X_tuning_tmp = sub_id_map[[t]][(sub_id_map[[t]]%in%tuning_index)]
                        Y_tuning_tmp = (sub_id_map[[t]]%in%tuning_index)
                        X_test_tmp = sub_id_map[[t]][(sub_id_map[[t]]%in%test_index)]
                        Y_test_tmp = (sub_id_map[[t]]%in%test_index)
			X_train[[t]] = apply(as.matrix(dose[X_train_tmp,]),2,as.numeric)
			Y_train[[t]] = Y[[t]][Y_train_tmp, 1]
			X_tune[[t]] = apply(as.matrix(dose[X_tuning_tmp,]),2,as.numeric)
                        Y_tune[[t]] = Y[[t]][Y_tuning_tmp, 1]
                        X_test[[t]] = apply(as.matrix(dose[X_test_tmp,]),2,as.numeric)
                        Y_test[[t]] = Y[[t]][Y_test_tmp, 1]
		}
		
		#------------- get initial est by elasticNet on single tissue --------------#
		single_initial_est = matrix(0, ncol(X_train[[1]]), length(tissue_gene))
		single_summary = list()
		
		##Do it to get a list of lambda##
		for(t in 1:length(tissue_gene)){
                        tt = cv.glmnet(X_train[[t]], Y_train[[t]], alpha = 0.5, nfolds = 5)
                        single_summary[[t]] = tt
                        single_initial_est[,t] = tt$glmnet.fit$beta[,which.min(tt$cvm)]
                }
			
		##Performance of Elastic net on tuning and testing data with various tuning parameters##		
		els_output = elastic_net_mse(single_summary, X_tune, Y_tune, X_test, Y_test)
                single_res_test[[f]] = els_output$mse
                single_lam[f,] = els_output$lam
                single_theta_est[[f]] = els_output$est
                remove(els_output)

		##Use elastic net ests row norm as weights##
                lam_range = minmax_lambda(single_summary)
                sig_norm = apply(single_initial_est, 1, function(x){sqrt(sum(x^2))})
                sig_norm[sig_norm==0] = rep(min(sig_norm[sig_norm>0]), sum(sig_norm==0))/2
                sig_norm = sig_norm/sum(sig_norm)
                weights2 = 1/sig_norm; weights2 = weights2/sum(weights2);

		tis_norm = apply(single_initial_est, 2, function(x){sum(abs(x))})
                tis_norm[tis_norm==0] = rep(min(tis_norm[tis_norm>0]), sum(tis_norm==0))/2
                tis_norm = tis_norm/sum(tis_norm)
                weights1 = 1/tis_norm; weights1 = weights1/sum(weights1);
                lam_V = seq(lam_range[1], lam_range[2], length.out = fold)
	
		initial_numeric = as.numeric(single_initial_est) ##List of weights of each tissue in continuous
                remove(single_summary); remove(single_initial_est);

		#-------------------------- train - validate - test ---------------------------#
		XY = grad_prep(X_train, Y_train)
                XX_train = lapply(X_train, function(x){t(x)%*%x/nrow(x)})
                spsz = unlist(lapply(X_train,nrow))
                res_tune[[f]] = array(-1, dim=c(fold, fold, length(tissue_gene)))
                rec_lamv[f,] = lam_V

		for(lam1 in 1:fold){
                        for(lam2 in 1:fold){
                                single_est = matrix(initial_numeric, ncol(dose), length(tissue_gene))
                                ans = glasso(X=X_train, Y=Y_train, X1=X_tune, Y1=Y_tune, XX=XX_train, XY=XY, Xnorm=Xnorm, lambda1=lam_V[lam1]/spsz, lambda2=lam_V[lam2], theta=single_est, verbose = if_verbose)
                                if(sum(ans$est!=0)>0){
                                        res_tune[[f]][lam1,lam2, ] = ans$tune_err
                                        if (if_verbose) { cat("lambda1=",lam_V[lam1], "; lambda2=", lam_V[lam2], "; avg tune err=", ans$avg_tune_err, '\n') }
                                        remove(single_est); remove(ans);
                                }else{
                                      	remove(single_est); remove(ans);
                                        break
                                }
                        }
                }
		
		#-------------------------- save results on test set for evaluation ---------------------------#
		avg_tune_res = apply(res_tune[[f]], c(1,2), mean)
                best.lam = which(avg_tune_res == min(avg_tune_res[avg_tune_res>=0]), arr.ind = TRUE)[1,]
                single_est = matrix(initial_numeric, ncol(dose), length(tissue_gene))
                ans = glasso(X=X_train, Y=Y_train, X1=X_tune, Y1=Y_tune, XX=XX_train, XY=XY, Xnorm=Xnorm, lambda1=lam_V[best.lam[1]]/spsz, lambda2=lam_V[best.lam[2]], theta=single_est, verbose = if_verbose)
                multi_res_tune[[f]] = multi_mse(ans$est, X_tune, Y_tune)
                multi_res_test[[f]] = multi_mse(ans$est, X_test, Y_test)
                multi_lam[f,] = lam_V[best.lam]
                multi_theta_est[[f]] = ans$est
                remove(single_est); remove(ans);
	}

###############################################################################################################################################################################################
	##Determine performance in internal validation##
        pear_folds = rep(0,fold)
        spear_folds = rep(0,fold)
        pear_zscore_folds = rep(0,fold)
        spear_zscore_folds = rep(0,fold)
	
	for (t in 1:length(tissue_gene)) {
                for (f in 1:fold) {
                        pear_folds[f] =  ifelse(sd(multi_res_test[[f]][[t]][,2]) != 0, cor( multi_res_test[[f]][[t]][,2], multi_res_test[[f]][[t]][,1], method = "pearson"), 0)
                        spear_folds[f] =  ifelse(sd(multi_res_test[[f]][[t]][,2]) != 0, cor( multi_res_test[[f]][[t]][,2], multi_res_test[[f]][[t]][,1], method = "spearman"), 0)

                        pear_zscore_folds[f] <- atanh(pear_folds[f])*sqrt(nrow(multi_res_test[[f]][[t]]) - 3) #Fisher transformation
                        spear_zscore_folds[f] <- atanh(spear_folds[f])*sqrt(nrow(multi_res_test[[f]][[t]]) - 3) #Fisher transformation

                }

                ##Calculate average of correlation across folds##
                pear_avg = mean(pear_folds)
                spear_avg = mean(spear_folds)

                ##Combine Z-scores via Stouffer's method##
                pear_zscore_est = sum(pear_zscore_folds) / sqrt(fold)
                pear_stouffer_pval <- 2*pnorm(abs(pear_zscore_est), lower.tail = FALSE)
                spear_zscore_est = sum(spear_zscore_folds) / sqrt(fold)
                spear_stouffer_pval <- 2*pnorm(abs(spear_zscore_est), lower.tail = FALSE)

		els_output = data.frame("tissue" = tissue_list[tissue_gene[t]], "gene_id" = map$gene_id[i], "type" = func_type, "region" = origin_id,
                                "pear_avg_t" = pear_avg, "pear_stouffer_pval" = pear_stouffer_pval, 
				"spear_avg_t" = spear_avg, "spear_stouffer_pval" = spear_stouffer_pval) %>% remove_rownames
                internal_validate[[tissue_gene[t]]] = rbind(internal_validate[[tissue_gene[t]]], els_output)
        
		filename_result <- paste(dir_path, "/utmost/internal_validation/", tissue_list[tissue_gene[t]], "/", tissue_list[tissue_gene[t]], 
					 "_chr", chr_num, "_GTEx_utmost_internal_validation_", total_file_num, ".", file_num, ".txt", sep="")
                write.table(internal_validate[[tissue_gene[t]]], filename_result, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
	}
	
###############################################################################################################################################################################################

		#------------ use tuning parameter chosen above to train model on entire dataset -------------#
        ## generate an estimate with whole data ##
	cat('training a model on entire data with parameters chosen from cv\n')
        X_all = list()
        Y_all = list()
        for(t in 1:T_num){
                X_all_tmp = sub_id_map[[t]]
                X_all[[t]] = apply(as.matrix(dose[X_all_tmp,]),2,as.numeric)
                Y_all[[t]] = Y[[t]][,1]
        }
	# initial values
        single_initial_est = matrix(0, ncol(X_train[[1]]), T_num)
        for(t in 1:T_num){
                tt = cv.glmnet(X_all[[t]], Y_all[[t]], alpha = 0.5, nfolds = 5)
                single_initial_est[,t] = tt$glmnet.fit$beta[,which.min(tt$cvm)]
        }

	sig_norm = apply(single_initial_est, 1, function(x){sqrt(sum(x^2))})
        sig_norm[sig_norm==0] = rep(min(sig_norm[sig_norm>0]), sum(sig_norm==0))/2
        sig_norm = sig_norm/sum(sig_norm)
        weights2 = 1/sig_norm; weights2 = weights2/sum(weights2);

        tis_norm = apply(single_initial_est, 2, function(x){sum(abs(x))})
        tis_norm[tis_norm==0] = rep(min(tis_norm[tis_norm>0]), sum(tis_norm==0))/2
        tis_norm = tis_norm/sum(tis_norm)
        weights1 = 1/tis_norm; weights1 = weights1/sum(weights1);

	spsz = unlist(lapply(X_all,nrow))
        initial_numeric = as.numeric(single_initial_est)
        #remove(single_initial_est)
        XY = grad_prep(X_all, Y_all)
        XX_all = lapply(X_all, function(x){t(x)%*%x/nrow(x)})
        tmp_res = rep(0, fold)

        for(f in 1:fold){
                ans = glasso_no_early_stopping(X=X_all, Y=Y_all, XX=XX_all, XY=XY, Xnorm=Xnorm, lambda1=multi_lam[f,1]/spsz, lambda2=multi_lam[f,2], theta=matrix(initial_numeric,ncol(dose), length(tissue_gene)), verbose = if_verbose)
                tmp_res[f] = ans$avg_train_err
        }
	final.lam = multi_lam[which.min(tmp_res),]
        ans = glasso_no_early_stopping(X=X_all, Y=Y_all, XX=XX_all, XY=XY, Xnorm=Xnorm, lambda1=final.lam[1]/spsz, lambda2=final.lam[2], theta=matrix(initial_numeric,ncol(dose), length(tissue_gene)), verbose = if_verbose)
	coef_temp = ans$est

###############################################################################################################################################################################################
	for (t in 1:length(tissue_gene)) {
		current_tissue_index = tissue_gene[t]
		current_tissue = tissue_list[current_tissue_index]

		if ( all( as.data.frame(coef_temp[,t]) == 0 ) ) {
			summary_temp = data.frame( "tissue" = current_tissue, "gene_id" = map$gene_id[i],
               	                           "func_type" = func_type, "region" = origin_id, fit = "ctimp",
                       	                   "input_snps" = ncol(dose), "nonzero_snps" = 0)
	               	summary[[current_tissue_index]] = rbind(summary[[current_tissue_index]], summary_temp)			
		} else {
			coef_nonzero = as.matrix(coef_temp[,t][ which(coef_temp[,t] != 0) ])
			rownames(coef_nonzero) = colnames(dose)[ which(coef_temp[,t] != 0) ]
			coef_output = data.frame( "tissue" = rep(current_tissue, nrow(coef_nonzero)), "gene_id" = rep(map$gene_id[i], nrow(coef_nonzero)),
                                          "snp" = rownames(coef_nonzero), "weight" = coef_nonzero[,1]) %>% remove_rownames
	                coef[[current_tissue_index]] = rbind(coef[[current_tissue_index]], coef_output)
		
			summary_temp = data.frame( "tissue" = current_tissue, "gene_id" = map$gene_id[i],
                                           "func_type" = func_type, "region" = origin_id, fit = "ctimp",
                                           "input_snps" = ncol(dose), "nonzero_snps" = nrow(coef_nonzero))
			summary[[current_tissue_index]] = rbind(summary[[current_tissue_index]], summary_temp)
		}
		
		filename_summ_result <- paste(dir_path, "/utmost/summary/", current_tissue, "/", current_tissue, "_chr", chr_num, "_GTEx_utmost_summary_", total_file_num, ".", file_num, ".txt", sep="")
                write.table(summary[[current_tissue_index]], filename_summ_result, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

                filename_coef_result <- paste(dir_path, "/utmost/output/", current_tissue, "/", current_tissue, "_chr", chr_num, "_GTEx_utmost_coef_", total_file_num, ".", file_num, ".txt", sep="")
                write.table(coef[[current_tissue_index]], filename_coef_result, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

	}

}

dummy="I am done"
filename_dummy = paste(dir_path, "/utmost/dummy/dummy_testing_chr", chr_num, "_", total_file_num, ".", file_num, ".txt", sep ="" )
write.table(dummy, filename_dummy, quote=FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
