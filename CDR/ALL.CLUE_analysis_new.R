library(data.table)
library(dplyr)
library(matrixStats)
library(stringr)
options(stringsAsFactors=F)

pheno_list = c("ast", "cd", "ecz", "ibd", "ra", "sle", "uc", "vit")
method_list = c("cauchy_combined_mashr", "cauchy_dapg_mashr", "cauchy_predixcan_mashr", 
		"cauchy_combined", "cauchy_dapg_utmost", "cauchy_predixcan_utmost",
		"mashr", "dapg", "multi-omics", "utmost", "predixcan")

##Drug with known indication##
df = data.frame()
for (pheno in pheno_list){
	print(pheno)
	
	##Load in clue result##
	result = as.data.frame(fread(paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/clue_sig_new/total/", pheno, "_all_clue_processed.txt", sep = "") , 
				header = FALSE, na.strings = c("NA", "NaN")))
	cell = result[1,]
	
	##Find method orders##
	meth = result[2, which(cell != "summary")]
        #meth = result[2, c(1,2,which(cell == "summary"))]
	meth = str_replace(meth[1,], paste(pheno, "_", sep =""), "")

	result = result[3:nrow(result), which(cell != "summary")]
	#result = result[3:nrow(result), c(1,2, which(cell == "summary"))]
	result[is.na(result)] <- 0

	##Load in drug indication file##
	drug_ref = read.table( paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/clue/drug_indication/", pheno, "_drug_indication.txt", sep = "") )

	##Identify score of drug with indication##
	score = data.frame()
	for (i in 1:nrow(drug_ref)){
		temp = dplyr::filter(result, grepl(drug_ref[i,1], V2))
		score = rbind(score, temp)
	}
	score$V1 = paste(toupper(pheno), score$V1, sep = "_")

	##Calculate mean/median of score across 9 cell lines of each method##
	for ( method_a in method_list ) {
		temp = cbind("method" = method_a, "trait" = pheno, "name" = score$V1, "name_d" = paste(pheno, score$V2, sep = "_"),
                        "mean" = rowMeans(as.data.frame(sapply(score[, which(meth == method_a)], as.numeric)), na.rm=TRUE),
                        "median" = rowMedians(as.matrix(sapply(score[, which(meth == method_a)], as.numeric)), na.rm=TRUE))
		df = rbind(df, temp)
	}
}
#df1 = df %>% filter(mean != "NaN") %>% filter(median != "NaN")	
#df2 = df1[df1$name %in% (df1 %>% group_by(name) %>% tally() %>% filter(n == 6))$name, ]
write.table(df, "/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/clue_sig_new/result_all_tissues_new.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
#write.table(df, "/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/clue_sig_new/result_all_tissues_summary.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

##Looking for novel drugs##
df = data.frame()
for (pheno in pheno_list){

	final_list = list()
	final_list[[1]] = data.frame()
	final_list[[2]] = data.frame()
	final_list[[3]] = data.frame()
	final_list[[4]] = data.frame()
	final_list[[5]] = data.frame()
	final_list[[6]] = data.frame()
	
	##Load in clue result##
	result = as.data.frame(fread(paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/clue_sig/total/", pheno, "_all_clue_processed.txt", sep = "") ,
                                header = FALSE, na.strings = c("NA", "NaN")))

	result$V1 = paste(toupper(pheno), result$V1, sep = "_")
	
	##Calculate mean/median of score across 9 cell lines of each method##
        temp_cauchy = cbind("method" = "cauchy", "trait" = pheno, "name" = result$V1, "mean" = rowMeans(result[,3:11], na.rm=TRUE), "median" = rowMedians(as.matrix(result[,3:11]), na.rm=TRUE))
        final_list[[1]] = rbind(final_list[[1]], temp_cauchy)
        temp_dapg = cbind("method" = "dapg", "trait" = pheno, "name" = result$V1, "mean" = rowMeans(result[,12:20], na.rm=TRUE), "median" = rowMedians(as.matrix(result[,12:20]), na.rm=TRUE))
	final_list[[2]] = rbind(final_list[[2]], temp_dapg)
        temp_mashr = cbind("method" = "mashr", "trait" = pheno, "name" = result$V1, "mean" = rowMeans(result[,21:29], na.rm=TRUE), "median" = rowMedians(as.matrix(result[,21:29]), na.rm=TRUE))
	final_list[[3]]	= rbind(final_list[[3]], temp_mashr)
        temp_multi = cbind("method" = "multi-omics", "trait" = pheno, "name" = result$V1, "mean" = rowMeans(result[,39:47], na.rm=TRUE), "median" = rowMedians(as.matrix(result[,39:47]), na.rm=TRUE))
	final_list[[4]]	= rbind(final_list[[4]], temp_multi)
        temp_predixcan = cbind("method" = "predixcan", "trait" = pheno, "name" = result$V1, "mean" = rowMeans(result[,48:56], na.rm=TRUE), "median" = rowMedians(as.matrix(result[,48:56]), na.rm=TRUE))
	final_list[[5]]	= rbind(final_list[[5]], temp_predixcan)
        temp_utmost = cbind("method" = "utmost", "trait" = pheno, "name" = result$V1, "mean" = rowMeans(result[,57:65], na.rm=TRUE), "median" = rowMedians(as.matrix(result[,57:65]), na.rm=TRUE))
	final_list[[6]]	= rbind(final_list[[6]], temp_utmost)

	##Identify drug overlaps##
	overlap = intersect(final_list[[6]]$name, intersect(final_list[[5]]$name, intersect(final_list[[1]]$name, intersect(final_list[[2]]$name, intersect(final_list[[3]]$name, final_list[[4]]$name)))))

        for (idx in 1:6){
                final_list[[idx]] = final_list[[idx]][final_list[[idx]]$name %in% overlap, ]
                final_list[[idx]] = final_list[[idx]][order(as.numeric(as.character(final_list[[idx]]$median))), ]
                final_list[[idx]]$Rank = 1:nrow(final_list[[idx]])
		final_list[[idx]] = final_list[[idx]] %>% select(name, median, Rank)
                colnames(final_list[[idx]])[2:3] = c(paste("Score", idx, sep = "_"), paste("Rank", idx, sep = "_"))
        }

	merge = merge(final_list[[1]], final_list[[2]], by.x = "name", by.y = "name")
        merge = merge(merge, final_list[[3]], by.x = "name", by.y = "name")
        merge = merge(merge, final_list[[4]], by.x = "name", by.y = "name")
	merge = merge(merge, final_list[[5]], by.x = "name", by.y = "name")
	merge = merge(merge, final_list[[6]], by.x = "name", by.y = "name")
	df = rbind(df, merge)
}

rownames(df) = df$name
df = df %>% select(-c("name"))

##Determine the sign based on the majority of the 6 methods##
neg_count = rowSums(df %>% select(-c("Rank_1", "Rank_2", "Rank_3", "Rank_4", "Rank_5", "Rank_6")) < -30)
pos_count = rowSums(df %>% select(-c("Rank_1", "Rank_2", "Rank_3", "Rank_4", "Rank_5", "Rank_6")) > 30)
zero_count = rowSums(df %>% select(-c("Rank_1", "Rank_2", "Rank_3", "Rank_4", "Rank_5", "Rank_6")) == 0)
count = cbind(neg_count, pos_count, zero_count)
x = neg_count > pos_count & neg_count > zero_count
#x = neg_count == 6
count = count[x,]

df1 = df[rownames(df) %in% rownames(count),]


