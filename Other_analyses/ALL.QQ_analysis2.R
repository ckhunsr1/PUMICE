library(dplyr)
library(stringr)
library(data.table)

##Import reference file (protein coding genes)##
#ref = read.table("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/id_position.txt", header = TRUE)
ref = read.table("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/id_conversion/ensembl_position.txt", header = TRUE)
mhc = ref$chr == 6 & ref$start > 26e6 & ref$end < 34e6
ref = ref[!mhc, ]

#method_list = c("multi-omics", "predixcan", "utmost", "dapg", "mashr", "combined", "combined_mashr")
#method_list = c("dapg_mashr", "dapg_utmost", "predixcan_mashr", "predixcan_utmost")
method_list = c("predixcan_utmost")

##Import p-values##
for (method in method_list){
	print(method)
	if ( method != "combined" && method != "combined_mashr" && method != "dapg_mashr" && method != "dapg_utmost" && method != "predixcan_mashr" && method != "predixcan_utmost" ){
		df = as.data.frame(fread(paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method, "_QQ_analysis.txt", sep = "")))
		df = df[!(is.na(df$V4)), ]
		df = df[df$V3 %in% ref$ensembl_gene_id, ]
		df$chisq = qchisq(df$V4, 1, lower.tail=FALSE)
		lambda = median(df$chisq)/qchisq(0.5,1)
		print(lambda)
		df$newchisq = df$chisq/lambda
		df$TWAS.P.new = pchisq(df$newchisq, df =1, lower.tail=FALSE)
		write.table(df %>% select(TWAS.P.new), paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method, "_QQ_analysis2.txt", sep =""), 
		    	    col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
	} else {
		df = as.data.frame(fread(paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method, "_QQ_analysis.txt", sep = "")))

		##Process minp##
                df1 = df[!(is.na(df$V4)), ] %>% select(V1, V2, V3, V4)
		df1 = df1[df1$V3 %in% ref$ensembl_gene_id, ]
		df1$chisq = qchisq(df1$V4, 1, lower.tail=FALSE)
                lambda = median(df1$chisq, na.rm = TRUE)/qchisq(0.5,1)
                print(lambda)
                df1$newchisq = df1$chisq/lambda
                df1$TWAS.P.new = pchisq(df1$newchisq, df =1, lower.tail=FALSE)
                write.table(df1 %>% select(TWAS.P.new), paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method, "_minp_QQ_analysis2.txt", sep =""),
                            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

		##Process cauchy##
		df2 = df[!(is.na(df$V5)), ] %>% select(V1, V2, V3, V5)
		df2 = df2[df2$V3 %in% ref$ensembl_gene_id, ]
		df2$chisq = qchisq(df2$V5, 1, lower.tail=FALSE)
                lambda = median(df2$chisq, na.rm = TRUE)/qchisq(0.5,1)
                print(lambda)
                df2$newchisq = df2$chisq/lambda
                df2$TWAS.P.new = pchisq(df2$newchisq, df =1, lower.tail=FALSE)
                write.table(df2 %>% select(TWAS.P.new), paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method, "_cauchy_QQ_analysis2.txt", sep =""),
                            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

	}

}
