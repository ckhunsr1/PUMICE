library(dplyr);
library(IRanges);
library(GenomicRanges);

tissue_list = c("Muscle_Skeletal", "Skin_Sun_Exposed_Lower_leg", "Thyroid", "Adipose_Subcutaneous", "Lung", "Artery_Tibial",
                        "Whole_Blood", "Nerve_Tibial", "Esophagus_Mucosa", "Skin_Not_Sun_Exposed_Suprapubic", "Esophagus_Muscularis", "Adipose_Visceral_Omentum",
                        "Cells_Transformed_fibroblasts", "Artery_Aorta", "Heart_Left_Ventricle", "Heart_Atrial_Appendage", "Breast_Mammary_Tissue", "Colon_Transverse",
                        "Stomach", "Testis", "Esophagus_Gastroesophageal_Junction", "Colon_Sigmoid", "Pancreas", "Pituitary",
                        "Adrenal_Gland", "Brain_Cerebellum", "Liver", "Brain_Caudate_basal_ganglia", "Artery_Coronary", "Brain_Cortex",
                        "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Spleen", "Prostate", "Brain_Frontal_Cortex_BA9", "Brain_Anterior_cingulate_cortex_BA24",
                        "Small_Intestine_Terminal_Ileum", "Brain_Hippocampus", "Brain_Putamen_basal_ganglia", "Brain_Hypothalamus", "Ovary", "Cells_EBV-transformed_lymphocytes",
                        "Vagina", "Uterus", "Brain_Amygdala", "Brain_Spinal_cord_cervical_c-1", "Minor_Salivary_Gland", "Brain_Substantia_nigra")

pheno_list = setdiff(dir("/gpfs/group/dxl46/default/private/poom/espresso/sumstats_clean_new"), "check.R")

#method_list = c("multi-omics", "predixcan", "utmost", "dapg", "mashr", "combined", "combined_mashr")
method_list = c("dapg_mashr", "dapg_utmost", "predixcan_mashr", "predixcan_utmost")

all_count = read.table("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/iv_all_count.txt", sep = "\t", header = TRUE)
all_count = apply(all_count[,2:4], 1, max)

##Import gene list and position##
ref = read.table("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/id_conversion/ensembl_position.txt", header=TRUE)
colnames(ref)[1] = "gene_id"
ref = ref %>% select(gene_id, start, end, chr)
ref$chr = as.numeric(as.character(ref$chr))
ref = ref[!(is.na(ref$chr)), ]
ref$gene_id = as.character(ref$gene_id)

##Import result##
result = list()

for (idx in 1:length(method_list)){
	method = method_list[idx]
	result_temp = data.frame()

	if ( method != "combined" && method != "combined_mashr" && method != "dapg_mashr" && method != "dapg_utmost" && method != "predixcan_mashr" && method != "predixcan_utmost" ){
		for (pheno in pheno_list){
			print(pheno)

			##Import GWAS clumping result##
	
			tryCatch({	
			##Import index variant##
			path = paste("/gpfs/group/dxl46/default/private/poom/FUSION/summary/clump/", pheno, "/index_gwas_allchr.txt", sep = "")
			id_var = read.table(path, header = FALSE)
			colnames(id_var) = c("id", "start", "end")
			id_var$id = as.character(id_var$id)
			id_var$chr = sapply(strsplit(id_var$id, "_"), function(x){as.numeric(as.character(x[1]))})
	
			##Identified genes at known loci##
			gene_range = makeGRangesFromDataFrame(ref, seqnames.field=c("chr"),start.field="start", end.field=c("end"), ignore.strand=TRUE)
			loci_range = makeGRangesFromDataFrame(id_var, seqnames.field=c("chr"),start.field="start", end.field=c("end"), ignore.strand=TRUE)
			fo = findOverlaps(query=gene_range, subject=loci_range, type = "within")
			known_genes = as.character(unique((as.data.frame(cbind(id_var[subjectHits(fo), 1], ref[queryHits(fo), 1])))$V2))
			},
			error = function(cond) {
			known_genes = c()
			known_genes
			})

			sig_gene = c()
			for (tissue in tissue_list){
				print(tissue)		
				path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/twas_assoc_new/", method, "/", pheno, "/", tissue, "/output", sep = "")
				filename = dir(path)
			
				df = data.frame()
				for (i in 1:length(filename)){
					temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE) %>% select(ID, TWAS.P)
					df = rbind(df, temp)
				}

				df = df %>% filter(TWAS.P < 0.05/all_count[which(tissue_list == tissue)])
				sig_gene_temp = as.character(df$ID)
				sig_gene = unique(c(sig_gene, sig_gene_temp))
			}
			uni_gene = setdiff(sig_gene, known_genes)
			temp = cbind("pheno" = pheno, "total_count" = length(sig_gene), "novel_count" = length(uni_gene))   
			result_temp = rbind(result_temp, temp)
		}
		result[[idx]] = result_temp

	} else {
		for (pheno in pheno_list){
                        print(pheno)

                        ##Import GWAS clumping result##

                        tryCatch({
                        ##Import index variant##
                        path = paste("/gpfs/group/dxl46/default/private/poom/FUSION/summary/clump/", pheno, "/index_gwas_allchr.txt", sep = "")
                        id_var = read.table(path, header = FALSE)
                        colnames(id_var) = c("id", "start", "end")
                        id_var$id = as.character(id_var$id)
                        id_var$chr = sapply(strsplit(id_var$id, "_"), function(x){as.numeric(as.character(x[1]))})

                        ##Identified genes at known loci##
                        gene_range = makeGRangesFromDataFrame(ref, seqnames.field=c("chr"),start.field="start", end.field=c("end"), ignore.strand=TRUE)
                        loci_range = makeGRangesFromDataFrame(id_var, seqnames.field=c("chr"),start.field="start", end.field=c("end"), ignore.strand=TRUE)
                        fo = findOverlaps(query=gene_range, subject=loci_range, type = "within")
                        known_genes = as.character(unique((as.data.frame(cbind(id_var[subjectHits(fo), 1], ref[queryHits(fo), 1])))$V2))
                        },
                        error = function(cond) {
                        known_genes = c()
                        known_genes
                        })

                        sig_gene1 = c()
			sig_gene2 = c()
			for (tissue in tissue_list){
                                print(tissue)
                                path = paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/twas_assoc_new/", method, "/", pheno, "/", tissue, "/output", sep = "")
                                filename = dir(path)

                                df = data.frame()
                                for (i in 1:length(filename)){
                                        temp = read.table(paste(path, "/", filename[i], sep = ""), header = TRUE) %>% select(gene.vec, pval.minp, pval.cauchy)
                                        df = rbind(df, temp)
                                }

				threshold = all_count[which(tissue_list == tissue)]
                                df1 = df %>% filter(pval.minp < 0.05/threshold)
                                sig_gene1_temp = as.character(df1$gene.vec)
                                sig_gene1 = unique(c(sig_gene1, sig_gene1_temp))

				df2 = df %>% filter(pval.cauchy < 0.05/threshold)
                                sig_gene2_temp = as.character(df2$gene.vec)
                                sig_gene2 = unique(c(sig_gene2, sig_gene2_temp))
                        }
                        uni_gene1 = setdiff(sig_gene1, known_genes)
			uni_gene2 = setdiff(sig_gene2, known_genes)
                        temp = cbind("pheno" = pheno, "total_count_minp" = length(sig_gene1), "novel_count_minp" = length(uni_gene1), 
					"total_count_cauchy" = length(sig_gene2), "novel_count_cauchy" = length(uni_gene2))
                        result_temp = rbind(result_temp, temp)
                }
                result[[idx]] = result_temp
	}
}

write.table(result[[1]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[1], "_gwas_loci_count.txt", sep =""), 
	    col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(result[[2]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[2], "_gwas_loci_count.txt", sep =""),  
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(result[[3]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[3], "_gwas_loci_count.txt", sep =""),  
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(result[[4]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[4], "_gwas_loci_count.txt", sep =""),
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
#write.table(result[[5]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[5], "_gwas_loci_count.txt", sep =""),
#            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
#write.table(result[[6]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[6], "_gwas_loci_count.txt", sep =""),
#            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
#write.table(result[[7]], paste("/gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/analyses/", method_list[7], "_gwas_loci_count.txt", sep =""),
#            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
