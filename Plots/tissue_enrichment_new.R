library(dplyr)
library(ggplot2)
library(reshape2)

pheno_list = c("ast", "ecz", "cs", "cd", "ibd", "pbc", "ra", "sle", "uc", "vit", "covid",
               "ai_ng", "ad", "als", "an", "anx_cc", "anx_fs", "adhd", "asd", "bip", "cpd_ng", "ds", "dpw_ng", "edu", "epl", "ext", "neurotic", "scz", "sc_ng", "si_ng", "swb",
               "hdl", "ldl", "tc", "tg", "cad", "mi",
               "gluc2", "homab", "fg", "fgt", "fg_female", "fg_male", "fi", "fi_female", "fi_male", "fpi", "hba1c", "homair", "isi", "t2d",
               "gest_dur", "gwg", "blength", "bweight", "mbweight", "bmi", "childbmi", "chob", "hc", "height", "hip", "hipadjbmi", "menarche", "menopause",
               "ts", "wc", "wcadjbmi", "whr", "whradjbmi",
               "fn_bmd", "fa_bmd", "ls_bmd",
               "gout", "amd", "su")

pheno_name = c("Asthma", "Atopic dermatitis", "Churg-Strauss syndorme", "Crohn's disease", "IBD",
               "Primary biliary cirrhosis", "Rheumatoid arthritis", "SLE", "Ulcerative colitis", "Vitiligo", "COVID19", 
               "AI-smoking", "Alzheimer's disease", "ALS", "Anorexia nervosa", "Anxiety-cc", "Anxiety-fs", "ADHD", "ASD",
               "Bipolar disorder", "Cigarettes/day", "Depressive symptoms", "Drinks/week", "Educational level", "Epilepsy", "Extraversion", "Neuroticism",
               "Schizophrenia", "Smoking cessation", "Smoking initiation", "Subjective well being", "HDL", "LDL", "Total cholesterol", "Triglyceride",
               "Coronary artery disease", "Myocardial infarction", "2hr glucose", "Beta-cell function", "Fasting glucose",
               "Fasting glucose over time", "Fasting glucose-female", "Fasting glucose-male", "Fasting insulin",
               "Fasting insulin-female", "Fasting insulin-male", "Fasting pro-insulin", "HbA1c", "Insulin resistance",
               "Insulin sensitivity", "Type 2 diabetes", "Gestational duration", "Gestational weight gain", "Birth length", "Birth weight", "Birth weight-maternal",
               "BMI", "Childhood BMI", "Childhood obesity", 
               "Head circumference", "Height", "Hip circumference", "Hip-adjusted BMI", "Mecharche", "Menopause",
               "Tanner stage", "Waist circumference", "WC-adjusted BMI", "Waist-to-hip-ratio",
               "WHR-adjusted BMI", "Femoral neck BMD", "Forearm BMD", "Lumbar spine BMD",
               "Gout", "Macular degeneration", "Serum urate")

tissue_list = c("Muscle_Skeletal", "Skin_Sun_Exposed_Lower_leg", "Thyroid", "Adipose_Subcutaneous", "Lung", "Artery_Tibial",
                "Whole_Blood", "Nerve_Tibial", "Esophagus_Mucosa", "Skin_Not_Sun_Exposed_Suprapubic", "Esophagus_Muscularis", "Adipose_Visceral_Omentum",
                "Cells_Transformed_fibroblasts", "Artery_Aorta", "Heart_Left_Ventricle", "Heart_Atrial_Appendage", "Breast_Mammary_Tissue", "Colon_Transverse",
                "Stomach", "Testis", "Esophagus_Gastroesophageal_Junction", "Colon_Sigmoid", "Pancreas", "Pituitary",
                "Adrenal_Gland", "Brain_Cerebellum", "Liver", "Brain_Caudate_basal_ganglia", "Artery_Coronary", "Brain_Cortex",
                "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Spleen", "Prostate", "Brain_Frontal_Cortex_BA9", "Brain_Anterior_cingulate_cortex_BA24",
                "Small_Intestine_Terminal_Ileum", "Brain_Hippocampus", "Brain_Putamen_basal_ganglia", "Brain_Hypothalamus", "Ovary", "Cells_EBV-transformed_lymphocytes",
                "Vagina", "Uterus", "Brain_Amygdala", "Brain_Spinal_cord_cervical_c-1", "Minor_Salivary_Gland", "Brain_Substantia_nigra")

method = "cauchy"
load(paste("/Users/chachritkhun/Desktop/TWAS_files/tissue_enrichment/", method, "_residual.RDat", sep = ""))
residual = adjp_residual
residual$tissue = rownames(residual)

x = melt(residual, id = c("tissue"))
x$variable = factor(x$variable, levels = pheno_list)
levels(x$variable) <- pheno_name

tiff("/Users/chachritkhun/Desktop/plot.tiff", units="in",  width=16, height=8, res=200)
ggplot(x, aes(variable, tissue)) +
  geom_point(pch=15, aes(color=value, size=abs(value))) +
  scale_size_continuous(range = c(0, 4)) +
  #scale_color_manual(values=jcolors) +
  scale_colour_gradient2(low = "#306EAF", mid = "#F8F8F8", high = "#B72333", limits=c(-8, 8)) +
  geom_vline(xintercept=seq(0, ncol(residual) - 1, 1)+.5,color="black") +
  geom_hline(yintercept=seq(0,nrow(residual), 1)+.5,color="black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.45, hjust=1)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  ggtitle("PUMICE+(UTMOST)") +
  theme(plot.title = element_text(face="bold", size = 16, family = "Helvetica", hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 8, family = "Helvetica")) +
  theme(axis.text.y = element_text(size = 8, family = "Helvetica")) 
dev.off()