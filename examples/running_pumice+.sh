#!/bin/bash

INPUT="<PATH TO INPUT FOLDER>"
OUTPUT="<PATH TO OUTPUT DIRECTORY AND FILENAME>"
BEDTOOLS="<PATH TO BEDTOOLS>"

Rscript PUMICE+.association_test.R \
	--geno ${INPUT}/EUR.chr22.final \
	--chr 22 \
	--gwas ${INPUT}/gwas_ss_chr22.txt \
  --pumice_weight ${INPUT}/PUMICE_Cells_EBV-transformed_lymphocytes.db \
  --utmost_weight ${INPUT}/UTMOST_Cells_EBV-transformed_lymphocytes.db \
  --out ${OUTPUT}
