#!/bin/bash

FILENUM=$1
MET=$2 ##3d,constant##
TYPE=$3 ##3d:loop,tad,domain,pchic; constant:250,1000##

INPUT="<PATH TO INPUT FOLDER>"
OUTPUT="<PATH TO OUTPUT FOLDER>"
BEDTOOLS="<PATH TO BEDTOOLS>"

Rscript PUMICE.nested_cv.R \
	--geno ${INPUT}/GEUVADIS_chr22_genotype_subset.traw \
	--chr 22 \
	--exp ${INPUT}/GEUVADIS_chr22_expression_subset.txt \
	--out ${OUTPUT} \
	--method ${MET} \
	--type ${TYPE} \
	--window_path ${INPUT}/LCL_chr22_${type}.txt \
	--bedtools_path ${BEDTOOLS} \
	--epi_path ${INPUT}/ENCFF028SGJ_chr22_screen.bed \
	--fold 5 \
	--total_file_num 10 \
	--file_num ${FILENUM}

