#!/bin/bash

INPUT="<PATH TO INPUT FOLDER>"
OUTPUT="<PATH TO OUTPUT FOLDER>"
BEDTOOLS="<PATH TO BEDTOOLS>"

##Running nested cross-validation in 3D genome windows##
for FILENUM in {1..10};
do
	for TYPE in "loop" "tad" "domain" "pchic";
	do

	Rscript PUMICE.nested_cv.R \
		--geno ${INPUT}/GEUVADIS_chr22_genotype_subset.traw \
		--chr 22 \
		--exp ${INPUT}/GEUVADIS_chr22_expression_subset.txt \
		--out ${OUTPUT} \
		--method 3d \
		--type ${TYPE} \
		--window_path ${INPUT}/LCL_chr22_${type}.txt \
		--bedtools_path ${BEDTOOLS} \
		--epi_path ${INPUT}/ENCFF028SGJ_chr22_screen.bed \
		--fold 5 \
		--total_file_num 10 \
		--file_num ${FILENUM}
	done
done

##Running nested cross-validation in constant windows##
for FILENUM in {1..10};
do
	for TYPE in "250" "1000";
	do

	Rscript PUMICE.nested_cv.R \
		--geno ${INPUT}/GEUVADIS_chr22_genotype_subset.traw \
		--chr 22 \
		--exp ${INPUT}/GEUVADIS_chr22_expression_subset.txt \
		--out ${OUTPUT} \
		--method constant \
		--type ${TYPE} \
		--bedtools_path ${BEDTOOLS} \
		--epi_path ${INPUT}/ENCFF028SGJ_chr22_screen.bed \
		--fold 5 \
		--total_file_num 10 \
		--file_num ${FILENUM}
	done
done

