#!/bin/bash

INPUT="<PATH TO INPUT FOLDER>"
OUTPUT="<PATH TO OUTPUT FOLDER>"
BEDTOOLS="<PATH TO BEDTOOLS>"

Rscript PUMICE.compute_weights.R \
	--geno ${INPUT}/GEUVADIS_chr22_genotype_subset.traw \
	--chr 22 \
	--exp ${INPUT}/GEUVADIS_chr22_expression_subset.txt \
	--out ${OUTPUT} \
	--pchic_path ${INPUT}/LCL_chr22_pchic.txt \
	--loop_path ${INPUT}/LCL_chr22_loop.txt \
	--tad_path ${INPUT}/LCL_chr22_tad.txt \
	--domain_path ${INPUT}/LCL_chr22_domain.txt \
	--bedtool_path ${BEDTOOLS} \
	--epi_path ${INPUT}/ENCFF028SGJ_chr22_screen.bed \
	--fold 5
