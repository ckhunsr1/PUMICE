#!/bin/bash

file_num=$1
met=$2 ##3d,constant##
type=$3 ##3d:loop,tad,domain,pchic; constant:250,1000##

Rscript PUMICE.nested_cv.R \
	--geno example_input/GEUVADIS_chr22_genotype_subset.traw \
	--chr 22 \
	--exp example_input/GEUVADIS_chr22_expression_subset.txt \
	--out "<PATH TO OUTPUT FOLDER>" \
	--method ${met} \
	--type ${type} \
	--window_path example_input/LCL_chr22_${type}.txt \
	--bedtool_path "<PATH TO BEDTOOLS>" \
	--epi_path example_input/ENCFF028SGJ_chr22_screen.bed \
	--fold 5 \
	--total_file_num 10 \
	--file_num ${file_num}

