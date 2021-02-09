#!/bin/bash

file_num=$1
met=$2 ##3d,constant##
type=$3 ##3d:loop,tad,domain,pchic; constant:250,1000##

cat <<EOS | qsub -
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=47:59:59
#PBS -l feature=rhel7
#PBS -l pmem=40gb
#PBS -m n
#PBS -N g."$type"."$file_num"
#PBS -A open
#PBS -e /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/log/log_"$met"_"$type"_"$file_num".e
#PBS -o /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/log/log_"$met"_"$type"_"$file_num".o

module load anaconda3/2020.07
source activate /storage/home/cxk502/work/conda_base

Rscript /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/PUMICE.nested_cv.R \
	--geno /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/GEUVADIS_chr22_genotype_subset.traw \
	--chr 22 \
	--exp /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/GEUVADIS_chr22_expression_subset.txt \
	--out /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/output \
	--method ${met} \
	--type ${type} \
	--window_path /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/LCL_chr22_${type}.txt \
	--bedtool_path /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/bedtools \
	--epi_path /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/ENCFF028SGJ_chr22_screen.bed \
	--fold 5 \
	--total_file_num 10 \
	--file_num ${file_num}

EOS
