#!/bin/bash

cat <<EOS | qsub -
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=47:59:59
#PBS -l feature=rhel7
#PBS -l pmem=40gb
#PBS -m n
#PBS -N g.total
#PBS -A open
#PBS -e /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/log/log_total.e
#PBS -o /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/log/log_total.o

module load anaconda3/2020.07
source activate /storage/home/cxk502/work/conda_base

Rscript /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/PUMICE.compute_weights.R \
	--geno /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/GEUVADIS_chr22_genotype_subset.traw \
	--chr 22 \
	--exp /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/GEUVADIS_chr22_expression_subset.txt \
	--out /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/output \
	--pchic_path /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/LCL_chr22_pchic.txt \
	--loop_path /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/LCL_chr22_loop.txt \
	--tad_path /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/LCL_chr22_tad.txt \
	--domain_path /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/LCL_chr22_domain.txt \
	--bedtool_path /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/bedtools \
	--epi_path /gpfs/group/dxl46/default/private/poom/GTEx/master_code/all/github/ENCFF028SGJ_chr22_screen.bed \
	--fold 5

EOS
