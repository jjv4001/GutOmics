k=$1
i=$2
j=$3
resources=/gpfs/commons/groups/zhu_lab/czhu/01.project_data/archieve/07.prev_Renlab/09.SNARE-seq/GWAS/resources
gs=/gpfs/commons/groups/zhu_lab/czhu/01.project_data/archieve/07.prev_Renlab/09.SNARE-seq/GWAS/resources/hg38.chrom.sizes
plink=/gpfs/commons/groups/zhu_lab/czhu/01.project_data/archieve/07.prev_Renlab/09.SNARE-seq/GWAS/resources/GRCh38/plink_files

python /gpfs/commons/groups/zhu_lab/czhu/01.project_data/archieve/07.prev_Renlab/09.SNARE-seq/GWAS/ldsc/ldsc/make_annot.py --bed-file $i --bimfile ${plink}/1000G.EUR.hg38.$j.bim --annot-file ld_score/${k}.${j}.annot.gz
python /gpfs/commons/groups/zhu_lab/czhu/01.project_data/archieve/07.prev_Renlab/09.SNARE-seq/GWAS/ldsc/ldsc/ldsc.py --l2 --bfile ${resources}/GRCh38/plink_files/1000G.EUR.hg38.${j} --ld-wind-cm 1 --annot ld_score/${k}.${j}.annot.gz --thin-annot --out ld_score/${k}.${j}
