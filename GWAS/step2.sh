#nohup sh run_cell_type_gwas.sh 2>&1 > log_run_cellType.log

### run for cell type first
path0=atac_gwas
path1=../peaks
alias cxrun='/gpfs/commons/groups/zhu_lab/czhu/scripts/scripts/01.general/cxrun_ne1'
mkdir $path0
cd $path0

##run ld_score
# mkdir ld_score
gs=/gpfs/commons/groups/zhu_lab/czhu/01.project_data/archieve/07.prev_Renlab/09.SNARE-seq/GWAS/resources/hg38.chrom.sizes
plink=/gpfs/commons/groups/zhu_lab/czhu/01.project_data/archieve/07.prev_Renlab/09.SNARE-seq/GWAS/resources/GRCh38/plink_files


# for i in ${path1}/*hg38.bed
# do
#   k=$(basename $i .bed)
#   cxrun 1 1 sh ../run_cell_type_individual.sh $k $i
# done



# ####

resources=/gpfs/commons/groups/zhu_lab/czhu/01.project_data/archieve/07.prev_Renlab/09.SNARE-seq/GWAS/resources

# for i in $(ls ld_score/*22.annot.gz | sed 's/.22.annot.gz//g')
# do
#   k=$(basename $i)
#   #for j in `seq 1 22` 
#   for j in `seq 1 22` 
#   do
#     cxrun 1 1 python /gpfs/commons/groups/zhu_lab/czhu/01.project_data/07.prev_Renlab/09.SNARE-seq/GWAS/ldsc/ldsc/ldsc.py --l2 --bfile ${resources}/GRCh38/plink_files/1000G.EUR.hg38.${j} --ld-wind-cm 1 --annot ld_score/${k}.${j}.annot.gz --thin-annot --out ld_score/${k}.${j} &
#   done
#   wait
# done

########################################################################################################################



# # ## run ldsc_test
mkdir results

ls ${path1}/*.bed | grep -v "merged" | while read i
do
  j=$(basename $i .bed)
  printf "${j}\tld_score/${j}.,ld_score/merged.\n"
done > merged.ldcts


cts_name=merged
sumstats=/gpfs/commons/groups/zhu_lab/czhu/01.project_data/archieve/07.prev_Renlab/09.SNARE-seq/GWAS/sumstats
resources=/gpfs/commons/groups/zhu_lab/czhu/01.project_data/archieve/07.prev_Renlab/09.SNARE-seq/GWAS/resources

ls $sumstats | grep -v merged | while read name
do
  i=$(basename $name .sumstats.gz)
  cxrun 1 1 python /gpfs/commons/groups/zhu_lab/czhu/01.project_data/archieve/07.prev_Renlab/09.SNARE-seq/GWAS/ldsc/ldsc/ldsc.py --h2-cts $sumstats/$i.sumstats.gz --ref-ld-chr $resources/GRCh38/baseline_v1.2/baseline. --out results/${i}.$cts_name --ref-ld-chr-cts $cts_name.ldcts --w-ld-chr $resources/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. ### /media/NAS1/Jinghui_BICCN/GWAS/w_hm3
      
done
wait
