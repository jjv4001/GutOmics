#nohup sh run_cell_type_gwas.sh 2>&1 > log_run_cellType.log

### run for cell type first
path0=atac_gwas
path1=../peaks

alias cxrun='/gpfs/commons/groups/zhu_lab/czhu/scripts/scripts/01.general/cxrun_ne1'

mkdir $path0
cd $path0

##run ld_score
mkdir ld_score
gs=/gpfs/commons/groups/zhu_lab/czhu/01.project_data/archieve/07.prev_Renlab/09.SNARE-seq/GWAS/resources/hg38.chrom.sizes
plink=/gpfs/commons/groups/zhu_lab/czhu/01.project_data/archieve/07.prev_Renlab/09.SNARE-seq/GWAS/resources/GRCh38/plink_files


for i in ${path1}/*.bed
do
  k=$(basename $i .bed)
  for j in `seq 1 22`
  do
    cxrun 1 1 sh ../run_cell_type_individual.sh $k $i $j
  done
done

