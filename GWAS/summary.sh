for set in merged
do
  for i in results/*.merged.cell_type_results.txt
  do
    j=$(basename $i .${set}.cell_type_results.txt)
    awk -v j=$j -v s=$set 'NR!=1{print $0"\t"j"\t"s}' $i
  done
done > res.summary.tsv
