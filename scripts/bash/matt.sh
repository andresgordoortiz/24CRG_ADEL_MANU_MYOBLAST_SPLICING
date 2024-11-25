mkdir results/tables/matt_out
sed 's/,/\t/g' results/tables/final_exon_list_unsorted_20241122.csv > results/tables/matt_out/exon_list.tab
matt add_cols results/tables/matt_out/exon_list.tab SCAFF START END STRAND < <(awk -F'\t' '{split($3,a,":"); split(a[2],b,"-"); print a[1], b[1], b[2], "+"}' results/tables/matt_out/exon_list.tab)