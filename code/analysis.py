# Variance

list=("RP" "species" "pathways" "taxo" "carbon" "blood" "anthro")
var='oralmsp'
for var in "${list[@]}"
do
splitter ${var} Type
#variance -df1 ${var}TypeFMT -df2 meta
stratify ${var}TypeFMT 'Timepoint'
change ${var}TypeFMTTimepoint
variance ${var}Visitwk6 meta --pval 1
variance ${var}Visitwk12 meta --pval 1
variance ${var}Visitwk26 meta --pval 1
variance ${var}Visityr4 meta --pval 1
done

merge $(ls ../results/*power.tsv| xargs -I{} basename {} .tsv) -f allpower
