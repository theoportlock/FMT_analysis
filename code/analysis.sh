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


### Correlation
filter -m 40 msp
corr -df mspfilter
describe mspfiltercorr --corr
filter -c 'sig' -lt 0.05 mspcorr
scale log msp
scale standard msplog
scale minmax msplogstandard
plot scatter msplogstandardminmax --logy --logx '{"x":"Streptococcus_salivarius","y":"Veillonella_dispar"}'
plot circos allcorr
plot heatmap allcorrsig


### Prediction
stratify pathwaysVisityr4 -l Treatment
predict classifier pathwaysVisityr4Treatment
explain SHAP_bin pathwaysVisityr4Treatment

list=("RP" "species" "pathways" "taxo" "carbon" "blood" "anthro")
for var in "${list[@]}"
do
stratify ${var}VisitBL -l Treatment
predict classifier ${var}VisitBLTreatment
stratify ${var}Visitwk6 -l Treatment
predict classifier ${var}Visitwk6Treatment
stratify ${var}Visitwk12 -l Treatment
predict classifier ${var}Visitwk12Treatment
stratify ${var}Visitwk26 -l Treatment
predict classifier ${var}Visitwk26Treatment
stratify ${var}Visityr4 -l Treatment
predict classifier ${var}Visityr4Treatment
done

merge $(ls ../results/*filter.tsv| xargs -I{} basename {} .tsv) -f allstandard
stratify metabstandardeegstandardbayleysstandardmspstandardanthrostandardpathwaysstandard -l Condition
predict classifier metabstandardeegstandardbayleysstandardmspstandardanthrostandardpathwaysstandardCondition
plot aucroc alldataStandardCondition
plot aucroc alldataStandardEF

# SHAP
explain SHAP_bin metabstandardeegstandardbayleysstandardmspstandardanthrostandardpathwaysstandardCondition
explain SHAP_interact metabstandardeegstandardbayleysstandardmspstandardanthrostandardpathwaysstandardCondition
expalain alldataStandardpredictEF
plot breakdown alldataStandardpredictCondition
plot breakdown alldataStandardpredictEF

### Network
list=("metab" "species" "pathways" "blood" "anthro")
for var in "${list[@]}"
do
filter ${var} --min_unique 5
stratify ${var}filter -l Treatment
change ${var}filterTreatment
filter -c 'mww_sig(FMT/Placebo)' -lt 0.05 ${var}filterTreatmentchange
filter -c 'log2(FMT/Placebo)' -absgt 1 ${var}filterTreatmentchangefilter
filter ${var}filter -fdf ${var}filterTreatmentchangefilter -fdfx 1
done

merge metabfilterfilter speciesfilterfilter pathwaysfilterfilter bloodfilter anthrofilterfilter
corr metabfilterfilterspeciesfilterfilterpathwaysfilterfilterbloodfilteranthrofilterfilter
filter -c 'sig' -lt 0.05 metabfilterfilterspeciesfilterfilterpathwaysfilterfilterbloodfilteranthrofilterfiltercorr
filter -c 'rho' -absgt 0.5 metabfilterfilterspeciesfilterfilterpathwaysfilterfilterbloodfilteranthrofilterfiltercorrfilter

merge metabfilterTreatmentchangefilter speciesfilterTreatmentchangefilter pathwaysfilterTreatmentchangefilter bloodfilterTreatmentchange anthrofilterTreatmentchange -a

### Yr4
#speciesVisityr4TreatmentFMTRespondermult
list=("metab" "species" "pathways" "blood" "anthro")
t='yr4'
val='Treatment'
val='Responder'
for var in "${list[@]}"
do
splitter ${var} Visit
filter ${var}Visit${t} --min_unique 2
stratify ${var}Visit${t}filter -l ${val}
change ${var}Visit${t}filter${val}
#filter -c 'mww_sig(FMT/Placebo)' -lt 0.05 ${var}Visit${t}filter${val}change
filter -c 'mww_sig(NonResponder/Responder)' -lt 0.05 ${var}Visit${t}filter${val}change
#filter -c 'log2(FMT/Placebo)' -absgt 1 ${var}Visit${t}filter${val}changefilter
filter -c 'log2(NonResponder/Responder)' -absgt 1 ${var}Visit${t}filter${val}changefilter
filter ${var}Visit${t}filter -fdf ${var}Visit${t}filter${val}changefilter -fdfx 1
done

merge metabVisit${t}filterfilter speciesVisit${t}filterfilter pathwaysVisit${t}filterfilter bloodVisit${t}filter anthroVisit${t}filterfilter
corr metabVisit${t}filterfilterspeciesVisit${t}filterfilterpathwaysVisit${t}filterfilterbloodVisit${t}filteranthroVisit${t}filterfilter
filter -c 'sig' -lt 0.05 metabVisit${t}filterfilterspeciesVisit${t}filterfilterpathwaysVisit${t}filterfilterbloodVisit${t}filteranthroVisit${t}filterfiltercorr
corr metabfilterfilterspeciesfilterfilterpathwaysfilterfilterbloodfilteranthrofilterfilter
filter -c 'sig' -lt 0.05 metabfilterfilterspeciesfilterfilterpathwaysfilterfilterbloodfilteranthrofilterfiltercorr
filter -c 'rho' -absgt 0.5 metabfilterfilterspeciesfilterfilterpathwaysfilterfilterbloodfilteranthrofilterfiltercorrfilter

merge metabVisit${t}filter${val}changefilter speciesVisit${t}filter${val}changefilter pathwaysVisit${t}filter${val}changefilter bloodVisit${t}filter${val}changefilter anthroVisit${t}filter${val}changefilter -a

# mechanism
plot scatter mspfilterfiltermetabfilterfilterpathwaysfilterfiltereegfilterfilternirsfilterfilteranthrominmaxfilterfilterbayleysfilterfilter '{"x":"ExpressiveCommunicationScore", "y":"P108-PWY: pyruvate fermentation to propanoate I", "figsize":(2,2), "scatter_kws":{"s":0.01}}'

stratify ${var}Visit${t}filter -l ${val}
plot box speciesVisit${t} '{"y":"s__Bacteroides_thetaiotaomicron","x":"figsize":(1,2)}'
plot box speciesVisityr4Responder '{"y":"s__Bacteroides_thetaiotaomicron","figsize":(1,2)}'

plot box speciesVisityr4TreatmentFMTResponder '{"y":"s__Bacteroides_thetaiotaomicron"}'
scale mult speciesVisityr4TreatmentFMTResponder
scale speciesVisityr4TreatmentFMTResponder
plot box speciesVisityr4TreatmentFMTRespondermult --logy '{"y":"s__Bacteroides_thetaiotaomicron"}'
plot box speciesVisityr4TreatmentFMTRespondermult --logy '{"y":"s__Prevotella_copri"}'


### Yr4
#speciesVisityr4TreatmentFMTRespondermult
list=("metab" "species" "pathways" "blood" "anthro")
t='yr4'
val='Responder'
for var in "${list[@]}"
do
splitter ${var} Visit
splitter ${var}Visit${t} Treatment
stratify ${var}Visit${t}TreatmentFMT -l Responder
change ${var}Visit${t}TreatmentFMTResponder
filter -c 'mww_sig(NonResponder/Responder)' -lt 0.05 ${var}Visit${t}TreatmentFMTResponderchange
filter -c 'log2(NonResponder/Responder)' -absgt 1 ${var}Visit${t}TreatmentFMTResponderchangefilter
filter ${var}Visit${t}TreatmentFMT -fdf ${var}Visit${t}TreatmentFMTResponderchangefilter -fdfx 1
done

merge metabVisityr4TreatmentFMTfilter speciesVisityr4TreatmentFMTfilter pathwaysVisityr4TreatmentFMTfilter bloodVisityr4TreatmentFMTfilter anthroVisityr4TreatmentFMTfilter
corr metabVisityr4TreatmentFMTfilterspeciesVisityr4TreatmentFMTfilterpathwaysVisityr4TreatmentFMTfilterbloodVisityr4TreatmentFMTfilteranthroVisityr4TreatmentFMTfilter
filter -c 'sig' -lt 0.05 metabVisityr4TreatmentFMTfilterspeciesVisityr4TreatmentFMTfilterpathwaysVisityr4TreatmentFMTfilterbloodVisityr4TreatmentFMTfilteranthroVisityr4TreatmentFMTfiltercorr
#filter -c 'rho' -absgt 0.5 metabfilterfilterspeciesfilterfilterpathwaysfilterfilterbloodfilteranthrofilterfiltercorrfilter

merge metabVisit${t}TreatmentFMTResponderchangefilter speciesVisit${t}TreatmentFMTResponderchangefilter pathwaysVisit${t}TreatmentFMTResponderchangefilter anthroVisit${t}TreatmentFMTResponderchangefilter -a
