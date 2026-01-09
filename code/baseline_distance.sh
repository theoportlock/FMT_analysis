#!/bin/bash
source env.sh

#input=results/beta_relabund_analysis/beta_filter.tsv
#meta=results/work/timemeta.tsv
#output=results/beta_relabund_analysis
input=$1
meta=$2
output=$3

mkdir -p $output

join.py $input \
	$meta \
	--left_on source \
	--right_on sampleID \
	-o $output/beta_source.tsv

join.py $output/beta_source.tsv \
	$meta \
	--left_on target \
	--right_on sampleID \
	-o $output/beta_source_target.tsv

filter.py \
	$output/beta_source_target.tsv \
	-q 'timepoint_x == "PCV2"' \
	-o $output/beta_source_target_baseline.tsv

filter.py \
	$output/beta_source_target_baseline.tsv \
	-q 'subjectID_x == subjectID_y' \
	-o $output/beta_source_target_baseline_samesubj.tsv
