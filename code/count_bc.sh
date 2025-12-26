#!/bin/bash
# count_bc.sh experiment_list.csv metadata.csv /outdir

for i in $(awk -F',' 'FNR>1' $1)
do
	python count_bc.py $3/data/MultiCode_tables $3/data/errorcorrectedcodes $3/data/bc_counts $2 $i
done