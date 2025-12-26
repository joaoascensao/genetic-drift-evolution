#!/bin/bash
# kt_errors.sh experiment_list.csv metadata.csv cfus.csv /outdir 

for i in $(awk -F',' 'FNR>1' $1)
do
	python kt_errors.py $4/data/bc_counts $4/data/kt_errors $2 $i $3
done