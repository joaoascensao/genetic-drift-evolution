#!/bin/bash
# single_error_correction.sh experiment_list.csv metadata.csv /outdir


for i in $(awk -F',' 'FNR>1' $1)
do
	python error_correction_1.py $3/data/MultiCode_tables $3/data/errorcorrectedcodes $2 $i
done


for i in $(awk -F',' 'FNR>1 { print $4 }' $2)
do
	echo "Starting levenshtein distance calculation for primer $i"
	./error_correction_levenshtein $3/data/errorcorrectedcodes/${i}_cc0.txt $3/data/errorcorrectedcodes/${i}_high_count_nodes.txt $3/data/errorcorrectedcodes/${i}_leven_out.txt
done

for i in $(awk -F',' 'FNR>1' $1)
do
	python error_correction_2.py $3/data/MultiCode_tables $3/data/errorcorrectedcodes $2 $i
done