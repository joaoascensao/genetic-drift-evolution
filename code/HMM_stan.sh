#!/bin/bash

for i in $(awk -F',' 'FNR>1' $1)
do
	python HMM_stan.py $4/data/bc_counts $4/data/drift $2 $i $3
done