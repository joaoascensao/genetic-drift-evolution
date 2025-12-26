#!/bin/bash

for i in $(awk -F',' 'FNR>1' $1)
do
	python HMM_stan_combinereps.py $i
done