#!/bin/bash
# run_MultiCodes.sh metadata.csv /outdir /rawdatadir
for i in $(awk -F',' 'FNR>1 { print $4 }' $1)
do
	ip=$(printf '%02d' $i)
	iw=$(echo $ip | sed 's/^0*//')
	x=`find $3 -name "*_${iw}_S*.fastq"  -type f`
	echo $iw
	echo $x
	perl feba_pipeline/bin/MultiCodes.pl -out $2/data/MultiCode_tables/IT0${ip} -index IT0${ip} -bs3 < $x
done