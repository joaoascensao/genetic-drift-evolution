#!/bin/bash

# make all necessary directories
sh mkdirs.sh ../E6/

# run MultiCodes.pl to parse barcodes from raw .fastq sequencing files
sh run_MultiCodes.sh ../E6/E6_meta.csv ../E6/ /media/joao/cff36999-03ac-4553-8b2d-a21de5bfec1f/tn7_final/E6/

# perform error correction to correct for sequencing errors in barcodes
sh error_correction.sh ../E6/E6_exps.csv ../E6/E6_meta.csv ../E6/

# count how many times we see each barcode for each sample
sh count_bc.sh ../E6/E6_exps.csv ../E6/E6_meta.csv ../E6/

# compute variance parameters for outlier detection
sh kt_errors.sh ../E6/E6_exps.csv ../E6/E6_meta.csv ../E6/E6_cfus.csv ../E6/

# detect outlier barcodes--those with trajectories that deviate significantly from neutrality
sh outliers.sh ../E6/E6_exps.csv ../E6/E6_meta.csv ../E6/

# filter out outlier barcodes, merge low-frequency barcodes together
sh filter_bcs.sh ../E6/E6_exps.csv ../E6/E6_meta.csv ../E6/E6_cfus.csv ../E6/

# MCMC inference of genetic drift parameters
sh HMM_stan.sh ../E6/E6_exps.csv ../E6/E6_meta.csv ../E6/E6_cfus.csv ../E6/

