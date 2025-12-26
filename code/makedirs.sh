#!/bin/bash
# sh makedirs.sh outdir
outdir=$1

mkdir -p $outdir
mkdir -p $outdir/data
mkdir -p $outdir/data/MultiCode_tables/
mkdir -p $outdir/data/errorcorrectedcodes/
mkdir -p $outdir/data/bc_counts/
mkdir -p $outdir/data/kt_errors/
mkdir -p $outdir/data/outliers/
mkdir -p $outdir/data/drift/
mkdir -p $outdir/data/drift/figures/
mkdir -p $outdir/data/bc_counts/figures/


echo 'Made directories in '$1