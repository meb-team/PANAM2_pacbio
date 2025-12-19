#!/bin/bash

# USAGE
# sh rawReads_filtering.sh /path/to/working/dir /path/to/dir/containing/bamFiles

# ---------------------------------------------------------------------------------------------- #

# ==================== #
# ===== Settings ===== #
# ==================== #

### DEPENDENCIES
# samtools version 1.15.1
# R version 4.1.2
# This script requires an additional R script (rawReads_filtering.R)

### Setting variables
# working directory
WKDIR=$1
# directory containing bam files
BAMDIR=$2

SCRIPTS=$3 # change me to the correct path containing all script needed
THREADS=12

### Output log file
exec > $WKDIR/rawReads_filtering.log 2>&1

# ---------------------------------------------------------------------------------------------- #

# =============================== #
# ===== Raw Reads filtering ===== #
# =============================== #

cd $WKDIR

### -- Convert bam to fastq files -------------------------------------------------------------
echo "1. Converting bam files to fastq files"
echo

# -- Create a directory for fastq files storage
mkdir -p $WKDIR/1.fastq

# -- and convert
for bam in $(ls $BAMDIR | grep 'demux' | grep bam$ | sed 's/\.bam//g')
do
    samtools bam2fq $BAMDIR/$bam.bam > $WKDIR/1.fastq/$bam.fastq
done

### -- Filtering raw reads with DADA2 ---------------------------------------------------------
echo "2. Filtering fastq files and converting to fasta files"
echo

# -- Filter
Rscript --vanilla $SCRIPTS/rawReads_filtering.R $WKDIR/1.fastq

# -- Create a directory for filtered fasta files storage
mkdir -p $WKDIR/2.fasta
mv $WKDIR/1.fastq/filtered/*.filtered.fasta $WKDIR/2.fasta
cd $WKDIR/2.fasta
# -- Rename fasta files (!!! dash (-) vs mothur...)
rename 's/\-\-/_/' *.filtered.fasta
rename 's/filtered.//' *.fasta
# -- Move fasta files to 2.fasta directory
for i in $(ls | grep fasta | sed 's/\.fasta//g')
do
    mkdir $i
    mv $i.fasta $i
done 

echo "DONE !"
echo
