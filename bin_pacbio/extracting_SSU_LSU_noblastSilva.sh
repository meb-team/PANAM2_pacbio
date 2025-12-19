#!/bin/bash

# USAGE
# filter_to_otus.sh <fasta file containing filtered sequences with full path> <number of thread>

# USAGE with parallel
# nohup ls -d -1 /path/to/dir/containing/fasta_file/*.fasta | parallel --jobs 4 'sh extracting_SSU_LSU.sh {}' &


# ---------------------------------------------------------------------------------------------- #

# ==================== #
# ===== Settings ===== #
# ==================== #

### DEPENDENCIES
# vsearch version 2.3.4
# mafft version 7.310
# mothur version 1.39.5
# ncbi-blast+ version 2.9.0-2
# barrnap version 0.9

# This script requires SILVA V138.1 (or newer version) SSURef NR99 fasta file.
# Before using it, make a blast database with the SILVA sequences :
# makeblastdb -in SILVA_138.1_SSURef_NR99_tax_silva.fasta -title "SILVA 138_1r SSURef N99" -dbtype nucl

### Setting variables
VAR=$1
DIR=$(dirname "${VAR}")
BASE=$(basename "${VAR}") 
SAMPLE="${BASE%.*}"
SCRIPTS=$(pwd)  # change me if needed
THREADS=4       # change me if needed
DB="/datater/terzian/Euks-PacBio-2021_redo/database/SILVA/SILVA_138.1_SSURef_NR99_tax_silva.fasta"  # change me to the correct path if needed

### Output log file
exec > $DIR/filter_to_otus.log 2>&1

cd $DIR

# ---------------------------------------------------------------------------------------------- #

# ========================================== #
# ===== Extracting SSU and LSU regions ===== #
# ========================================== #

### -- Preclustering fasta filtered sequences -------------------------------------------------------------------
echo "1. Preclustering all sequences (99% id) with vsearch"
echo

mkdir preclusters
vsearch --cluster_fast $BASE --id 0.99 --threads "${THREADS}" --clusters preclusters/precluster.c_ --uc precluster.uc

### -- Separate preclusters with more than 2 sequences 
mkdir preclusters/large_preclusters
cd preclusters
for i in precluster.c_*
do
    #echo $i
    number=$(grep -c ">" $i)
    if (( $number > 2 ))
    then
        mv $i large_preclusters/
    fi
done
cd large_preclusters

### -- construct count table for mothur
for i in precluster.c_*
do
    x=$(echo $i | cut -d '.' -f2)
    count=$(grep -c ">" $i)
    echo -e ""$x"_conseq\t"$count"" >> $SAMPLE.count_table
done
mv $SAMPLE.count_table ./..

### -- Align large preclusters with mafft -----------------------------------------------------------------------
echo "2. Aligning large preclusters (>2 seq) with mafft..."
echo
for i in precluster.c_*
do
    mafft --quiet --thread "${THREADS}" $i > $i.aligned.fasta
done

### -- Generate majority rule consensus sequences
echo "...and generate majority rule consensus sequences with mothur"
for i in *.aligned.fasta
do
    mothur -q "#consensus.seqs(fasta=$i, cutoff=51)"
    #mothur "#consensus.seqs(fasta=$i, cutoff=51)"
done
rm mothur*  # rm mothur log file...

### -- Header name
for i in *.aligned.cons.fasta
do
    x=$(echo $i | cut -d '.' -f2)
    sed -E "s/(>)(.*)/\1${x}_\2/" $i > tmp.$i
    mv tmp.$i $i
done

mv *.aligned.cons.fasta ../
cd ../
cat precluster.c_* > $SAMPLE.preclusters.fasta

### -- Add small-preclusters to count table
grep ">m" $SAMPLE.preclusters.fasta | while read line; do header=$(echo $line | sed -E 's/>(.*)/\1/')
echo -e "$header\t1" >> $SAMPLE.count_table; done
echo -e "Representative_Sequence\ttotal" | cat - $SAMPLE.count_table > temp && mv temp $SAMPLE.count_table

mv $SAMPLE.count_table ../
mv $SAMPLE.preclusters.fasta ../
cd ../

### -- Gap Curation with mothur
mothur "#degap.seqs(fasta=$SAMPLE.preclusters.fasta)"

### -- Decontamination - BLAST against SILVA --------------------------------------------------------------------
# echo "3. Remove prokaroytic sequences - BLAST againt SILVA"


# blastn -query $SAMPLE.preclusters.ng.fasta -db $DB -evalue 1e-10 -num_threads "${THREADS}" -out $SAMPLE.preclusters.ng_vs_Silva.blastn -outfmt '6 std salltitles'

# cat $SAMPLE.preclusters.ng_vs_Silva.blastn | sort -k1,1 -k12,12nr |	awk '!seen[$1]++' > tophit.$SAMPLE.preclusters.ng_vs_Silva.blastn

# cat tophit.$SAMPLE.preclusters.ng_vs_Silva.blastn | grep -v "Eukaryota" | cut -f 1 > prok.list


# mothur -q "#remove.seqs(fasta=$SAMPLE.preclusters.ng.fasta, count=$SAMPLE.count_table, accnos=prok.list)"


### -- Detecting and extracting 18S and 28S wit barrnap ---------------------------------------------------------
echo "4. Extracting SSU and LSU wit barrnap"

barrnap --threads "${THREADS}" --reject 0.4 --kingdom euk -o barrnap.fasta < $SAMPLE.preclusters.ng.fasta > barrnap.gff

grep -A1 ">18S_rRNA" barrnap.fasta | sed 's/>.*::\(.*\):.*/>\1/g' > all18S.fasta
grep -A1 ">28S_rRNA" barrnap.fasta | sed 's/>.*::\(.*\):.*/>\1/g' > all28S.fasta

### lists of 18S / 28S / 18S+28S
for i in $(grep '>' $SAMPLE.preclusters.ng.fasta | sed 's/>//g')
do
    grep -w $i barrnap.gff |grep "18S" | cut -f 1 | sort | uniq -u
done > barrnap.18S.list

sort barrnap.18S.list > barrnap.18S.list.sorted

for i in $(grep '>' $SAMPLE.preclusters.ng.fasta | sed 's/>//g')
do
    grep -w $i barrnap.gff |grep "28S" | cut -f 1 | sort | uniq -u
done > barrnap.28S.list

sort barrnap.28S.list > barrnap.28S.list.sorted

comm -12 barrnap.18S.list.sorted barrnap.28S.list.sorted > barrnap.comm.list




