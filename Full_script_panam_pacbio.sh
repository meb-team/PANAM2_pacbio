#!/bin/bash

helpFunction()
{
   echo ""
   echo "Usage: $0 -w binpath -x dirpath -y bamdir -z amplicons"
   echo -e "\t-w Full path to bin folder containing scripts"
   echo -e "\t-x Full path to working directory"
   echo -e "\t-y Full path to bam folder"
   echo -e "\t-z Full path to AMPLICONS file"
   exit 1 # Exit script after printing help
}

while getopts "w:x:y:z:" opt
do
   case "$opt" in
      w ) binpath="$OPTARG" ;;
      x ) dirpath="$OPTARG" ;;
      y ) bamdir="$OPTARG" ;;
      z ) amplicons="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$binpath" ] || [ -z "$dirpath" ] || [ -z "$bamdir" ] || [ -z "$amplicons" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

#Begin script in case all parameters are correct
echo "$binpath"
echo "$dirpath"
echo "$bamdir"
echo "$amplicons"


mkdir -p $dirpath/exp/1.rawReads_Filtering
mkdir -p $dirpath/exp/2.SSU-LSU_Extraction/
mkdir -p $dirpath/exp/3.OTU_table/quality_output

bash $binpath/rawReads_filtering.sh $dirpath/exp/1.rawReads_Filtering $bamdir $binpath


### need to create 2.SSU-LSU_Extraction/ dir
cp -r $dirpath/exp/1.rawReads_Filtering/2.fasta/* $dirpath/exp/2.SSU-LSU_Extraction/

### Need to keep double quote so $binpath is interpreted as a variable and not a string as with single quote
nohup ls -d -1 $dirpath/exp/2.SSU-LSU_Extraction/demux*/demux*.fasta | parallel --jobs 4 "bash $binpath/extracting_SSU_LSU_noblastSilva.sh {}"

for i in $(ls $dirpath/exp/2.SSU-LSU_Extraction/ | grep demux); do mkdir $dirpath/exp/3.OTU_table/quality_output/$i; \
cp $dirpath/exp/2.SSU-LSU_Extraction/$i/all18S.fasta $dirpath/exp/3.OTU_table/quality_output/$i; \
cp $dirpath/exp/2.SSU-LSU_Extraction/$i/$i.count_table $dirpath/exp/3.OTU_table/quality_output/$i; done \


for i in $(ls $dirpath/exp/3.OTU_table/quality_output | grep demux);
do vsearch --derep_full $dirpath/exp/3.OTU_table/quality_output/$i/all18S.fasta --output $dirpath/exp/3.OTU_table/quality_output/$i/all18S.derep.fasta --sizeout;
vsearch --uchime_denovo $dirpath/exp/3.OTU_table/quality_output/$i/all18S.derep.fasta --chimeras $dirpath/exp/3.OTU_table/quality_output/$i/chimera.out.fasta --nonchimeras $dirpath/exp/3.OTU_table/quality_output/$i/non.chimera.out.fasta;
vsearch -search_exact $dirpath/exp/3.OTU_table/quality_output/$i/all18S.fasta -db $dirpath/exp/3.OTU_table/quality_output/$i/non.chimera.out.fasta -uc $dirpath/exp/3.OTU_table/quality_output/$i/current_mapping.uc;
grep '100.0' $dirpath/exp/3.OTU_table/quality_output/$i/current_mapping.uc | cut -f9 > $dirpath/exp/3.OTU_table/quality_output/$i/list_extract.txt;
seqtk subseq $dirpath/exp/3.OTU_table/quality_output/$i/all18S.fasta $dirpath/exp/3.OTU_table/quality_output/$i/list_extract.txt > $dirpath/exp/3.OTU_table/quality_output/$i/all18S.filtered.fasta; done

cd $dirpath/exp/3.OTU_table/quality_output/

for i in $(ls $dirpath/exp/3.OTU_table/quality_output/ | grep demux); do perl $binpath/seqAll_pacbio_chimeras.pl $i $amplicons ;
done

cat $dirpath/exp/3.OTU_table/quality_output/demux*/*.final.fasta > $dirpath/exp/3.OTU_table/quality_output/seqAll.fasta