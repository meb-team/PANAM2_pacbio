PANAM2 is a pipeline for depicting the microbial community by implementing a phylogenetic approach for the annotation of the gene coding for 16S and 18S rRNA. It is based on PANAM (Taib et al. 2013) with some improvements and focusing on Illumina sequencing.
https://github.com/panammeb/PANAM2


**These script are part of the preprocessing of reads for the PANAM2 pipeline**

* Only the `Full_script_panam_pacbio.sh` script needs to be launch to perform preprocessing. 

* An example of command line :

```
nohup bash Full_script_panam_pacbio.sh -w <BIN_FOLDER> -x <WORKING_DIRECTORY> -y <RAW_READS_FOLDER> -z <AMPLICONS_file.txt> &> nohup.preprocessing_panam_pacbio.out &
```
	
* A brief description of `Full_script_panam_pacbio.sh` parameters :

	+ -w <BIN_FOLDER> : path to the scripts
	+ -x <WORKING_DIRECTORY> : path to the folder where all results and intermediate file will be stored
	+ -y <RAW_READS_FOLDER> : path to the raw fastq reads
	+ -z <AMPLICONS_file.txt> : path to the amplicons file (one "samplename\tsampleOrigin" per line)

Example of a amplicons content :

```
sample	sampleOrigin
Sample_1	demux.Location_xxx
Sample_2	demux.Location_xxx
Sample_3	demux.Location_xxx
```


* A conda env can be build using the Rpanam_env.yml (note: This environment is for preprocessing only at the moment)

	+ An example of conda command line to create the env : `conda env create -n Rpanam --file Rpanam_env.yml`
	+ Then, enter the environment : `conda activate Rpanam`