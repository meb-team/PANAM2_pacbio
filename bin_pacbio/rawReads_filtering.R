#!/usr/bin/env Rscript --vanilla

# ---------------------------------------------------------------------------------------------- #

# ==================== #
# ===== Settings ===== #
# ==================== #

### Packages
pkg <- c("dada2", "Biostrings", "ggplot2", "reshape2", "ShortRead")
# Install packages not yet installed
installed_pkg <- pkg %in% rownames(installed.packages())
if (any(installed_pkg == FALSE)) {
  install.packages(pkg[!installed_pkg])
}
# Packages loading
invisible(lapply(pkg, library, character.only = TRUE))

### Input argument - path to working directory
args = commandArgs(trailingOnly=TRUE)
# setwd(args[1])

### Set path and construct path to fastq files
path <- args[1]
print(paste("path: ", path))

fastqfiles <- list.files(path, pattern=".fastq")

pfastqfiles <- file.path(path, basename(fastqfiles))

print(pfastqfiles)

sample <- tools::file_path_sans_ext(basename(fastqfiles))

print(sample)

### Set ggplot theme
theme_set(theme_bw())

# ---------------------------------------------------------------------------------------------- #

# ============================== #
# ===== Filter fastq files ===== #
# ============================== #

### Plot quality profile of the fastq files
# pdf("QualityProfile.pdf")
# plotQualityProfile(fn)
# dev.off()

### Plot length distribution of read sequences
# pdf("length.pdf")
# hist(nchar(getSequences(fn)), 100)
# dev.off()

# print("This is ok until here : 1")


### Filter fastq sequences
filteredfastq <- file.path(path, "filtered", basename(fastqfiles))
track <- filterAndTrim(pfastqfiles, filteredfastq, minQ=3, minLen=3000, maxLen=6000, maxN=0, rm.phix=FALSE, maxEE=4, verbose=TRUE)
track

# print("This is ok until here : 2")


### Export filtered fastq seqences as fasta file
for (i in sample) {
  sample_filtered <- paste0(path, "/filtered/", i, ".filtered.fasta")
  filteredfastq <- paste0(path, "/filtered/", i, ".fastq")
  writeFasta(readFastq(filteredfastq), sample_filtered, width=20000L)
}

