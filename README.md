# Diet_Intervention_Sperm methodology steps from scratch
This repository comprises of the steps and codes needed to reproduce the results of the Diet intervention paper

## 1. Download the fastq files 
Download the raw fastq files from GSE159752 Use SRA RUN Selector (https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA670419&o=acc_s%3Aa) and the SRA toolkit.

## 2. Run FASTQC
To check the quality of the reads in each sample use: _auto-run-fastqc.pl list-fastq-files.txt_. The perl script runs FASTQC on eash fastq file in the list.

## 3. Run Trimmomatic 
To remove adapter sequences and low quality reads use: _auto-run-trim-orig.pl list-fastq-files.txt_
The _auto-run-trim-orig.pl_ program runs Trimmomatic with the following parameters
> /Trimmomatic-0.39/trimmomatic-0.39.jar SE $file $name-trim-orig.fastq ILLUMINACLIP:~//Downloads/Trimmomatic-0.39/adapters/TruSeq3-SE-mod.fa:2:30:10 LEADING:3 TRAILING:3 CROP:46 SLIDINGWINDOW:4:15 MINLEN:10>

## 4. Rerun fastqc 
To ensure adapter removal use: _auto-run-fastqc.pl list-trimmed-fastq-files.txt_
