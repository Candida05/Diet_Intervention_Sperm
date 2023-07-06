# Diet_Intervention_Sperm methodology steps from scratch
This repository comprises of the steps and codes needed to reproduce the results of the Diet intervention paper

## 1. Download the fastq files 
Download the raw fastq files from GSE159752 Use SRA RUN Selector (https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA670419&o=acc_s%3Aa) and the SRA toolkit.

## 2. Run FASTQC
To check the quality of the reads in each sample use: _auto-run-fastqc.pl list-fastq-files.txt_. The perl script runs FASTQC on eash fastq file in the list.

## 3. Run Trimmomatic 
To remove adapter sequences and low quality reads use: _auto-run-trim-orig.pl list-fastq-files.txt_
The _auto-run-trim-orig.pl_ program runs Trimmomatic with the following parameters
> /Trimmomatic-0.39/trimmomatic-0.39.jar SE $file $name-trim-orig.fastq ILLUMINACLIP:~//Downloads/Trimmomatic-0.39/adapters/TruSeq3-SE-mod.fa:2:30:10 LEADING:3 TRAILING:3 CROP:46 SLIDINGWINDOW:4:15 MINLEN:10

## 4. Rerun fastqc 
To ensure adapter removal use: _auto-run-fastqc.pl list-trimmed-fastq-files.txt_

## 5. Run Genboree
Using the trimmed fastq files run: 
### Genboree 4th Generation pipeline (v.6.2) ###
_Genome: hg38_
_Minimum read length: 18nt default_
_order of matching -> miRNA > tRNA> piRNA > Gencode > circRNA_

Run the _get-percentage.pl_ code on the pre and post intervention _exceRpt_readMappingSummary.txt_ files

input files: exceRpt-Pipeline-2022-9-22-18-50-25-23-Sept-J1_exceRpt_readMappingSummary.txt, exceRpt-Pipeline-2022-9-23-16-57-54-24-Sept-J2_exceRpt_readMappingSummary.txt
output files: J1_exceRpt_readMappingSummary-percent.txt (Pre-intervention) and J2_exceRpt_readMappingSummary-percent.txt

Using the following:
failed_quality_filter, failed_homopolymer_filter, UniVec_contaminants, rRNA, genome, not_mapped_to_genome_or_libs
miRNA_sense, miRNA_antisense, tRNA_sense, tRNA_antisense, piRNA_sense, piRNA_antisense, gencode_sense, gencode_antisense

Not using the following as expression counts are minimal:
All sRNA antisense, circularRNA_sense, circularRNA_antisense, repetitiveElements, endogenous_gapped
input_to_exogenous_miRNA, exogenous_miRNA, input_to_exogenous_rRNA, exogenous_rRNA, input_to_exogenous_genomes, exogenous_genomes

