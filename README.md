# Diet_Intervention_Sperm methodology steps from scratch
This repository comprises of the steps and codes needed to reproduce the results of the Diet intervention paper

## 1. Download the fastq files 
Download the raw fastq files from GSE159752 Use SRA RUN Selector (https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA670419&o=acc_s%3Aa) and the SRA toolkit.

## 2. Run FASTQC
To check the quality of the reads in each sample use: 
>auto-run-fastqc.pl list-fastq-files.txt.
The perl script runs FASTQC on eash fastq file in the list.

## 3. Run Trimmomatic 
To remove adapter sequences and low quality reads use: auto-run-trim-orig.pl list-fastq-files.txt
The auto-run-trim-orig.pl program runs Trimmomatic with the following parameters
> /Trimmomatic-0.39/trimmomatic-0.39.jar SE $file $name-trim-orig.fastq ILLUMINACLIP:~//Downloads/Trimmomatic-0.39/adapters/TruSeq3-SE-mod.fa:2:30:10 LEADING:3 TRAILING:3 CROP:46 SLIDINGWINDOW:4:15 MINLEN:10

## 4. Rerun FASTQC
To ensure adapter removal use: 
>auto-run-fastqc.pl list-trimmed-fastq-files.txt

## 5. Run Genboree
Using the trimmed fastq files run: 
### Genboree 4th Generation pipeline (v.6.2) ###
>Genome: hg38
>Minimum read length: 18nt default
>order of matching -> miRNA > tRNA> piRNA > Gencode > circRNA

>Run the get-percentage.pl code on the pre and post intervention exceRpt_readMappingSummary.txt files

Input files: exceRpt-Pipeline-2022-9-22-18-50-25-23-Sept-J1_exceRpt_readMappingSummary.txt, exceRpt-Pipeline-2022-9-23-16-57-54-24-Sept-J2_exceRpt_readMappingSummary.txt
Output files: J1_exceRpt_readMappingSummary-percent.txt (Pre-intervention) and J2_exceRpt_readMappingSummary-percent.txt

Using the following:
failed_quality_filter, failed_homopolymer_filter, UniVec_contaminants, rRNA, genome, not_mapped_to_genome_or_libs
miRNA_sense, miRNA_antisense, tRNA_sense, tRNA_antisense, piRNA_sense, piRNA_antisense, gencode_sense, gencode_antisense

Not using the following as expression counts are minimal:
All sRNA antisense, circularRNA_sense, circularRNA_antisense, repetitiveElements, endogenous_gapped
input_to_exogenous_miRNA, exogenous_miRNA, input_to_exogenous_rRNA, exogenous_rRNA, input_to_exogenous_genomes, exogenous_genomes

## 6. Run MiRDeep2 pipeline
## 7. Run Kumar et al tRF quantification pipeline

## 8. Baseline Correlation Analysis
The baseline (pre-intervention) high confidence (expressed with a CPM >= 1 in all 17 pre-intervention samples) sncRNA (miRNA, piRNA tRF) are correlated with:
(i) Age (ii) BMI (iii) Sperm Concentration (iv) Sperm Motility 
To determine the effect of these factors on the sncRNA as well as to determine the confounding factors to be adjusted for in the differential expression analysis.

Input files: The TMM normalized expression data in CPM for the baseline sncRNA for pre-intervention (MV1) samples. The sample factor file comprising of age/BMI/sperm concentration and sperm motility.
miRNA: MV1-sample-factors-ed.txt, MV1-miRNA-CPM-TMM.txt
tRF: MV1-sample-factors-ed.txt, MV1-tRF-CPM-TMM.txt
piRNA: MV1-sample-factors-ed.txt, MV1-piRNA-CPM-TMM.txt

Output files: correlation test output files comprising of the correlation test estimates and p-values. 
miRNA: mirs-age-corr1811.txt, mirs-bmi-corr1811.txt, mirs-spconc-corr1811.txt, mirs-spmot-corr1811.txt
tRF: trfs-age-corr1811.txt, trfs-bmi-corr1811.txt, trfs-spconc-corr1811.txt, trfs-spmot-corr1811.txt
piRNA: pirs-age-corr1811.txt, pirs-bmi-corr1811.txt, pirs-spconc-corr1811.txt, pirs-spmot-corr1811.txt

## 9. Differential Expression Analysis
To obtain the sncRNA (miRNA, piRNA and tRF) getting altered due to the diet intervention, differential expression analysis is carried out using Limma double voom. In addition, blocking for ID was used to account for paired analysis. 
The sncRNA raw expression data was first filtered (CPM >=1 in any 8 samples) and normalised (TMM normalisation). Since age, BMI, sperm concentration and sperm motility were detected as confounding factors by the baseline correlation analysis, these were adjusted for in the design model. 
The following contrasts were carried out to determine the effect of the intervention over the control:

interventionEffect = (intervention.MV4 - intervention.MV1),
controlEffect = (control.MV4 - control.MV1), 
intervention_over_control = (intervention.MV4 - intervention.MV1) - (control.MV4 - control.MV1)



