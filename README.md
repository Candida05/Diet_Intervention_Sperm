# Diet_Intervention_Sperm methodology steps from scratch
This repository comprises of the steps and codes needed to reproduce the results of the Diet intervention paper

## 1. Download the fastq files 
Download the raw fastq files from GSE159752 Use SRA RUN Selector (https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA670419&o=acc_s%3Aa) and the SRA toolkit.

To download SRA toolkit:
- wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.2/sratoolkit.3.0.2-ubuntu64.tar.gz
- tar -zxvf sratoolkit.3.0.2-ubuntu64.tar.gz          # extract file 
- export PATH=$PATH:~/sratoolkit.3.0.2-ubuntu64/bin   # add binaries to path using export path or editing ~/.bashrc file
- which fastq-dump                                    # verify the binaries added to the system path

To download fastq files:
- prefetch SRR12858025
- fastq-dump SRR12858025

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

 - Input files: exceRpt-Pipeline-2022-9-22-18-50-25-23-Sept-J1_exceRpt_readMappingSummary.txt, exceRpt-Pipeline-2022-9-23-16-57-54-24-Sept-J2_exceRpt_readMappingSummary.txt
 - Output files: J1_exceRpt_readMappingSummary-percent.txt (Pre-intervention) and J2_exceRpt_readMappingSummary-percent.txt (Post-intervention)

Using the following:
failed_quality_filter, failed_homopolymer_filter, UniVec_contaminants, rRNA, genome, not_mapped_to_genome_or_libs
miRNA_sense, miRNA_antisense, tRNA_sense, tRNA_antisense, piRNA_sense, piRNA_antisense, gencode_sense, gencode_antisense

Not using the following as expression counts are minimal:
All sRNA antisense, circularRNA_sense, circularRNA_antisense, repetitiveElements, endogenous_gapped
input_to_exogenous_miRNA, exogenous_miRNA, input_to_exogenous_rRNA, exogenous_rRNA, input_to_exogenous_genomes, exogenous_genomes

## 6. Run MiRDeep2 pipeline
For quantification and expression profiling of miRNA, the reads were separately mapped to human miRNA sequences using miRDeep2 pipeline. Mapper and quantifier modules of miRDeep2 pipeline were used for this purpose as follows:
 
>mapper.pl SRR12858025-trim-orig.fastq -e -g S25 -h -i -j -l 15 -m -s SRR12858025-trim-orig-col.fa -v

SRR12858025-trim-orig.fastq is the input fastq file and SRR12858025-trim-orig-col.fa is the output file with the collapsed reads

-e	input file is fastq format
-g	three-letter prefix for reads (by default 'seq')
-h	parse to fasta format
-i	convert rna to dna alphabet (to map against genome)
-j	remove all entries that have a sequence that contains letters other than a,c,g,t,u,n,A,C,G,T,U,N
-l	int discard reads shorter than int nts, default = 18
-m	collapse reads
-s	print processed reads to this file
-v	outputs progress report


>quantifier.pl -p hsa-release-22.1-prec-miRNA.fa -m hsa-release-22.1-mat-miRNA.fa -r SRR12858025-trim-orig-col.fa -t hsa -y SRR12858025-160922

files for mapping: hsa-release-22.1-prec-miRNA.fa (miRNA precursor file) and hsa-release-22.1-mat-miRNA.fa (mature miRNA file)

input file for mapping: SRR12858025-trim-orig-col.fa (collapsed file output of mapper module)

-p	precursor.fa miRNA precursor sequences from miRbase
-m	mature.fa miRNA sequences from miRbase
-r	reads.fa read sequences
-t	species

## 7. Run tRF & piRNA quantification pipeline
For the tRF quantification, the in-house developed pipeline by Kumar et al 2016 [PMID: 27263052] was used.
To run this pipeline permission is needed from the authors. 
For the piRNA quantification pipeline, Genboree's piRNA mapped data was used.   

## 8. Baseline Correlation Analysis
The baseline (pre-intervention) high confidence (expressed with a CPM >= 1 in all 17 pre-intervention samples) sncRNA (miRNA, piRNA tRF) are correlated with:
(i) Age (ii) BMI (iii) Sperm Concentration (iv) Sperm Motility 
To determine the effect of these factors on the sncRNA as well as to determine the confounding factors to be adjusted for in the differential expression analysis.

- Code: Script_Figure2_miRNA.R, Script_Figure2_tRFs.R, Script_Figure2_piRNA.R
  
- Input files: The TMM normalized expression data in CPM for the baseline sncRNA for pre-intervention (MV1) samples. The sample factor file comprising of age/BMI/sperm concentration and sperm motility.
  - miRNA: MV1-sample-factors-ed.txt, MV1-miRNA-CPM-TMM.txt
  - tRF: MV1-sample-factors-ed.txt, MV1-tRF-CPM-TMM.txt
  - piRNA: MV1-sample-factors-ed.txt, MV1-piRNA-CPM-TMM.txt

- Output files: correlation test output files comprising of the correlation test estimates and p-values. The figures comprise of the top two highly correlated sncRNA with age, BMI, sperm concentration and sperm motility.
  - miRNA: mirs-age-corr1811.txt, mirs-bmi-corr1811.txt, mirs-spconc-corr1811.txt, mirs-spmot-corr1811.txt, Figure2_miRNA.png
  - tRF: trfs-age-corr1811.txt, trfs-bmi-corr1811.txt, trfs-spconc-corr1811.txt, trfs-spmot-corr1811.txt, Figure2_piRNA.png
  - piRNA: pirs-age-corr1811.txt, pirs-bmi-corr1811.txt, pirs-spconc-corr1811.txt, pirs-spmot-corr1811.txt, Figure2_tRFs.png

## 9. Differential Expression Analysis
To obtain the sncRNA (miRNA, piRNA and tRF) getting altered due to the diet intervention, differential expression analysis is carried out using Limma double voom. In addition, blocking for ID was used to account for paired analysis. 
The sncRNA raw expression data was first filtered (CPM >=1 in any 8 samples) and normalised (TMM normalisation). Since age, BMI, sperm concentration and sperm motility were detected as confounding factors by the baseline correlation analysis, these were adjusted for in the design model. 
The following contrasts were carried out to determine the effect of the intervention over the control:

interventionEffect = (intervention.MV4 - intervention.MV1),
controlEffect = (control.MV4 - control.MV1), 
intervention_over_control = (intervention.MV4 - intervention.MV1) - (control.MV4 - control.MV1)
where MV4 is the time point 2 (T2) after the 6 weeks of intervention and MV1 is the time point 1 (T1) of the first visit

- Code: adjusted-miRDeep2-miRNA-ebayes.R, adjusted-PK-tRF-ebayes.R, adjusted-Genboree-piRNA-ebayes.R
  
- Input files: The raw and non-normalized count based expression data for all the samples. The sample factor file comprising of age/BMI/sperm concentration and sperm motility for all samples.
  - miRNA: All-sample-factors-ed.txt, miRDeep2-miRNA-expression-profile-orig.txt
  - tRF: All-sample-factors-ed.txt, PK-tRF-expression-profile-orig.txt
  - piRNA: All-sample-factors-ed.txt, Genboree-piRNA-expression-profile-orig.txt
 
 - Output files: Differentially expressed sncRNA with the logFC and p-values. These files need to be sorted in descending order of p-value to obtain the significant ones. The volcano plots depict the most significant differentially expressed sncRNA (p-val <0.01).   
   - miRNA: DE-miRDeep2-miRNAs.txt, volcano-DE-miRDeep2-miRNAs.tiff
   - tRF: DE-PK-trfs.txt, volcano-DE-PK-trfs.tiff
   - piRNA: DE-Genboree-piRNAs.txt, volcano-DE-Genboree-piRNAs.tiff

## 10. Genomic annotation enrichment of the sperm baseline sncRNAs Figure 1F 
To ascertian the genomic regions of all the mapped locations of the baseline sncRNA, the Bioconductor package regioneR was used. The package regioneR offers a statistical framework based on customizable permutation tests to assess the association between genomic region sets and other genomic features.

- Code: running-Region-miRNA.R, running-Region-tRF.R, running-Region-piRNA.R
- Input files: The bed file format of the genomic coordinates of the all the mapped locations of baseline sncRNA on the genome (some sncRNA can originate from multiple locations on the genome).
  - BL-mirs-all-matches-hg38-coord-6c-2608.bed, Bl-pirs-all-matches-hg18-6c-lifted-orig-hg38-2408-rem-random.bed, Bl-tRFs-all-matches-hg19-6c-lifted-orig-hg38.bed
  - The bed files of the following genomic region sets.
   chr1-22.gene.bed, chr1-22.intron.bed, chr1-22.codingexon.bed, chr1-22.exonplus.bed, chr1-22.CpG.bed, 
   chr1-22.5utr.bed, chr1-22.3utr.bed, chr1-22.Gene-2000Down.bed, chr1-22.Gene-2000UP.bed

 - Output results: z-scores and p-values for the enrichment in the different regions. 

## 11. Effects of dietary intervention on vitamin D and omega-3 fatty acid levels in circulation Figure 3A,B,C
To study the effects of dietary intervention, we measured the levels of vitamin D and omega-3 fatty acids in blood samples from all male participants (N=102) collected before and after 6 weeks of intervention. To also assess the sampling bias, these blood measures for the 17 subjects included in sncRNA analysis were compared with the full dataset.

- Code: Script_Figure3abc.R
- Input files: Data_SpermDietIntervention_102s.csv, Data_SpermDietIntervention_17s.csv
- Output files: Figure3A_change_in_conc_vitD_in_serum_nmol_L.png, Figure3B_change_in_percent_EPA_in_RBC.png, Figure3C_change_in_percent_DHA_in_RBC.png

