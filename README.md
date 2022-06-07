## Genomic_data_Science_Class

This is DNA methylation tutorial for Genomic Data Science class in 1st semester, 2022 (Graduate class, Soongsil University)

This tutorial consists of preprocessing step, alignment step, differentially methylation analysis step and downstrem analysis step.

It is based on the reference paper. ("Multi-tissue DNA methylation age predictor in mouse", Thomas M. Stubbs et al, Genome Biology, 2017)

#### We used RRBS methylation dataset for DNA methylation tutorial.
This dataset is male C57BL/6-BABR mice.
Tissues (Liver, Lung, Heart and Cortex) were isolated from mice at Newborn, 14weeks, 27weeks and 41weeks.
our dataset is mice cortex tissues.
* SRR5195625
* SRR5195641
* SRR5195657
* SRR5195667

#### We uploaded these files:
* Coverage : coverage files after preprocessing step
* Differentially_methylation_analysis.R : Differentially methylation analysis code
* Preprocessing_code.pdf : Preprocessing code infomation
* readBismarkCoverage.R : The function code which it reads bismark coverage file.
* cpgi.mm10.bed.txt : The mm10 CpG island annotation file
* mm10_gencode25_enschrconversion.bed.txt : The mm10 transcript annotation file
* KEGG_pathway_DMC.csv : KEGG pathway analysis result from DAVID (Top 20 DMC (pval_ratio > 1))

#### These are the analysis tools we used.
* R >= 4.1.1
* FastQC 0.11.9
* samtools 1.15.1
* trim_galore 0.6.7
* UmiBam 0.2.0 (https://github.com/FelixKrueger/Umi-Grinder)
