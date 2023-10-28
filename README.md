## Implementation of DNA methylation analysis in the paper below.

"Multi-tissue DNA methylation age predictor in mouse", Thomas M. Stubbs et al, Genome Biology, 2017)

#### Project time : 2022.03.09 - 2022.04.30

#### Description :
* The project was conducted as part of "geomic data science" class in the Department of Bioinformatics at Soongsil University, Korea.
* The purpose of this project is to learn how to analyze DNA methylation datasets.
* According to this project, we try to replicate result from the paper and know how the genome data science process works.
* It mainly deals with preprocessing step, alignment step and differentially methylation analysis

#### We used RRBS methylation dataset for DNA methylation tutorial.
This dataset is male C57BL/6-BABR mice.
Tissues (Liver, Lung, Heart and Cortex) were isolated from mice at Newborn, 14weeks, 27weeks and 41weeks.

#### our dataset is mice cortex tissues. (Download from [SRA Explorer](https://sra-explorer.info/))
* SRR5195663
* SRR5195667
* SRR5195679
* SRR5195653
* SRR5195657
* SRR5195660
* SRR5195671
* SRR5195637
* SRR5195641
* SRR5195645
* SRR5195649
* SRR5195675
* SRR5195621
* SRR5195625
* SRR5195629
* SRR5195633

#### We uploaded these files:
* Coverage : coverage files after preprocessing step
* Differentially_methylation_analysis.R : Differentially methylation analysis code
* Preprocessing_code.pdf : Preprocessing code information
* readBismarkCoverage.R : The function code which it reads bismark coverage file.
* cpgi.mm10.bed.txt : The mm10 CpG island annotation file
* mm10_gencode25_enschrconversion.bed.txt : The mm10 transcript annotation file
* KEGG_pathway_DMC.csv : KEGG pathway analysis result from DAVID (Top 20 DMC (pval_ratio > 1))

#### These are the analysis tools we used.
* R >= 4.1.1
* FastQC 0.11.9
* samtools 1.15.1
* trim_galore 0.6.7
* Bismark 0.23.1
* UmiBam 0.2.0 (https://github.com/FelixKrueger/Umi-Grinder)
