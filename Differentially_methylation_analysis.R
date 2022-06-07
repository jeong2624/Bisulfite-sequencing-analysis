### Load methylKit package.
library(methylKit)
library(gplots)
library(plotrix)
library(genomation)

##### Reading the methylation coverage file

### Assign methylation coverage files to file.list variable.
# Note) You should write the absolute path of methylation coverage files.

samples = read.table("qFile_Cortex.txt", header = TRUE)
x = 1:length(samples[,1])
file.list = as.list(sapply(x, function(x) toString(samples[x, 1])))
sample.list = as.list(sapply(x, function(x) toString(samples[x, 2])))

### Read the files to a methylRawList object: myobj

# Note)
# To read bismark coverage files, you should download "readBismarkFiles.R" Rscript file from github.
# github homepage directory : https://gist.github.com/al2na/4839e615e2401d73fe51


# Note)
# treatment vector consists of 0 or 1. ("0" is defined as "Control", otherwise, "1" is defined as "Test")
source("readBismarkCoverage.R")

myobj = readBismarkCoverage(file.list,
                            sample.id = sample.list,
                            assembly = "mm10",
                            treatment = c(rep(1, 16)),
                            context = "CpG", min.cov = 5)

### Merge samples.

# Note)
# destrand : if TRUE, reads covering both strands of a CpG dinucleotide will be merged, do not set to TRUE if not only interested in CpGs (default: FALSE). 
# If the methylRawList object contains regions rather than bases setting destrand to TRUE will have no effect.

methunfil = unite(myobj, min.per.group = 16L)

### Filtering samples based on read coverage.
methunfil = filterByCoverage(myobj, lo.count = 5, lo.perc = NULL,
                             hi.count = NULL, hi.perc = 99.9)

# Get sample ids.
getSampleID(methunfil)

#### test vs control condition setting

### A object - 27week vs. 1week

### B object - 41week vs. 1week

meth_A = reorganize(methunfil, sample.ids = c("SRR5195653_cortex_1", "SRR5195657_cortex_1", "SRR5195660_cortex_1",
                                              "SRR5195671_cortex_1", "SRR5195637_cortex_27", "SRR5195641_cortex_27",
                                              "SRR5195645_cortex_27", "SRR5195649_cortex_27", "SRR5195675_cortex_27"),
                    treatment=c(rep(0, 4), rep(1, 5)))

meth_B = reorganize(methunfil, sample.ids = c("SRR5195653_cortex_1", "SRR5195657_cortex_1", "SRR5195660_cortex_1",
                                              "SRR5195671_cortex_1", "SRR5195621_cortex_41", "SRR5195625_cortex_41",
                                              "SRR5195629_cortex_41", "SRR5195633_cortex_41"), 
                    treatment=c(rep(0, 4), rep(1, 4)))

### Merge samples and get sample id.
meth_A = unite(meth_A)
getSampleID(meth_A)

### Merge samples and get sample id.
meth_B = unite(meth_B)
getSampleID(meth_B)

### Find differential Methylation site. (27week vs. 1week)
myDiff_A = calculateDiffMeth(meth_A)
head(myDiff_A)

# Count CpG sites if FDR < 1%. (27week vs. 1week)
myDiff25p_A = getMethylDiff(myDiff_A, difference = 25, qvalue = 0.01)
head(myDiff25p_A)

dim(myDiff25p_A)

# Select Top 20 CpG sites by methylation difference (27week vs. 1week)
head(myDiff25p_A[order(myDiff25p_A$pvalue), ], 20)

### Find differential Methylation site. (41week vs. 1week)
myDiff_B = calculateDiffMeth(meth_B)
head(myDiff_B)

# Count CpG sites if FDR < 1%. (41week vs. 1week)
myDiff25p_B = getMethylDiff(myDiff_B, difference = 25, qvalue = 0.01)
head(myDiff25p_B)

dim(myDiff25p_B)

# Select Top 20 CpG sites if FDR < 5% (41week vs. 1week)
head(myDiff25p_B[order(myDiff25p_B$pvalue), ], 20)

A_fdr_5 = subset(getData(myDiff_A), subset = qvalue < 0.01)
head(A_fdr_5[order(A_fdr_5$pvalue), ], 20)

B_fdr_5 = subset(getData(myDiff_B), subset = qvalue < 0.01)
head(B_fdr_5[order(B_fdr_5$pvalue), ], 20)

### Annotate differentially methylated bases/regions###

# read-in transcript locations to be used in annotation
# IMPORTANT: annotation files that come with the package (version >=0.5) are a subset of full annotation files.
# Download appropriate annotation files from UCSC (or other sources) in BED format

gene.obj = readTranscriptFeatures("mm10_gencode25_enschrconversion.bed.txt")

# annotate differentially methylated Cs with promoter/ exon / intron using annotation data (27week vs. 1week)
annotateWithGeneParts(as(myDiff25p_A,"GRanges"), gene.obj, intersect.chr = TRUE)

# plot the percentage of differentially methylated bases overlapping with exon/intron/promoters (27week vs. 1week)
diffAnn_A = annotateWithGeneParts(as(myDiff25p_A,"GRanges"), gene.obj, intersect.chr = TRUE)

plotTargetAnnotation(diffAnn_A, precedence = TRUE,
                     main="differential methylation annotation (27week vs. 1week)")

# annotate differentially methylated Cs with promoter/ exon / intron using annotation data (41week vs. 1week)
annotateWithGeneParts(as(myDiff25p_B,"GRanges"), gene.obj, intersect.chr = TRUE)

# plot the percentage of differentially methylated bases overlapping with exon/intron/promoters (41week vs. 1week)
diffAnn_B = annotateWithGeneParts(as(myDiff25p_B,"GRanges"), gene.obj, intersect.chr = TRUE)

plotTargetAnnotation(diffAnn_B, precedence = TRUE,
                     main="differential methylation annotation (41week vs. 1week)")

# get the distance to TSS and nearest gene name. (27week vs. 1week)
diffAnn_A = annotateWithGeneParts(as(myDiff25p_A, "GRanges"), gene.obj, intersect.chr = TRUE)
head(getAssociationWithTSS(diffAnn_A))

TSS_A = getAssociationWithTSS(diffAnn_A)
head(TSS_A[order(abs(TSS_A$dist.to.feature)), ], 10)

# get the distance to TSS and nearest gene name. (41week vs. 1week)
diffAnn_B = annotateWithGeneParts(as(myDiff25p_B, "GRanges"), gene.obj, intersect.chr = TRUE)
head(getAssociationWithTSS(diffAnn_B))

TSS_B = getAssociationWithTSS(diffAnn_B)
head(TSS_B[order(abs(TSS_B$dist.to.feature)), ], 10)

# read the shores and flanking regions and name CpG islands as CpGi. (27week vs. 1week)
cpg.obj_A = readFeatureFlank("cpgi.mm10.bed.txt", feature.flank.name = c("CpGi", "shores"))

# convert methylDiff object to GRanges and annotate
diffCpGann_A = annotateWithFeatureFlank(as(myDiff25p_A,"GRanges"),
                                        cpg.obj_A$CpGi, cpg.obj_A$shores,
                                        feature.name = "CpGi", flank.name = "shores")

diffCpGann_A

plotTargetAnnotation(diffCpGann_A, col=c("green","gray","white"),
                     main="differential methylation annotation (27week vs. 1week)")

# read the shores and flanking regions and name CpG islands as CpGi. (41week vs. 1week)
cpg.obj_B = readFeatureFlank("cpgi.mm10.bed.txt", feature.flank.name = c("CpGi", "shores"))

# convert methylDiff object to GRanges and annotate
diffCpGann_B = annotateWithFeatureFlank(as(myDiff25p_B,"GRanges"),
                                        cpg.obj_B$CpGi, cpg.obj_B$shores,
                                        feature.name = "CpGi", flank.name = "shores")

diffCpGann_B

plotTargetAnnotation(diffCpGann_B, col=c("green","gray","white"),
                     main="differential methylation annotation (41week vs. 1week)")

# Select DMC if q-value < 0.05 (27week vs. 1week)
myDiff_A = getMethylDiff(myDiff_A, difference = 25, qvalue = 0.05)

# Select DMC if q-value < 0.05 (41week vs. 1week)
myDiff_B = getMethylDiff(myDiff_B, difference = 25, qvalue = 0.05)

# merge these data
all_diff = merge(getData(myDiff_A), getData(myDiff_B), by = c("chr", "start", "end"), all = TRUE)

head(all_diff)

# get p-value ratio
all_diff$pval_ratio = log(all_diff$pvalue.x / all_diff$pvalue.y)

head(all_diff)

# select pval_ratio > 1
select_all_diff = subset(all_diff, subset = pval_ratio > 1)

dim(select_all_diff)

# Top 20 DMC (pval_ratio > 1)
head(select_all_diff[order(-select_all_diff$pval_ratio), ], 20)

# Get transcript ID related to DMC if p-value ratio > 1. 
all_diff = merge(getData(myDiff_A), getData(myDiff_B), by = c("chr", "start", "end"), all = TRUE)
all_diff$pval_ratio = log(all_diff$pvalue.x / all_diff$pvalue.y)

select_all_diff = subset(all_diff, subset = pval_ratio > 1)
pval_1 = annotateWithGeneParts(as(select_all_diff, "GRanges"), gene.obj, intersect.chr = TRUE)

TSS_pval_1 = getAssociationWithTSS(pval_1)
head(TSS_pval_1)

uniq_gene = unique(TSS_pval_1$feature.name)

