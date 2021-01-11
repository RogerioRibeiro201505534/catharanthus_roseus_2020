#### SCRIPT TO RUN DESeq2 AND STORE RESULTS AND QC PLOTS ####

#misc
rm(list = ls())
setwd("C:/Users/Faculdade/Desktop/Dissertação/results 2020/mapping_to_genome/07_DESeq2/featureCounts")

#load Libraries
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library('DESeq2')
library('ggplot2')
library('RColorBrewer')
library('dplyr')
library('pheatmap')
library('ggrepel')
library('EnhancedVolcano')


#Only load one of the following set of variables !!!!!!!

#### Idioblasto EXPERIMENT ####
####    idio vs meso       ####
condition <- "meso_vs_idio"
sampleNames <- c('meso_1', 'meso_2', 'meso_3','idio_1', 'idio_2', 'idio_3')
sampleCondition <- c("meso","meso","meso","idio","idio","idio")
levels <- c("meso","idio")
sampleFiles <- c("meso_1_genome_counts.txt","meso_2_genome_counts.txt","meso_3_genome_counts.txt", "idio_1_genome_counts.txt", "idio_2_genome_counts.txt", "idio_3_genome_counts.txt")

#### LEAF IN/OUT EXPERIMENT ####
####       F1_in/out        ####
condition <- "F1in_vs_F1out"
sampleNames <- c('P1F1in', 'P2F1in', 'P3F1in','P1F1out', 'P2F1out', 'P3F1out')
sampleCondition <- c("f1in_f1in","f1in_f1in","f1in_f1in","f1out_f1out","f1out_f1out","f1out_f1out")
levels <- c("f1in_f1in","f1out_f1out")
sampleFiles <- c('F1P1_INT_genome_counts.txt', 'F1P2_INT_genome_counts.txt', 'F1P3_INT_genome_counts.txt','F1P1_EXT_genome_counts.txt', 'F1P2_EXT_genome_counts.txt', 'F1P3_EXT_genome_counts.txt')


###        F4_in/out        ####
condition <- "F4in_vs_F4out"
sampleNames <- c('P1F4in', 'P2F4in', 'P3F4in','P1F4out', 'P2F4out', 'P3F4out')
sampleCondition <- c("f4in_f4in","f4in_f4in","f4in_f4in","f4out_f4out","f4out_f4out","f4out_f4out")
levels <- c("f4in_f4in","f4out_f4out")
sampleFiles <- c('F4P1_INT_genome_counts.txt', 'F4P2_INT_genome_counts.txt', 'F4P3_INT_genome_counts.txt','F4P1_EXT_genome_counts.txt', 'F4P2_EXT_genome_counts.txt', 'F4P3_EXT_genome_counts.txt')



## create directories
# condition directory where all outputs are stored
out.dir <- condition
if (!(file.exists(out.dir))) {
  dir.create(out.dir)
}

# to store all csv files
results.dir <- paste(out.dir,"/results", sep="")
if (!(file.exists(results.dir))) {
  dir.create(results.dir)
}

# to store all plots
plots.dir <- paste(out.dir,"/plots" ,sep="")
if (!(file.exists(plots.dir))) {
  dir.create(plots.dir)
}

#### run deseq2 first on all data and explore data and normalized reads ####
# create dataframe that becomes a deseq table.
directory <- "counts"
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = as.factor(sampleCondition))
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ condition)

# The levels in colData are important because they are used in the log calculations; 
# it makes sense to set untreated or control first 
# so that the direction of the logs fold changes does not confuse everyone (typically we do comparisons to the control)
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition, levels = levels)

#pre filterning 
#filtering the dataset and removing 0 counts will improve computing times downstream

nrow(ddsHTSeq) # 37320
keep <- rowSums(counts(ddsHTSeq)) > 0
ddsHTSeq <- ddsHTSeq[keep,]
nrow(ddsHTSeq)

# Create DESeq2Dataset object
dds <- DESeq(ddsHTSeq)

####Count Normalization####
sizeFact <- sizeFactors(dds)
write.table(sizeFact, file = paste(results.dir,"/sizeFactors.txt", sep=""), sep="\t", quote=F, col.names=NA)

# To retrieve the normalized counts matrix from dds
normalizedCounts <- counts(dds, normalized=TRUE) 

# We can save this normalized data matrix to file for later use
write.table(normalizedCounts, file = paste(results.dir,"/normalizedCountsTable.txt", sep=""), sep="\t", quote=F, col.names=NA) # acrescentei

# Plot dispersion estimates
png(paste(plots.dir, "/_DispersionEstimates.png", sep=""), width=12,height=8, units='in', res=300)
print(plotDispEsts(dds))
dev.off()


### Sample-level QC
# To improve the distances/clustering for the PCA and hierarchical clustering visualization methods,
# we need to moderate the variance across the mean by applying the rlog transformation to the normalized counts
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE) # betaPriorVar:  if missing is estimated from the data

## PCA
# Plot PCA
png(paste(plots.dir, "/_PCA.png", sep=""), width=12,height=8, units='in', res=300)
print(plotPCA(rld, intgroup="condition"))
dev.off()

## Hierarchical Clustering Heatmap
# pheatmap() function requires a matrix/dataframe of numeric values as input,
# and so the first thing we need to is retrieve that information from the rld object
# Extract the rlog matrix from the object
rld_mat <- assay(rld)

# Compute pairwise correlation values
rld_cor <- cor(rld_mat)
head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames

# Plot heatmap
png(paste(plots.dir, "/_HC_heatmap.png", sep=""), width=12,height=8, units='in', res=300)
print(pheatmap(rld_cor))
dev.off()


#### Differential expression analysis ####
res <- results(dds, alpha = 0.05)
res <- res[order(res$padj),]
summary(res)
capture.output(summary(res), file = paste(results.dir,"/summaryRes.txt", sep=""))


#Get MA plots (plots log2foldchanges vs mean normalized counts)
png(paste(plots.dir, "/_MAplot.png", sep=""), width=12,height=8, units='in', res=300)
plotMA(res, ylim = c(-4,4), alpha = 0.05, main = "normal log2foldChanges")
dev.off()

#plot a volcano plot 
png(paste(plots.dir, "/_Volcano_plot.png", sep=""), width=12,height=8, units='in', res=300)
EnhancedVolcano(res, 
                lab = rownames(res), 
                x = "log2FoldChange", y = "padj",
                pCutoff = 0.05,
                FCcutoff = 1,
                title = condition, 
                legendLabels=c('Not sig.','Log (base 2) FC','adjusted p-value',
                               'ajudsted p-value & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 3.0)

dev.off()

### Set thresholds
padj.cutoff <- 0.05
actual.fold.change <- 2 # set the actual fold change threshold
lfc.cutoff <- round(log2(actual.fold.change), 2) 

#### Significant hits ####
#Get all genes with the actual fold change > 2 and p.adjusted cutoff < 0.05 from the results table 
sig <- subset(data.frame(res), padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
sig <- sig[order(sig$padj),]
n.sig <- nrow(sig)

# UP genes
n.sig.UP <- nrow(subset(sig, log2FoldChange > 0)) # number of UP genes
names.sig.UP <- rownames(subset(sig, log2FoldChange > 0))
length(names.sig.UP) == n.sig.UP

# DOWN genes
n.sig.DOWN <- nrow(subset(sig, log2FoldChange < 0)) # number of DOWN genes
names.sig.DOWN <- rownames(subset(sig, log2FoldChange < 0))
length(names.sig.DOWN)  == n.sig.DOWN

#### non significant hits ####
nonsig <- subset(data.frame(res), (padj >= padj.cutoff | (padj < padj.cutoff & abs(log2FoldChange) <= lfc.cutoff) | is.na(padj) == TRUE))
n.nonsig <- nrow(nonsig)
names.nonsig <- rownames(nonsig)

total.genes <- n.sig + n.nonsig


cat(paste(
paste("padj cutoff =",padj.cutoff),
paste("actual fold change cutoff =",actual.fold.change),
paste("number of significantly differentially expressed genes = ", n.sig, sep = ""),
paste("OVER expressed genes = ", n.sig.UP, sep = ""),
paste("UNDER expressed genes = ", n.sig.DOWN, sep = ""),
paste("OVER + UNDER expressed genes = number of diff expressed genes => ", n.sig.UP + n.sig.DOWN == n.sig, sep = ""),
paste("non significantly differentially expressed genes = ", n.nonsig, sep = ""),
paste("All genes are accounted for =>  ", total.genes == length(res@rownames), sep = "")
,sep = "\n"), file = paste(results.dir, "/_summaryOfHits.txt", sep=""))

#### write result to csv file ####
# save full results to file
write.csv(as.data.frame(res),file = paste(results.dir, "/_Res.csv", sep=""), quote = F)
# save significant hits to file
write.csv(as.data.frame(sig),file = paste(results.dir, "/_Sig.csv", sep=""), quote = F)
# save significantUP Gene names hits to file
write.csv(as.data.frame(names.sig.UP), file = paste(results.dir,"/_sigUP.csv", sep=""), quote = F)
# save significantDOWN Gene names to file
write.csv(as.data.frame(names.sig.DOWN), file = paste(results.dir, "/_sigDOWN.csv", sep=""), quote = F)
# save non-significant Gene names to file
write.csv(as.data.frame(names.nonsig), file = paste(results.dir, "/_nonsig.csv", sep=""), quote = F)


