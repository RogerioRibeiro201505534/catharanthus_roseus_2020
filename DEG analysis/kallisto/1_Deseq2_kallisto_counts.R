#### SCRIPT TO RUN DESEQ2 WITH KALLISTO COUNTS AND STORE RESULTS AND QC PLOTS ####
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library('DESeq2')
library('ggplot2')
library('RColorBrewer')
library('dplyr')
library('pheatmap')
library('ggrepel')
library('tximport')
library(EnhancedVolcano)

#Only load one of the following set of variables 

#### Idioblasto EXPERIMENT ####
####    idio vs meso       ####
condition <- "meso_vs_idio"
sampleNames <- c('meso_1', 'meso_2', 'meso_3','idio_1', 'idio_2', 'idio_3')
sampleCondition <- c("meso","meso","meso","idio","idio","idio")
levels <- c("meso","idio")
sampleFiles <- "abundance.tsv"


#### LEAF IN/OUT EXPERIMENT ####
####       F1_in/out        ####
condition <- "F1in_vs_F1out"
sampleNames <- c('F1P1_Int', 'F1P2_Int', 'F1P3_Int','F1P1_Ext', 'F1P2_Ext', 'F1P3_Ext')
sampleCondition <- c("f1in_f1in","f1in_f1in","f1in_f1in","f1out_f1out","f1out_f1out","f1out_f1out")
levels <- c("f1in_f1in","f1out_f1out")
sampleFiles <- "abundance.tsv"

###        F4_in/out        ####
condition <- "F4in_vs_F4out"
sampleNames <- c('F4P1_Int', 'F4P2_Int', 'F4P3_Int','F4P1_Ext', 'F4P2_Ext', 'F4P3_Ext')
sampleCondition <- c("f4in_f4in","f4in_f4in","f4in_f4in","f4out_f4out","f4out_f4out","f4out_f4out")
levels <- c("f4in_f4in","f4out_f4out")
sampleFiles <- "abundance.tsv"

## create directories
# condition directory where all outputs are stored
out.dir <- paste(condition, sep="")
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
#load the tx2geneTable to associate each transcript Id to gene id. Collum name is irelevant, but this must be 
#the order the datarame is loaded 

tx2gene <- read.csv(file = "transcript2gene.tsv", sep = "\t", header = F)

# create dataframe that becomes a deseq table.
directory <- "counts"
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = as.factor(sampleCondition))
files <- paste("counts", sampleTable$sampleName, sampleTable$fileName, sep = "/")
names(files) <- paste0(sampleTable$sampleName)
kallisto_Counts <- tximport(files, type="kallisto", dropInfReps = TRUE, tx2gene = tx2gene)

ddskallisto <- DESeqDataSetFromTximport(kallisto_Counts, colData = sampleTable, design = ~ condition)
#using counts and average transcript lengths from tximport

# The levels in colData are important because they are used in the log calculations; 
# it makes sense to set untreated or control first !!!!!!
# so that the direction of the logs fold changes does not confuse everyone (typically we do comparisons to the control)! 
colData(ddskallisto)$condition <- factor(colData(ddskallisto)$condition, levels = levels)


nrow(ddskallisto) #37320
keep <- rowSums(counts(ddskallisto)) > 0
ddskallisto <- ddskallisto[keep,]
nrow(ddskallisto)


# Create DESeq2Dataset object
dds <- DESeq(ddskallisto)

####Count Normalization####
#since in kallisto the size factor are computed by gene and by samples a matrix of gene x samples size factor are created 
#To get an idea of size factor the arithemic mean is calculated instead 

normalization_matrix <- normalizationFactors(dds)
estimated_size_fact <- colMeans(normalization_matrix)

write.table(estimated_size_fact, file = paste(results.dir,"/estsizeFactors.txt", sep=""), sep="\t", quote=F, col.names=NA)

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
                title = "meso vs idio", 
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
sig <- subset(data.frame(res), padj < padj.cutoff)
# Narrow the results using a log2 foldchange threshold as well
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
write.csv(as.data.frame(res),file = paste(results.dir, "/_Res.csv", sep=""))
# save significant hits to file
write.csv(as.data.frame(sig),file = paste(results.dir, "/_Sig.csv", sep=""))
# save significantUP Gene names hits to file
write.csv(as.data.frame(names.sig.UP), file = paste(results.dir,"/_sigUP.csv", sep=""), quote = F)
# save significantDOWN Gene names to file
write.csv(as.data.frame(names.sig.DOWN), file = paste(results.dir, "/_sigDOWN.csv", sep=""), quote = F)
# save non-significant Gene names to file
write.csv(as.data.frame(names.nonsig), file = paste(results.dir, "/_nonsig.csv", sep=""), quote = F)


