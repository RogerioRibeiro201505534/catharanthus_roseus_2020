#### SCRIPT TO RUN DESEQ2 AND STORE RESULTS AND QC PLOTS ####
rm(list = ls())
setwd("C:/Users/Faculdade/Desktop/Dissertação/results 2020/nanopore/06_DE")


#load library 
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library('DESeq2')
library('ggplot2')
library('RColorBrewer')
library('dplyr')
library('pheatmap')
library('ggrepel')
library("Rtsne")
library("tximport")


sampleNames <- c('P1F2Int', 'P2F2Int', 'P3F2Int','P1F2Ext', 'P2F2Ext', 'P3F2Ext', 'P1F3Int', 'P2F3Int', 'P3F3Int','P1F3Ext', 'P2F3Ext', 'P3F3Ext')
sampleCondition <- c("F2Int","F2Int","F2Int","F2ext","F2ext","F2ext", "f3Int","f3Int","f3Int","f3ext","f3ext","f3ext")
sampleFiles <- c("quant.sf")

tx2gene <- read.csv(file = "annotation/transcript2geneID.tsv", sep = "\t", header = F)


sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
sampleTable$condition <- as.factor(sampleTable$condition)
files <- paste("salmonQuant", sampleTable$sampleName, sampleTable$fileName, sep = "/")
names(files) <- paste0(sampleTable$sampleName)
salmon_quant <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut = F)

ddsalmon <- DESeqDataSetFromTximport(salmon_quant, colData = sampleTable, design = ~ condition)


nrow(ddsalmon) #37447
keep <- rowSums(counts(ddsalmon)) > 0
ddsalmon <- ddsalmon[keep,]
nrow(ddsalmon) #19397

dds <- DESeq(ddsalmon)


### Sample-level QC
# To improve the distances/clustering for the PCA and hierarchical clustering visualization methods,
# we need to moderate the variance across the mean by applying the rlog transformation to the normalized counts
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE) # betaPriorVar:  if missing is estimated from the data

homedir = "PCA_all_samples"
if (!dir.exists(homedir)){
  dir.create(homedir)
}


#Set a colour pallete 
colour_blind_pallete <- c("#ffcc33", "#0072B2", "#D55E00", "#CC79A7")
names(colour_blind_pallete) <-  c("F2in","F2out", "F3in","F3out")

# Plot PCA
pca_data <- plotPCA(rld, intgroup="condition", returnData = T)
#added to match the plots from featureCounts
pca_data$group <- c("F2in", "F2in", "F2in", "F2out", "F2out", "F2out", "F3in", "F3in", "F3in", "F3out", "F3out", "F3out")

p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +   
  geom_point(size = 5, aes(shape = group), position = "jitter") + 
  scale_color_manual(breaks = names(colour_blind_pallete), values=colour_blind_pallete) + 
  scale_shape_manual(values = c(0,1,2,3,4,6,7,8)) + 
  theme_classic() + 
  theme(legend.title=element_blank(), legend.text = element_text(size = 12), plot.title = element_text(size=22)) + 
  ggtitle("Kallisto PCA Plot")


png(paste(homedir, "/_PCA.png", sep=""), width=12,height=8, units='in', res=300)
p
dev.off()

p_2 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +   
  geom_point(size = 5, alpha = 0.7, position = "jitter") + 
  scale_color_manual(breaks = names(colour_blind_pallete), values=colour_blind_pallete) + 
  theme_classic() + 
  theme(legend.title=element_blank(), legend.text = element_text(size = 12), plot.title = element_text(size=22)) + 
  ggtitle("Salmon PCA Plot")

png(paste(homedir, "/_PCA_2.png", sep=""), width=12,height=8, units='in', res=300)
p_2
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
row.names(rld_cor) <- c("P1F2in", "P2F2in", "P3F2in", "P1F2out", "P2F2out", "P3F2out", 
                        "P1F3in", "P2F3in", "P3F3in", "P1F3out", "P2F3out", "P3F3out")

colnames(rld_cor) <- row.names(rld_cor)
png(paste(homedir, "/_HC_heatmap.png", sep=""), width=12,height=8, units='in', res=300)
print(pheatmap(rld_cor))
dev.off()


#Tsne 

set.seed(1)

tsne_results <- Rtsne(t(rld_mat), normalize = FALSE, num_threads = 0, check_duplicates = FALSE, perplexity = 2)

tsne_table <- as.data.frame(tsne_results$Y)
tsne_table$sample_type <- c("F2in", "F2in", "F2in", "F2out", "F2out", "F2out", "F3in", "F3in", "F3in", "F3out", "F3out", "F3out")
colnames(tsne_table) <- c("tsne1", "tsne2", "sample")
tsne_table$leaf <- rep(c("Int", "Int", "Int", "Ext", "Ext", "Ext"), 2)

tsne_1 <- ggplot(tsne_table, aes(x = tsne1, y = tsne2, colour = sample)) + 
  geom_point(size = 5, aes(shape = sample), position = "jitter") + 
  scale_color_manual(breaks = names(colour_blind_pallete), values=colour_blind_pallete) + 
  scale_shape_manual(values = c(0,1,2,4,3,6,7,8)) + 
  theme_classic() + 
  theme(legend.title=element_blank(), legend.text = element_text(size = 12), plot.title = element_text(size=22), axis.title = element_text(size = 16)) + 
  ggtitle("Salmon Tsne Plot")



png(paste(homedir, "/_tsne.png", sep=""), width=12,height=8, units='in', res=300)
tsne_1
dev.off()


tsne_2 <- ggplot(tsne_table, aes(x = tsne1, y = tsne2, colour = sample)) +   
  geom_point(size = 5, alpha = 0.9, aes(shape = leaf)) + 
  scale_color_manual(breaks = names(colour_blind_pallete), values=colour_blind_pallete) +
  scale_shape_manual(values = c(15,19)) + 
  theme_classic() + 
  theme(legend.title=element_blank(), legend.text = element_text(size = 12), plot.title = element_text(size=22), axis.title = element_text(size = 16), axis.text = element_text(size = 13)) + 
  ggtitle("tSNE plot all Nanopore samples\n")

png(paste(homedir, "/_tsne_2.png", sep=""), width=12,height=8, units='in', res=300)
tsne_2
dev.off()



##Plot prettier heatmap, more details
row.names(rld_cor) <- c("P1F2in", "P2F2in", "P3F2in", "P1F2out", "P2F2out", "P3F2out", 
                        "P1F3in", "P2F3in", "P3F3in", "P1F3out", "P2F3out", "P3F3out")



colnames(rld_cor) <- row.names(rld_cor)
png(paste(homedir, "/_HC_heatmap_2.png", sep=""), width=14,height=10, units='in', res=300)
print(pheatmap(rld_cor, main ="Heatmap Nanopore samples\n", fontsize = 15))
dev.off()
