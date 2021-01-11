#### SCRIPT TO RUN DESEQ2 AND STORE RESULTS AND QC PLOTS ####
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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


homedir = "PCA_all_samples"
if (!dir.exists(homedir)){
  dir.create(homedir)
}




sampleNames <- c('F1P1_Int', 'F1P2_Int', 'F1P3_Int','F1P1_Ext', 'F1P2_Ext', 'F1P3_Ext', 'F4P1_Int', 'F4P2_Int', 'F4P3_Int','F4P1_Ext', 'F4P2_Ext', 'F4P3_Ext', 'idio_1', 'idio_2', 'idio_3', 'leaves_1', 'leaves_2', 'meso_1', 'meso_2', 'meso_3', 'totalptt_1', 'totalptt_2', 'totalptt_3')
sampleCondition <- c("f1in_f1in","f1in_f1in","f1in_f1in","f1out_f1out","f1out_f1out","f1out_f1out", "f4in_f4in","f4in_f4in","f4in_f4in","f4out_f4out","f4out_f4out","f4out_f4out", "idio", "idio", "idio", "leaf", "leaf", "meso", "meso", "meso", "totalptt", "totalptt", "totalptt")
sampleFiles <- "abundance.tsv"

tx2gene <- read.csv(file = "transcript2gene.tsv", sep = "\t", header = F)


#### run deseq2 first on all data and explore data and normalized reads ####
# create dataframe that becomes a deseq table.
directory <- "counts/"
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = as.factor(sampleCondition))
files <-  paste("counts", sampleTable$sampleName, sampleTable$fileName, sep = "/")
names(files) <- paste0(sampleTable$sampleName)

kallisto_Counts <- tximport(files, type = "kallisto", tx2gene = tx2gene, dropInfReps = TRUE)


ddkallisto <- DESeqDataSetFromTximport(kallisto_Counts, colData = sampleTable, design = ~ condition)


#filter 0 count rows
nrow(ddkallisto) #37320
keep <- rowSums(counts(ddkallisto)) > 0
ddkallisto <- ddkallisto[keep,]
nrow(ddkallisto) #28332


# Create DESeq2Dataset object
dds <- DESeq(ddkallisto)
# To retrieve the normalized counts matrix from dds
normalizedCounts <- counts(dds, normalized=TRUE) 

# We can save this normalized data matrix to file for later use
write.table(normalizedCounts, file = paste(homedir,"/normalizedCountsTable.txt", sep=""), sep="\t", quote=F, col.names=NA) # acrescentei


### Sample-level QC
# To improve the distances/clustering for the PCA and hierarchical clustering visualization methods,
# we need to moderate the variance across the mean by applying the rlog transformation to the normalized counts
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE) # betaPriorVar:  if missing is estimated from the data



#Set a colour pallete 
colour_blind_pallete <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#ffcc33", "#0072B2", "#D55E00", "#CC79A7")
names(colour_blind_pallete) <-  c("f1in","f1out", "f4in","f4out", "idio", "leaf",  "meso", "totalptt")

# Plot PCA
pca_data <- plotPCA(rld, intgroup="condition", returnData = T)
#added to match the plots from featureCounts
pca_data$group <- c("f1in", "f1in", "f1in", "f1out", "f1out", "f1out", "f4in", "f4in", "f4in", "f4out", "f4out", "f4out", "idio", "idio", "idio", "leaf", "leaf", "meso", "meso", "meso", "totalptt", "totalptt", "totalptt")

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
  ggtitle("Kallisto PCA Plot")

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
row.names(rld_cor) <- c("P1F1in", "P2F1in", "P3F1in", "P1F1out", "P2F1out", "P3F1out", 
                        "P1F4in", "P2F4in", "P3F4in", "P1F4out", "P2F4out", "P3F4out", 
                        "idio_1", "idio_2", "idio_3", "leaf_1", "leaf_2",
                        "meso_1", "meso_2", "meso_3", "totalptt_1", "totalptt_2", "totalptt_3")
colnames(rld_cor) <- row.names(rld_cor)
png(paste(homedir, "/_HC_heatmap.png", sep=""), width=12,height=8, units='in', res=300)
print(pheatmap(rld_cor))
dev.off()


#Tsne 

set.seed(1)

tsne_results <- Rtsne(t(rld_mat), normalize = FALSE, num_threads = 0, check_duplicates = FALSE, perplexity = 2)

tsne_table <- as.data.frame(tsne_results$Y)
tsne_table$sample_type <- c("f1in", "f1in", "f1in", "f1out", "f1out", "f1out", "f4in", "f4in", "f4in", "f4out", "f4out", "f4out", "idio", "idio", "idio", "leaf", "leaf", "meso", "meso", "meso", "totalptt", "totalptt", "totalptt")
colnames(tsne_table) <- c("tsne1", "tsne2", "sample")

tsne_1 <- ggplot(tsne_table, aes(x = tsne1, y = tsne2, colour = sample)) + 
  geom_point(size = 5, aes(shape = sample), position = "jitter") + 
  scale_color_manual(breaks = names(colour_blind_pallete), values=colour_blind_pallete) + 
  scale_shape_manual(values = c(0,1,2,4,3,6,7,8)) + 
  theme_classic() + 
  theme(legend.title=element_blank(), legend.text = element_text(size = 12), plot.title = element_text(size=22), axis.title = element_text(size = 16)) + 
  ggtitle("kallisto Tsne Plot")



png(paste(homedir, "/_tsne.png", sep=""), width=12,height=8, units='in', res=300)
tsne_1
dev.off()


tsne_2 <- ggplot(tsne_table, aes(x = tsne1, y = tsne2, colour = sample)) +   
  geom_point(size = 5, alpha = 0.9) + 
  scale_color_manual(breaks = names(colour_blind_pallete), values=colour_blind_pallete) + 
  theme_classic() + 
  theme(legend.title=element_blank(), legend.text = element_text(size = 12), plot.title = element_text(size=22), axis.title = element_text(size = 16), axis.text = element_text(size = 13)) + 
  ggtitle("tSNE plot all Illumina samples\n")

png(paste(homedir, "/_tsne_2.png", sep=""), width=12,height=8, units='in', res=300)
tsne_2
dev.off()



##Plot prettier heatmap, more details
row.names(rld_cor) <- c("P1F1in", "P2F1in", "P3F1in", "P1F1out", "P2F1out", "P3F1out", 
                        "P1F4in", "P2F4in", "P3F4in", "P1F4out", "P2F4out", "P3F4out", 
                        "idio 1", "idio 2", "idio 3", "leaf 1", "leaf 2",
                        "meso 1", "meso 2", "meso 3", "totalptt 1", "totalptt 2", "totalptt 3")


colnames(rld_cor) <- row.names(rld_cor)
png(paste(homedir, "/_HC_heatmap_2.png", sep=""), width=14,height=10, units='in', res=300)
print(pheatmap(rld_cor, main ="Heatmap Illumina samples\n", fontsize = 15))
dev.off()



#plot PRX1 
