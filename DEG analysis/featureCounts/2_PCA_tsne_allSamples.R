#### SCRIPT TO RUN DESEQ2 AND STORE RESULTS AND QC PLOTS ####
rm(list = ls())
setwd("C:/Users/Faculdade/Desktop/Dissertação/results 2020/mapping_to_genome/07_DESeq2/featureCounts")
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

#library load 
library('DESeq2')
library('ggplot2')
library('RColorBrewer')
library('dplyr')
library('pheatmap')
library('ggrepel')
library("Rtsne")



sampleNames <- c('P1F1in', 'P2F1in', 'P3F1in','P1F1out', 'P2F1out', 'P3F1out', 'P1F4in', 'P2F4in', 'P3F4in','P1F4out', 'P2F4out', 'P3F4out', 'idio_1', 'idio_2', 'idio_3', 'leaf_1', 'leaf_2', 'meso_1', 'meso_2', 'meso_3', 'totalptt_1', 'totalptt_2', 'totalptt_3')
sampleCondition <- c("f1in_f1in","f1in_f1in","f1in_f1in","f1out_f1out","f1out_f1out","f1out_f1out", "f4in_f4in","f4in_f4in","f4in_f4in","f4out_f4out","f4out_f4out","f4out_f4out", "idio", "idio", "idio", "leaf", "leaf", "meso", "meso", "meso", "totalptt", "totalptt", "totalptt")
sampleFiles <- c('F1P1_INT_genome_counts.txt', 'F1P2_INT_genome_counts.txt', 'F1P3_INT_genome_counts.txt','F1P1_EXT_genome_counts.txt', 'F1P2_EXT_genome_counts.txt', 'F1P3_EXT_genome_counts.txt','F4P1_INT_genome_counts.txt', 'F4P2_INT_genome_counts.txt', 'F4P3_INT_genome_counts.txt','F4P1_EXT_genome_counts.txt', 'F4P2_EXT_genome_counts.txt', 'F4P3_EXT_genome_counts.txt', 'idio_1_genome_counts.txt', 'idio_2_genome_counts.txt', 'idio_3_genome_counts.txt', 'leaves_1_genome_counts.txt', 'leaves_2_genome_counts.txt', 'meso_1_genome_counts.txt', 'meso_2_genome_counts.txt', 'meso_3_genome_counts.txt', 'totalptt_1_genome_counts.txt', 'totalptt_2_genome_counts.txt', 'totalptt_3_genome_counts.txt')


#### run deseq2 first on all data and explore data and normalized reads ####
# create dataframe that becomes a deseq table.
directory <- "counts"
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = as.factor(sampleCondition))
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ condition)


#filter 0 count rows
nrow(ddsHTSeq) #37320
keep <- rowSums(counts(ddsHTSeq)) > 0
ddsHTSeq <- ddsHTSeq[keep,]
nrow(ddsHTSeq) #24003


# Create DESeq2Dataset object
dds <- DESeq(ddsHTSeq)

rld <- rlog(dds, blind=TRUE) # betaPriorVar:  if missing is estimated from the data

homedir = "PCA_all_samples"
if (!dir.exists(homedir)){
  dir.create(homedir)
}

#Set a colour pallete 
colour_blind_pallete <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#ffcc33", "#0072B2", "#D55E00", "#CC79A7")
names(colour_blind_pallete) <-  c("f1in","f1out", "f4in","f4out", "idio", "leaf",  "meso", "totalptt")


## PCA
# Plot PCA
pca_data <- plotPCA(rld, intgroup="condition", returnData = T)
pca_data$group <- c("f1in", "f1in", "f1in", "f1out", "f1out", "f1out", "f4in", "f4in", "f4in", "f4out", "f4out", "f4out", "idio", "idio", "idio", "leaf", "leaf", "meso", "meso", "meso", "totalptt", "totalptt", "totalptt")

p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +   
  geom_point(size = 5, aes(shape = group), position = "jitter") + 
  scale_color_manual(breaks = names(colour_blind_pallete), values=colour_blind_pallete) + 
  scale_shape_manual(values = c(0,1,2,3,4,6,7,8)) + 
  theme_classic() + 
  theme(legend.title=element_blank(), legend.text = element_text(size = 12), plot.title = element_text(size=22)) + 
  ggtitle("FeatureCounts PCA Plot")
 

png(paste(homedir, "/_PCA.png", sep=""), width=12,height=8, units='in', res=500)
p
dev.off()

p_2 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +   
  geom_point(size = 5, alpha = 0.9) + 
  scale_color_manual(breaks = names(colour_blind_pallete), values=colour_blind_pallete) + 
  theme_classic() + 
  theme(legend.title=element_blank(), legend.text = element_text(size = 12), plot.title = element_text(size=22)) + 
  ggtitle("FeatureCounts PCA Plot")

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
  theme(legend.title=element_blank(), legend.text = element_text(size = 12), plot.title = element_text(size=22)) + 
  ggtitle("FeatureCounts Tsne Plot")



png(paste(homedir, "/_tsne.png", sep=""), width=12,height=8, units='in', res=300)
tsne_1
dev.off()


tsne_2 <- ggplot(tsne_table, aes(x = tsne1, y = tsne2, colour = sample)) +   
  geom_point(size = 5, alpha = 0.9) + 
  scale_color_manual(breaks = names(colour_blind_pallete), values=colour_blind_pallete) + 
  theme_classic() + 
  theme(legend.title=element_blank(), legend.text = element_text(size = 12), plot.title = element_text(size=22)) + 
  ggtitle("FeatureCounts tsne Plot")

png(paste(homedir, "/_tsne_2.png", sep=""), width=12,height=8, units='in', res=300)
tsne_2
dev.off()

