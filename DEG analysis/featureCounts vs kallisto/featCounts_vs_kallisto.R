#Use this script to compare kallisto and featureCounts analaysis and generate 
#figures 

#Misc 
setwd("C:/Users/Faculdade/Desktop/Dissertação/results 2020/mapping_to_genome/07_DESeq2/_featureCounts vs kallisto")
rm(list = ls())



#Load libraries
library('DESeq2')
library("tximport")
library("ggplot2")
library("gridExtra")

#######################################
###Load the data featureCounts Data ###
#######################################

sampleNames <- c('P1F1in_featCounts', 'P2F1in_featCounts', 'P3F1in_featCounts','P1F1out_featCounts', 'P2F1out_featCounts', 'P3F1out_featCounts', 'P1F4in_featCounts', 'P2F4in_featCounts', 'P3F4in_featCounts','P1F4out_featCounts', 'P2F4out_featCounts', 'P3F4out_featCounts', 'idio_1_featCounts', 'idio_2_featCounts', 'idio_3_featCounts', 'leaf_1_featCounts', 'leaf_2_featCounts', 'meso_1_featCounts', 'meso_2_featCounts', 'meso_3_featCounts', 'totalptt_1_featCounts', 'totalptt_2_featCounts', 'totalptt_3_featCounts')
sampleCondition <- c("f1in_f1in","f1in_f1in","f1in_f1in","f1out_f1out","f1out_f1out","f1out_f1out", "f4in_f4in","f4in_f4in","f4in_f4in","f4out_f4out","f4out_f4out","f4out_f4out", "idio", "idio", "idio", "leaf", "leaf", "meso", "meso", "meso", "totalptt", "totalptt", "totalptt")
sampleFiles <- c('F1P1_INT_genome_counts.txt', 'F1P2_INT_genome_counts.txt', 'F1P3_INT_genome_counts.txt','F1P1_EXT_genome_counts.txt', 'F1P2_EXT_genome_counts.txt', 'F1P3_EXT_genome_counts.txt','F4P1_INT_genome_counts.txt', 'F4P2_INT_genome_counts.txt', 'F4P3_INT_genome_counts.txt','F4P1_EXT_genome_counts.txt', 'F4P2_EXT_genome_counts.txt', 'F4P3_EXT_genome_counts.txt', 'idio_1_genome_counts.txt', 'idio_2_genome_counts.txt', 'idio_3_genome_counts.txt', 'leaves_1_genome_counts.txt', 'leaves_2_genome_counts.txt', 'meso_1_genome_counts.txt', 'meso_2_genome_counts.txt', 'meso_3_genome_counts.txt', 'totalptt_1_genome_counts.txt', 'totalptt_2_genome_counts.txt', 'totalptt_3_genome_counts.txt')


directory <- "featureCounts_counts"
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = as.factor(sampleCondition))
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ condition)

dds <- DESeq(ddsHTSeq)

#save the normalized table 
normalizedCounts <- counts(dds, normalized=TRUE) 
write.table(normalizedCounts, "FeatureCounts_normalizedCountsTable.txt", sep="\t", quote=F, col.names=NA)
#save sizeFactors
sizeFact <- sizeFactors(dds)
write.table(sizeFact, file = "featureCounts_sizeFactors.txt", sep="\t", quote=F, col.names=NA) # acrescentei

###########################
####Load kallisto Data#####
##########################
sampleNames <- c('F1P1_Int', 'F1P2_Int', 'F1P3_Int','F1P1_Ext', 'F1P2_Ext', 'F1P3_Ext', 'F4P1_Int', 'F4P2_Int', 'F4P3_Int','F4P1_Ext', 'F4P2_Ext', 'F4P3_Ext', 'idio_1', 'idio_2', 'idio_3', 'leaves_1', 'leaves_2', 'meso_1', 'meso_2', 'meso_3', 'totalptt_1', 'totalptt_2', 'totalptt_3')
sampleCondition <- c("f1in_f1in","f1in_f1in","f1in_f1in","f1out_f1out","f1out_f1out","f1out_f1out", "f4in_f4in","f4in_f4in","f4in_f4in","f4out_f4out","f4out_f4out","f4out_f4out", "idio", "idio", "idio", "leaf", "leaf", "meso", "meso", "meso", "totalptt", "totalptt", "totalptt")
sampleFiles <- "abundance.tsv"

tx2gene <- read.csv(file = "transcript2gene.tsv", sep = "\t", header = F)


#### run deseq2 first on all data and explore data and normalized reads ####
# create dataframe that becomes a deseq table.
directory <- "kallisto_counts/"
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = as.factor(sampleCondition))
files <-  paste("kallisto_counts", sampleTable$sampleName, sampleTable$fileName, sep = "/")
names(files) <- paste0(sampleTable$sampleName)

#Import counts using tximport 
kallisto_Counts <- tximport(files, type = "kallisto", tx2gene = tx2gene, dropInfReps = TRUE)

ddskallisto <- DESeqDataSetFromTximport(kallisto_Counts, colData = sampleTable, design = ~ condition)

dds <- DESeq(ddskallisto)

####Count Normalization####
#since in kallisto the size factor are computed by gene and by samples a matrix of gene x samples size factor are created 
#To get an idea of size factor the arithemic mean is calculated instead 
normalization_matrix <- normalizationFactors(dds)
estimated_size_fact <- colMeans(normalization_matrix)
write.table(estimated_size_fact, file = "kallisto_estsizeFactors.txt", sep="\t", quote=F, col.names=NA)

# To retrieve the normalized counts matrix from dds
normalizedCounts <- counts(dds, normalized=TRUE) 
# We can save this normalized data matrix to file for later use 
kallcolnames <- c('F1P1in_kallisto', 'F1P2in_kallisto', 'F1P3in_kallisto','F1P1out_kallisto', 'F1P2out_kallisto', 'F1P3out_kallisto', 'F4P1in_kallisto', 'F4P2in_kallisto', 'F4P3in_kallisto','F4P1out_kallisto', 'F4P2out_kallisto', 'F4P3out_kallisto', 'idio_1_kallisto', 'idio_2_kallisto', 'idio_3_kallisto', 'leaves_1_kallisto', 'leaves_2_kallisto', 'meso_1_kallisto', 'meso_2_kallisto', 'meso_3_kallisto', 'totalptt_1_kallisto', 'totalptt_2_kallisto', 'totalptt_3_kallisto')
write.table(normalizedCounts, "kallisto_normalizedCountsTable.txt", sep="\t", quote=F, col.names=kallcolnames) 
 
##########################################
###KALLISTO VS FEATURE COUNTS ANALYSIS####
##########################################
#Start by clear the work eviroment. If normalized count tables have already been generated 
#skip to this step 
rm(list = ls())

#Create results directories
results_dir <- "./results"
if (!(file.exists(results_dir))) {
  dir.create(results_dir)
}



#analyse size factors 
feacounts_size_factors <- read.table(file = "featureCounts_sizeFactors.txt")
kallisto_size_factors <- read.table(file = "kallisto_estsizeFactors.txt")

size_factors <- data.frame("kallisto" = kallisto_size_factors$x, "FeatCounts" = feacounts_size_factors$x ,row.names = rownames(kallisto_size_factors))

write.table(size_factors,file = paste(results_dir, "geral_size_factors.txt", sep = "/"), sep = "\t", quote = FALSE)

#compute correlation
size_factor_pearson <- cor(size_factors$kallisto, size_factors$FeatCounts, method = "pearson")
#0.9995009


#Analyse counts correlations
featCounts_normalized <- read.table(file = "FeatureCounts_normalizedCountsTable.txt")
kallisto_normalized <- read.table(file = "kallisto_normalizedCountsTable.txt")

merged_data_set <- merge(featCounts_normalized, kallisto_normalized, by = 0)
row.names(merged_data_set) <- merged_data_set$Row.names
merged_data_set <- merged_data_set[,-1]

##Compute correlations between samples analised by kallisto and by featureCounts 
#Leaf pair 1 
P1F1incor <- cor(merged_data_set$P1F1in_featCounts, merged_data_set$F1P1in_kallisto, method = "spearman")
P2F1incor <- cor(merged_data_set$P2F1in_featCounts, merged_data_set$F1P2in_kallisto, method = "spearman")
P3F1incor <- cor(merged_data_set$P3F1in_featCounts, merged_data_set$F1P3in_kallisto, method = "spearman")
P1F1extcor <- cor(merged_data_set$P1F1out_featCounts, merged_data_set$F1P1out_kallisto, method = "spearman")
P2F1extcor <- cor(merged_data_set$P2F1out_featCounts, merged_data_set$F1P2out_kallisto, method = "spearman")
P3F1extcor <- cor(merged_data_set$P3F1out_featCounts, merged_data_set$F1P3out_kallisto, method = "spearman")
#Leaf pair 4 
P1F4incor <- cor(merged_data_set$P1F4in_featCounts, merged_data_set$F4P1in_kallisto, method = "spearman")
P2F4incor <- cor(merged_data_set$P2F4in_featCounts, merged_data_set$F4P2in_kallisto, method = "spearman")
P3F4incor <- cor(merged_data_set$P3F4in_featCounts, merged_data_set$F4P3in_kallisto, method = "spearman")
P1F4extcor <- cor(merged_data_set$P1F4out_featCounts, merged_data_set$F4P1out_kallisto, method = "spearman")
P2F4extcor <- cor(merged_data_set$P2F4out_featCounts, merged_data_set$F4P2out_kallisto, method = "spearman")
P3F4extcor <- cor(merged_data_set$P3F4out_featCounts, merged_data_set$F4P3out_kallisto, method = "spearman")
#Idio
idio_1cor <- cor(merged_data_set$idio_1_featCounts, merged_data_set$idio_1_kallisto, method = "spearman")
idio_2cor <- cor(merged_data_set$idio_2_featCounts, merged_data_set$idio_2_kallisto, method = "spearman")
idio_3cor <- cor(merged_data_set$idio_3_featCounts, merged_data_set$idio_3_kallisto, method = "spearman")
#leaf 
leaf_1cor <- cor(merged_data_set$leaf_1_featCounts, merged_data_set$leaves_1_kallisto, method = "spearman")
leaf_2cor <- cor(merged_data_set$leaf_2_featCounts, merged_data_set$leaves_2_kallisto, method = "spearman")
#meso 
meso_1cor <- cor(merged_data_set$meso_1_featCounts, merged_data_set$meso_1_kallisto, method = "spearman")
meso_2cor <- cor(merged_data_set$meso_2_featCounts, merged_data_set$meso_2_kallisto, method = "spearman")
meso_3cor <- cor(merged_data_set$meso_3_featCounts, merged_data_set$meso_3_kallisto, method = "spearman")
#totalptt 
totalptt_1cor <- cor(merged_data_set$totalptt_1_featCounts, merged_data_set$totalptt_1_kallisto, method = "spearman")
totalptt_2cor <- cor(merged_data_set$totalptt_2_featCounts, merged_data_set$totalptt_2_kallisto, method = "spearman")
totalptt_3cor <- cor(merged_data_set$totalptt_3_featCounts, merged_data_set$totalptt_3_kallisto, method = "spearman")

#make spearman correlation table: 
samples <-  c("P1F1incor", "P2F1incor", "P3F1incor", "P1F1extcor", "P2F1extcor", "P3F1extcor", "P1F4incor", "P2F4incor", "P3F4incor", "P1F4extcor", "P2F4extcor", "P3F4extcor", "idio_1cor", "idio_2cor", "idio_3cor", "leaf_1cor", "leaf_2cor", "meso_1cor", "meso_2cor", "meso_3cor", "totalptt_1cor", "totalptt_2cor", "totalptt_3cor")
cor <- c(P1F1incor, P2F1incor, P3F1incor, P1F1extcor, P2F1extcor, P3F1extcor, P1F4incor, P2F4incor, P3F4incor, P1F4extcor, P2F4extcor, P3F4extcor, idio_1cor, idio_2cor, idio_3cor, leaf_1cor, leaf_2cor, meso_1cor, meso_2cor, meso_3cor, totalptt_1cor, totalptt_2cor, totalptt_3cor)

cor_data_frame <- data.frame("sample" = samples, "cor"  = cor)

write.table(cor_data_frame, file = paste(results_dir, "spearman_correlation_table.txt", sep = "/"), sep = "\t", quote = FALSE)

########################
###Create scatterplot###
########################
#for experience 1###

plot_1 <- ggplot(merged_data_set, aes(x = idio_1_featCounts, y= idio_1_kallisto)) + 
  geom_point(size = 0.1, colour = "red") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("idio 1 (featCounts)") + ylab("idio 1 (Kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12))

plot_2 <- ggplot(merged_data_set, aes(x = idio_2_featCounts, y= idio_2_kallisto)) + 
  geom_point(size = 0.1, colour = "red") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("idio 2 (FeatCounts)") + ylab("idio 2 (Kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

plot_3 <- ggplot(merged_data_set, aes(x = idio_3_featCounts, y= idio_3_kallisto)) + 
  geom_point(size = 0.1, colour = "red") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("idio 3 (FeatCounts)") + ylab("idio 3 (Kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

plot_4 <- ggplot(merged_data_set, aes(x = leaf_1_featCounts, y= leaves_1_kallisto)) + 
  geom_point(size = 0.1, colour = "red") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("leaf 1 (FeatCounts)") + ylab("leaf 1 (Kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

plot_5 <- ggplot(merged_data_set, aes(x = leaf_2_featCounts, y= leaves_2_kallisto)) + 
  geom_point(size = 0.1, colour = "red") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("leaf 2 (FeatCounts)") + ylab("leaf 2 (Kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

plot_6 <- ggplot(merged_data_set, aes(x = meso_1_featCounts, y= meso_1_kallisto)) + 
  geom_point(size = 0.1, colour = "red") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("meso 1 (FeatCounts)") + ylab("meso 1 (Kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

plot_7 <- ggplot(merged_data_set, aes(x = meso_2_featCounts, y= meso_2_kallisto)) + 
  geom_point(size = 0.1, colour = "red") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("meso 2 (FeatCounts)") + ylab("meso 2 (Kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

plot_8 <- ggplot(merged_data_set, aes(x = meso_3_featCounts, y= meso_3_kallisto)) + 
  geom_point(size = 0.1, colour = "red") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("meso 3 (FeatCounts)") + ylab("meso 3 (Kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

plot_9 <- ggplot(merged_data_set, aes(x = totalptt_1_featCounts, y= totalptt_1_kallisto)) + 
  geom_point(size = 0.1, colour = "red") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("totalptt 1 (FeatCounts)") + ylab("totalptt 1 (Kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

plot_10 <- ggplot(merged_data_set, aes(x = idio_2_featCounts, y= idio_2_kallisto)) + 
  geom_point(size = 0.1, colour = "red") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("totalptt 2 (FeatCounts)") + ylab("totalptt 2 (Kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

plot_11 <- ggplot(merged_data_set, aes(x = idio_3_featCounts, y= idio_3_kallisto)) + 
  geom_point(size = 0.1, colour = "red") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("totalptt 3 (FeatCounts)") + ylab("totalptt 3 (Kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

png(paste(results_dir, "scatterplots_experience_1_idio_leaf.png", sep = "/"), res = 300, width = 2100, height = 2100)
grid.arrange(plot_1, plot_4, plot_2, plot_5, plot_3, ncol = 2)
dev.off()

png(paste(results_dir, "scatterplots_experience_1_meso_totalptt.png", sep = "/"), res = 300, width = 2100, height = 2100)
grid.arrange(plot_6, plot_9, plot_7, plot_10, plot_8, plot_11, ncol = 2)
dev.off()


#for experience 2###

plot_1 <- ggplot(merged_data_set, aes(x = P1F1in_featCounts, y= F1P1in_kallisto)) + 
  geom_point(size = 0.1, colour = "green") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("P1F1I (FeatCounts)") + ylab("P1F1I (kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

plot_2 <- ggplot(merged_data_set, aes(x = P2F1in_featCounts, y= F1P2in_kallisto)) + 
  geom_point(size = 0.1, colour = "green") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("P2F1I (FeatCounts)") + ylab("P2F1I (kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

plot_3 <- ggplot(merged_data_set, aes(x = P3F1in_featCounts, y= F1P3in_kallisto)) + 
  geom_point(size = 0.1, colour = "green") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("P3F1I (FeatCounts)") + ylab("P3F1I (kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

plot_4 <- ggplot(merged_data_set, aes(x = P1F1out_featCounts, y= F1P1out_kallisto)) + 
  geom_point(size = 0.1, colour = "green") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("P1F1E (FeatCounts)") + ylab("P1F1E (kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

plot_5 <- ggplot(merged_data_set, aes(x = P2F1out_featCounts, y= F1P2out_kallisto)) + 
  geom_point(size = 0.1, colour = "green") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("P2F1E (FeatCounts)") + ylab("P2F1E (kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

plot_6 <- ggplot(merged_data_set, aes(x = P3F1out_featCounts, y= F1P3out_kallisto)) + 
  geom_point(size = 0.1, colour = "green") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("P3F1E (FeatCounts)") + ylab("P3F1E (kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

plot_7 <- ggplot(merged_data_set, aes(x = P1F4in_featCounts, y= F4P1in_kallisto)) + 
  geom_point(size = 0.1, colour = "green") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("P1F4I (FeatCounts)") + ylab("P1F4I (kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

plot_8 <- ggplot(merged_data_set, aes(x = P2F4in_featCounts, y= F4P2in_kallisto)) + 
  geom_point(size = 0.1, colour = "green") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("P2F4I (FeatCounts)") + ylab("P2F4I (kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

plot_9 <- ggplot(merged_data_set, aes(x = P3F4in_featCounts, y= F4P3in_kallisto)) + 
  geom_point(size = 0.1, colour = "green") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("P3F4I (FeatCounts)") + ylab("P3F4I (kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

plot_10 <- ggplot(merged_data_set, aes(x = P1F4out_featCounts, y= F4P1out_kallisto)) + 
  geom_point(size = 0.1, colour = "green") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("P1F4E (FeatCounts)") + ylab("P1F4E (kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

plot_11 <- ggplot(merged_data_set, aes(x = P2F4out_featCounts, y= F4P2out_kallisto)) + 
  geom_point(size = 0.1, colour = "green") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("P2F4E (FeatCounts)") + ylab("P2F4E (kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

plot_12 <- ggplot(merged_data_set, aes(x = P3F4out_featCounts, y= F4P3out_kallisto)) + 
  geom_point(size = 0.1, colour = "green") + 
  geom_smooth(method = "lm", colour = "black") + xlim (0,2500) + ylim(0,2500) + 
  xlab("P3F4E (FeatCounts)") + ylab("P3F4E (kallisto)") + theme_classic() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) 

png(paste(results_dir, "scatterplots_experience_2_leaf_pair_1.png", sep = "/"), res = 300, width = 2100, height = 2100)
grid.arrange(plot_1, plot_4, plot_2, plot_5, plot_3, plot_6, ncol = 2)
dev.off()

png(paste(results_dir, "scatterplots_experience_2_leaf_pair_4.png", sep = "/"), res = 300, width = 2100, height = 2100)
grid.arrange(plot_7, plot_10, plot_8, plot_11, plot_9, plot_12, ncol = 2)
dev.off()

################################
### foldChange scatterplots ####
################################

###Meso vs idio###
meso_vs_idio_featCounts_table <- read.table(file = "meso_vs_idio_featureCounts_Res.csv", sep = ",", header = TRUE, row.names = 1)
meso_vs_idio_kallisto_table <- read.table(file = "meso_vs_idio_kallisto_Res.csv", sep = ",", header = TRUE, row.names = 1)

#remove row with padj.value = NA and merge datasets 
meso_vs_idio_featCounts_table <- meso_vs_idio_featCounts_table[complete.cases(meso_vs_idio_featCounts_table), ]
meso_vs_idio_kallisto_table <- meso_vs_idio_kallisto_table[complete.cases(meso_vs_idio_kallisto_table), ]
colnames(meso_vs_idio_featCounts_table) <- c("baseMean_featCounts", "log2FoldChange_featCounts", "lfcSE_featCounts", "stat_featCounts", "pvalue_featCounts", "padj_featCounts")
colnames(meso_vs_idio_kallisto_table) <- c("baseMean_kallisto", "log2FoldChange_kallisto", "lfcSE_kallisto", "stat_kallisto", "pvalue_kallisto", "padj_kallisto")

merged_data_set_meso_vs_idio <- merge(meso_vs_idio_featCounts_table, meso_vs_idio_kallisto_table, by = 0)
row.names(merged_data_set_meso_vs_idio) <- merged_data_set_meso_vs_idio$Row.names
merged_data_set_meso_vs_idio <- merged_data_set_meso_vs_idio[,c(-1)]

#upregulated_interesect_3_featureCounts <- read.csv("Venn_geral_UP_list_featureCounts.csv", row.names= 1) 
#upregulated_interesect_3_kallisto <- read.csv("Venn_geral_UP_list_kallisto.csv", row.names = 1) 
#label_set <- union(upregulated_interesect_3_featureCounts$x, upregulated_interesect_3_kallisto$x)

#calculated wald statistic and logfold change correlations across samples 
cor_lfc_meso_vs_idio <- cor(merged_data_set_meso_vs_idio$log2FoldChange_featCounts, merged_data_set_meso_vs_idio$log2FoldChange_kallisto, method = "spearman")
cor_wald_meso_vs_idio <- cor(merged_data_set_meso_vs_idio$stat_featCounts, merged_data_set_meso_vs_idio$stat_kallisto, method = "spearman")



#plot stat (wald statistic in this case)
#add this comand in order to show some genes
#geom_point(data = subset(merged_data_set_meso_vs_idio, rownames(merged_data_set_meso_vs_idio) %in% label_set), size = 1, colour = "red")

plot_wald_stat_meso_vs_idio <- ggplot(merged_data_set_meso_vs_idio, aes(x = stat_featCounts, y = stat_kallisto)) + 
  geom_point(size = 0.1, colour = "blue") + 
  xlab("Wald Statistic Feature Counts") + ylab("Wald Statistic Kallisto") + 
  theme_classic() + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16)) 
  
  
plot_fold_chage_meso_vs_idio <- ggplot(merged_data_set_meso_vs_idio, aes(x = log2FoldChange_featCounts, y = log2FoldChange_kallisto)) + 
  geom_point(size = 0.1, colour = "blue") + 
  xlab("log2FoldChange Feature Counts") + ylab("log2FoldChange Kallisto") + 
  theme_classic() + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16)) 


png(paste(results_dir, "scatterplots_meso_vs_idio.png", sep = "/"), res = 300, width = 1800, height = 1800)
grid.arrange(plot_wald_stat_meso_vs_idio, plot_fold_chage_meso_vs_idio, ncol = 1)
dev.off()


###F1in vs F1out### 
F1in_vs_F1out_featCounts_table <- read.table(file = "F1in_vs_F1out_featureCounts_Res.csv", sep = ",", header = TRUE, row.names = 1)
F1in_vs_F1out_kallisto_table <- read.table(file = "F1in_vs_F1out_kallisto_Res.csv", sep = ",", header = TRUE, row.names = 1)

#remove row with padj.value = NA and merge datasets 
F1in_vs_F1out_featCounts_table <- F1in_vs_F1out_featCounts_table[complete.cases(F1in_vs_F1out_featCounts_table), ]
F1in_vs_F1out_kallisto_table <- F1in_vs_F1out_kallisto_table[complete.cases(F1in_vs_F1out_kallisto_table), ]
colnames(F1in_vs_F1out_featCounts_table) <- c("baseMean_featCounts", "log2FoldChange_featCounts", "lfcSE_featCounts", "stat_featCounts", "pvalue_featCounts", "padj_featCounts")
colnames(F1in_vs_F1out_kallisto_table) <- c("baseMean_kallisto", "log2FoldChange_kallisto", "lfcSE_kallisto", "stat_kallisto", "pvalue_kallisto", "padj_kallisto")

merged_data_set_F1in_vs_F1out <- merge(F1in_vs_F1out_featCounts_table, F1in_vs_F1out_kallisto_table, by = 0)
row.names(merged_data_set_F1in_vs_F1out) <- merged_data_set_F1in_vs_F1out$Row.names
merged_data_set_F1in_vs_F1out <- merged_data_set_F1in_vs_F1out[,c(-1)]

#upregulated_interesect_3_featureCounts <- read.csv("Venn_geral_UP_list_featureCounts.csv", row.names= 1) 
#upregulated_interesect_3_kallisto <- read.csv("Venn_geral_UP_list_kallisto.csv", row.names = 1) 
#label_set <- union(upregulated_interesect_3_featureCounts$x, upregulated_interesect_3_kallisto$x)

#calculated wald statistic and logfold change correlations across samples 
cor_lfc_f1in_vs_f1out <- cor(merged_data_set_F1in_vs_F1out$log2FoldChange_featCounts, merged_data_set_F1in_vs_F1out$log2FoldChange_kallisto, method = "spearman")
cor_wald_f1in_vs_f1out <- cor(merged_data_set_F1in_vs_F1out$stat_featCounts, merged_data_set_F1in_vs_F1out$stat_kallisto, method = "spearman")


#plot stat (wald statistic in this case)
plot_wald_stat_F1in_vs_F1out <- ggplot(merged_data_set_F1in_vs_F1out, aes(x = stat_featCounts, y = stat_kallisto)) + 
  geom_point(size = 0.1, colour = "blue") + 
  xlab("Wald Statistic Feature Counts") + ylab("Wald Statistic Kallisto") +
  theme_classic() + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12)) 

plot_fold_chage_F1in_vs_F1out <- ggplot(merged_data_set_F1in_vs_F1out, aes(x = log2FoldChange_featCounts, y = log2FoldChange_kallisto)) + 
  geom_point(size = 0.1, colour = "blue") + 
  xlab("log2FoldChange Feature Counts") + ylab("log2FoldChange Kallisto") +
  theme_classic() + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12)) 


png(paste(results_dir, "scatterplots_F1in_vs_F1out.png", sep = "/"), res = 300, width = 1800, height = 1800)
grid.arrange(plot_wald_stat_F1in_vs_F1out, plot_fold_chage_F1in_vs_F1out, ncol = 1)
dev.off()



###F4in vs F4out### 
F4in_vs_F4out_featCounts_table <- read.table(file = "F4in_vs_F4out_featureCounts_Res.csv", sep = ",", header = TRUE, row.names = 1)
F4in_vs_F4out_kallisto_table <- read.table(file = "F4in_vs_F4out_kallisto_Res.csv", sep = ",", header = TRUE, row.names = 1)

#remove row with padj.value = NA and merge datasets 
F4in_vs_F4out_featCounts_table <- F4in_vs_F4out_featCounts_table[complete.cases(F4in_vs_F4out_featCounts_table), ]
F4in_vs_F4out_kallisto_table <- F4in_vs_F4out_kallisto_table[complete.cases(F4in_vs_F4out_kallisto_table), ]
colnames(F4in_vs_F4out_featCounts_table) <- c("baseMean_featCounts", "log2FoldChange_featCounts", "lfcSE_featCounts", "stat_featCounts", "pvalue_featCounts", "padj_featCounts")
colnames(F4in_vs_F4out_kallisto_table) <- c("baseMean_kallisto", "log2FoldChange_kallisto", "lfcSE_kallisto", "stat_kallisto", "pvalue_kallisto", "padj_kallisto")

merged_data_set_F4in_vs_F4out <- merge(F4in_vs_F4out_featCounts_table, F4in_vs_F4out_kallisto_table, by = 0)
row.names(merged_data_set_F4in_vs_F4out) <- merged_data_set_F4in_vs_F4out$Row.names
merged_data_set_F4in_vs_F4out <- merged_data_set_F4in_vs_F4out[,c(-1)]

#upregulated_interesect_3_featureCounts <- read.csv("Venn_geral_UP_list_featureCounts.csv", row.names= 1) 
#upregulated_interesect_3_kallisto <- read.csv("Venn_geral_UP_list_kallisto.csv", row.names = 1) 
#label_set <- union(upregulated_interesect_3_featureCounts$x, upregulated_interesect_3_kallisto$x)

#calculated wald statistic and logfold change correlations across samples 
cor_lfc_F4in_vs_F4out <- cor(merged_data_set_F4in_vs_F4out$log2FoldChange_featCounts, merged_data_set_F4in_vs_F4out$log2FoldChange_kallisto, method = "spearman")
cor_wald_F4in_vs_F4out <- cor(merged_data_set_F4in_vs_F4out$stat_featCounts, merged_data_set_F4in_vs_F4out$stat_kallisto, method = "spearman")


#plot stat (wald statistic in this case)
plot_wald_stat_F4in_vs_F4out <- ggplot(merged_data_set_F4in_vs_F4out, aes(x = stat_featCounts, y = stat_kallisto)) + 
  geom_point(size = 0.1, colour = "blue") + 
  xlab("Wald Statistic Feature Counts") + ylab("Wald Statistic Kallisto") +
  theme_classic() + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12)) 

plot_fold_chage_F4in_vs_F4out <- ggplot(merged_data_set_F4in_vs_F4out, aes(x = log2FoldChange_featCounts, y = log2FoldChange_kallisto)) + 
  geom_point(size = 0.1, colour = "blue") + 
  xlab("log2FoldChange Feature Counts") + ylab("log2FoldChange Kallisto") +
  theme_classic() + theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12)) 


png(paste(results_dir, "scatterplots_F4in_vs_F4out.png", sep = "/"), res = 300, width = 1800, height = 1800)
grid.arrange(plot_wald_stat_F4in_vs_F4out, plot_fold_chage_F4in_vs_F4out, ncol = 1)
dev.off()



###Make file with spearman correlation of wald stat and log2foldchange 
waldstatcor <- c(cor_wald_meso_vs_idio, cor_wald_f1in_vs_f1out, cor_wald_F4in_vs_F4out)
lfccor <- c(cor_lfc_meso_vs_idio, cor_lfc_f1in_vs_f1out, cor_lfc_F4in_vs_F4out)

wald_lf_co.df <- data.frame("WaldStat" = waldstatcor, "lfccor" = lfccor, row.names = c("meso_vs_idio", "F1in vs f1out", "F4in vs F4out"))
write.csv(wald_lf_co.df, file = "AFterDeSeq2Cor.txt", quote = F)

