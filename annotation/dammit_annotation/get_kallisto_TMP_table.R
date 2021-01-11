#Set working directory
setwd("C:/Users/Faculdade/Desktop/Dissertação/results 2020/mapping_to_genome/annotation/03_dammit/kallisto_TPM")
#misc 
rm(list= ls())

#Load libraries 
library(tximport)


#Parse data from the tsv file with tximport 
sampleNames <- c('F1P1_Int', 'F1P2_Int', 'F1P3_Int','F1P1_Ext', 'F1P2_Ext', 'F1P3_Ext','F4P1_Int', 'F4P2_Int', 'F4P3_Int','F4P1_Ext', 'F4P2_Ext', 'F4P3_Ext', 'idio_1', 'idio_2', 'idio_3', 'meso_1', 'meso_2', 'meso_3', 'leaves_1', 'leaves_2', 'totalptt_1', 'totalptt_2', 'totalptt_3')
sampleCondition <- c("f1in_f1in","f1in_f1in","f1in_f1in","f1out_f1out","f1out_f1out","f1out_f1out", "f4in_f4in","f4in_f4in","f4in_f4in","f4out_f4out","f4out_f4out","f4out_f4out", 'idio', 'idio', 'idio', 'meso', 'meso', 'meso', 'leaves', 'leaves', 'totalptt', 'totalptt', 'totalptt')
sampleFiles <- "abundance.tsv"



tx2gene <- read.csv(file = "transcript2gene.tsv", sep = "\t", header = F)
directory <- "C:/Users/Faculdade/Desktop/Dissertação/results 2020/mapping_to_genome/annotation/03_dammit/kallisto_TPM/"
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
files <- paste("counts", sampleTable$sampleName, sampleTable$fileName, sep = "/")
names(files) <- paste0(sampleTable$sampleName)

kallisto_Counts <- tximport(files, type="kallisto", tx2gene = tx2gene, dropInfReps = TRUE)


abundances.df <- data.frame(kallisto_Counts[["abundance"]])

write.table(abundances.df, "kallisto_TPM.tsv", quote = FALSE, sep = "\t")


