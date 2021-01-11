#misc
setwd("C:/Users/Faculdade/Desktop/Dissertação/results 2020/mapping_to_genome/annotation/03_dammit")
rm(list = ls())

#Load libraries 
library(xlsx)
library(gdata)
options(java.parameters = "-Xmx2048m")  

#Load datasets to make the table ---- 
dammit_annotation <- read.delim("StepsToFinalAnnotation/dammit_annotation_final.tsv", sep = "\t", row.names = 1)
meso_vs_idio_up_list <- read.delim("DESeq2_tables/meso_vs_idio_kallisto_Res.csv", sep = ",", row.names = 1)
f1in_vs_f1out_up_list <- read.delim("DESeq2_tables/f1in_vs_f1out_kallisto_Res.csv", sep = ",", row.names = 1)
f4in_vs_f4out_up_list <- read.delim("DESeq2_tables/f4in_vs_f4out_kallisto_Res.csv", sep = ",", row.names = 1)

idio_up_gene_list = as.character(read.delim(file = "DESeq2_tables/meso_vs_idio_kallisto_sigUP.csv", sep = ",")[,2])
f1_up_gene_list = as.character(read.delim(file = "DESeq2_tables/f1in_vs_f1out_kallisto_sigUP.csv", sep = ",")[,2])
f4_up_gene_list = as.character(read.delim(file = "DESeq2_tables/f4in_vs_f4out_kallisto_sigUP.csv", sep = ",")[,2])

#get expression table and make a an average dataset
kallisto_expression_table <- read.delim(file = "DESeq2_tables/kallisto_normalizedCountsTable.txt")
kallisto_expression_table <- data.frame("avgF1in" = rowMeans(kallisto_expression_table[1:3]), "avgF1out" = rowMeans(kallisto_expression_table[4:6]), "avgF4in" = rowMeans(kallisto_expression_table[7:9]), "avgF4out" = rowMeans(kallisto_expression_table[10:12]), "avgMeso" = rowMeans(kallisto_expression_table[18:20]), "avgIdio" = rowMeans(kallisto_expression_table[13:15]))


#get new gene list 
new_gene_bool_vector <- ifelse(startsWith(row.names(dammit_annotation), "CATHA."), T, F)
new_gene_list <- row.names(dammit_annotation[new_gene_bool_vector,])

#Make table  ----- 
#subset expression levels 
new_gene_table <- subset(kallisto_expression_table, row.names(kallisto_expression_table) %in% new_gene_list)

#add collums 
new_gene_table$gene_name <- c()
new_gene_table$LogFolchangeIdio <- c()
new_gene_table$LogFolchangeF1out <- c()
new_gene_table$LogFolchangeF4out<- c()
new_gene_table$FolchangeIdio <- c()
new_gene_table$FolchangeF1out <- c()
new_gene_table$FolchangeF4out<- c()



for (i in 1:length(new_gene_table$avgF1in)){
  new_gene_table$gene_name[i] <- as.character(dammit_annotation[row.names(new_gene_table)[i], 3])
  new_gene_table$LogFolchangeIdio[i] <- meso_vs_idio_up_list[row.names(new_gene_table)[i], 2]
  new_gene_table$LogFolchangeF1Out[i] <- f1in_vs_f1out_up_list[row.names(new_gene_table)[i], 2]
  new_gene_table$LogFolchangeF4Out[i] <- f4in_vs_f4out_up_list[row.names(new_gene_table)[i], 2]
  new_gene_table$FoldChange_idio[i] <- new_gene_table$avgIdio[i] / new_gene_table$avgMeso[i]
  new_gene_table$FolchangeF1out[i] <- new_gene_table$avgF1in[i] / new_gene_table$avgF1out[i]
  new_gene_table$FolchangeF4out[i] <- new_gene_table$avgF4in[i] / new_gene_table$avgF4out[i]
}


new_gene_table$UpIdio <- ifelse(row.names(new_gene_table) %in% idio_up_gene_list, T, F)
new_gene_table$UpF1out <- ifelse(row.names(new_gene_table) %in% f1_up_gene_list, T, F)
new_gene_table$UpF4out <- ifelse(row.names(new_gene_table) %in% f4_up_gene_list, T, F)

new_gene_table$proteome_location <- dammit_annotation[as.character(row.names(new_gene_table)),15]
new_gene_table$cathacyc_annotation <- dammit_annotation[as.character(row.names(new_gene_table)),16]

#save the table 
write.xlsx(new_gene_table, file = "results/new_gene_annotated.xlsx", sheetName = "Folha 1")
