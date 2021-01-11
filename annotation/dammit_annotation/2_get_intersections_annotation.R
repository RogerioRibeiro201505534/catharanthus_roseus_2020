#get working directory 
setwd("C:/Users/Faculdade/Desktop/Dissertação/results 2020/mapping_to_genome/annotation/03_dammit")

#misc 
rm(list = ls())

#Load libraries 
library(bio3d)
library("xlsx")


##output dir  
results_dir <- "./results"
if (!(file.exists(results_dir))) {
  dir.create(results_dir)
}

#Read dammit annot output 
dammit_annot <- read.delim("StepsToFinalAnnotation/dammit_annotation_final.tsv", sep = "\t", row.names = 1)

#Read p.values and expression tables 
idio_res_kallisto <- read.csv(file = "DESeq2_tables/meso_vs_idio_kallisto_Res.csv", row.names = 1)
f1out_res_kallisto <- read.csv(file = "DESeq2_tables/F1in_vs_F1out_kallisto_Res.csv", row.names = 1)
f4out_res_kallisto <- read.csv(file = "DESeq2_tables/F4in_vs_F4out_kallisto_Res.csv", row.names = 1)


#Get pvalues and logFoldChange from the DESeq2_tables 
idio_kallisto_foldchange <- idio_res_kallisto$log2FoldChange
idio_kallisto_padj <- idio_res_kallisto$padj
names(idio_kallisto_foldchange) <- row.names(idio_res_kallisto)
names(idio_kallisto_padj) <- row.names(idio_res_kallisto)

f1out_kallisto_foldchange <- f1out_res_kallisto$log2FoldChange
f1out_kallisto_padj <- f1out_res_kallisto$padj
names(f1out_kallisto_foldchange) <- row.names(f1out_res_kallisto)
names(f1out_kallisto_padj) <- row.names(f1out_res_kallisto)

f4out_kallisto_foldchange <- f4out_res_kallisto$log2FoldChange
f4out_kallisto_padj <- f4out_res_kallisto$padj
names(f4out_kallisto_foldchange) <- row.names(f4out_res_kallisto)
names(f4out_kallisto_padj) <- row.names(f4out_res_kallisto)

####################################
## Intersection of 3 conditions ####
####################################

###create output dir  
results_dir_2 <- paste(results_dir, "intersect_3", sep = "/")
if (!(file.exists(results_dir_2))) {
  dir.create(results_dir_2)
}


#### get annotation for the upregulated gene ####
Venn_geral_genes_kallisto <- read.delim("intersection_tables/Venn_geral_UP_list.csv", sep = ",", row.names = 1)
Venn_geral_dammit <- subset(dammit_annot, row.names(dammit_annot) %in% Venn_geral_genes_kallisto$x)


#Add the logfold change and p.value to the table 

Venn_geral_dammit$foldChange_idio_kallisto <- c()
Venn_geral_dammit$padj_idio_kallisto <- c()
Venn_geral_dammit$foldChange_f1_kallisto <- c()
Venn_geral_dammit$padj_f1_kallisto <- c()
Venn_geral_dammit$foldChange_f4_kallisto <- c()
Venn_geral_dammit$padj_f4_kallisto <- c()


j <- 0
for (i in row.names(Venn_geral_dammit)){
  j <- j + 1
  Venn_geral_dammit$foldChange_idio_kallisto[j] <- idio_kallisto_foldchange[i]
  Venn_geral_dammit$padj_idio_kallisto[j] <- idio_kallisto_padj[i]
  Venn_geral_dammit$foldChange_f1_kallisto[j] <- f1out_kallisto_foldchange[i]
  Venn_geral_dammit$padj_f1_kallisto[j] <- f1out_kallisto_padj[i]
  Venn_geral_dammit$foldChange_f4_kallisto[j] <- f4out_kallisto_foldchange[i]
  Venn_geral_dammit$padj_f4_kallisto[j] <- f4out_kallisto_padj[i]
}


#write the output table
write.xlsx(Venn_geral_dammit, file = paste(results_dir_2, "dammit_annotation_geral_kallisto.xlsx", sep = "/"), showNA = F)



##############################
### F1 vs F4 intersection ####
##############################

###create output dir  
results_dir_2 <- paste(results_dir, "F1_F4_intersect", sep = "/")
if (!(file.exists(results_dir_2))) {
  dir.create(results_dir_2)
}




#### get annotation for the upregulated genes ####
Venn_F1_F4_genes_kallisto <- read.delim("intersection_tables/Venn_F1out_F4Out_UP_list.csv", row.names = 1, sep = ",")


Venn_F1_F4_dammit <- subset(dammit_annot, row.names(dammit_annot) %in% Venn_F1_F4_genes_kallisto$x)


#Add the logfold change and p.value to the table 

Venn_F1_F4_dammit$foldChange_idio_kallisto <- c()
Venn_F1_F4_dammit$padj_idio_kallisto <- c()
Venn_F1_F4_dammit$foldChange_f1_kallisto <- c()
Venn_F1_F4_dammit$padj_f1_kallisto <- c()
Venn_F1_F4_dammit$foldChange_f4_kallisto <- c()
Venn_F1_F4_dammit$padj_f4_kallisto <- c()


j <- 0
for (i in row.names(Venn_F1_F4_dammit)){
  j <- j + 1
  Venn_F1_F4_dammit$foldChange_idio_kallisto[j] <- idio_kallisto_foldchange[i]
  Venn_F1_F4_dammit$padj_idio_kallisto[j] <- idio_kallisto_padj[i]
  Venn_F1_F4_dammit$foldChange_f1_kallisto[j] <- f1out_kallisto_foldchange[i]
  Venn_F1_F4_dammit$padj_f1_kallisto[j] <- f1out_kallisto_padj[i]
  Venn_F1_F4_dammit$foldChange_f4_kallisto[j] <- f4out_kallisto_foldchange[i]
  Venn_F1_F4_dammit$padj_f4_kallisto[j] <- f4out_kallisto_padj[i]
  
}


#write the output table
write.xlsx(Venn_F1_F4_dammit, file = paste(results_dir_2, "dammit_annotation_F1_F4.xlsx", sep = "/"), showNA = T)




################################
### idio vs F4 intersection ####
################################

###create output dir  
results_dir_2 <- paste(results_dir, "idio_F4_intersect", sep = "/")
if (!(file.exists(results_dir_2))) {
  dir.create(results_dir_2)
}


#### get annotation for the upregulated genes ####
Venn_idio_F4_genes_kallisto <- read.delim("intersection_tables/Venn_Idio_F4Out_UP_list.csv", row.names = 1, sep = ",")

Venn_idio_F4_dammit <- subset(dammit_annot, row.names(dammit_annot) %in% Venn_idio_F4_genes_kallisto$x)


#Add the logfold change and p.value to the table 

Venn_idio_F4_dammit$foldChange_idio_kallisto <- c()
Venn_idio_F4_dammit$padj_idio_kallisto <- c()
Venn_idio_F4_dammit$foldChange_f1_kallisto <- c()
Venn_idio_F4_dammit$padj_f1_kallisto <- c()
Venn_idio_F4_dammit$foldChange_f4_kallisto <- c()
Venn_idio_F4_dammit$padj_f4_kallisto <- c()


j <- 0
for (i in row.names(Venn_idio_F4_dammit)){
  j <- j +1 
  Venn_idio_F4_dammit$foldChange_idio_kallisto[j] <- idio_kallisto_foldchange[i]
  Venn_idio_F4_dammit$padj_idio_kallisto[j] <- idio_kallisto_padj[i]
  Venn_idio_F4_dammit$foldChange_f1_kallisto[j] <- f1out_kallisto_foldchange[i]
  Venn_idio_F4_dammit$padj_f1_kallisto[j] <- f1out_kallisto_padj[i]
  Venn_idio_F4_dammit$foldChange_f4_kallisto[j] <- f4out_kallisto_foldchange[i]
  Venn_idio_F4_dammit$padj_f4_kallisto[j] <- f4out_kallisto_padj[i]
  
}


#write the output table
write.xlsx(Venn_idio_F4_dammit, file = paste(results_dir_2, "dammit_annotation_idio_F4.xlsx", sep = "/"), showNA = T)


################################
### idio vs F1 intersection ####
################################

###create output dir  
results_dir_2 <- paste(results_dir, "idio_F1_intersect", sep = "/")
if (!(file.exists(results_dir_2))) {
  dir.create(results_dir_2)
}


#### get annotation for the upregulated genes ####
Venn_idio_F1_genes_kallisto <- read.delim("intersection_tables/Venn_Idio_F1Out_UP_list.csv", row.names = 1, sep = ",")


Venn_idio_F1_dammit <- subset(dammit_annot, row.names(dammit_annot) %in% Venn_idio_F1_genes_kallisto$x)


#Add the logfold change and p.value to the table 
Venn_idio_F1_dammit$foldChange_idio_kallisto <- c()
Venn_idio_F1_dammit$padj_idio_kallisto <- c()
Venn_idio_F1_dammit$foldChange_f1_kallisto <- c()
Venn_idio_F1_dammit$padj_f1_kallisto <- c()
Venn_idio_F1_dammit$foldChange_f4_kallisto <- c()
Venn_idio_F1_dammit$padj_f4_kallisto <- c()


j <- 0
for (i in row.names(Venn_idio_F1_dammit)){
  j <- j + 1
  Venn_idio_F1_dammit$foldChange_idio_kallisto[j] <- idio_kallisto_foldchange[i]
  Venn_idio_F1_dammit$padj_idio_kallisto[j] <- idio_kallisto_padj[i]
  Venn_idio_F1_dammit$foldChange_f1_kallisto[j] <- f1out_kallisto_foldchange[i]
  Venn_idio_F1_dammit$padj_f1_kallisto[j] <- f1out_kallisto_padj[i]
  Venn_idio_F1_dammit$foldChange_f4_kallisto[j] <- f4out_kallisto_foldchange[i]
  Venn_idio_F1_dammit$padj_f4_kallisto[j] <- f4out_kallisto_padj[i]
  
}


#write the output table
write.xlsx(Venn_idio_F1_dammit, file = paste(results_dir_2, "dammit_annotation_idio_F1.xlsx", sep = "/"), showNA = T)


