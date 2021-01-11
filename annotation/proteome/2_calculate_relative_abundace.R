#Script used to calculate the relative abundance of each protein

#misc 
setwd("C:/Users/Faculdade/Desktop/Dissertação/results 2020/mapping_to_genome/annotation/02_proteome_annotation")
rm(list = ls())

#load libraries
library(VennDiagram)


#load the minimap output table 
minimap_out.final <- read.delim(file = "minimap_final.tsv")

gene_vector <- minimap_out.final$gene_id
names(gene_vector) <- minimap_out.final$Query
annotation_vector <- minimap_out.final$cathacyc
names(annotation_vector) <- minimap_out.final$Query

#Load samples files  ---- 
vac_1_table <- read.delim("protein_data_per_sample/proteome_table_Vac#1_not_treated_v2.tsv")
vac_2_table <- read.delim("protein_data_per_sample/proteome_table_Vac#2_not_treated_v2.tsv")
vac_1_table <- vac_1_table[,-13]
vac_2_table <- vac_2_table[,-13]

ton_1_table <- read.delim("protein_data_per_sample/proteome_table_Ton#1_not_treated_v2.tsv")
ton_2_table <- read.delim("protein_data_per_sample/proteome_table_Ton#2_not_treated_v2.tsv")
ton_1_table <- ton_1_table[,-13]
ton_2_table <- ton_2_table[,-13]

###For vacule 1 ----
##calculate relative abundance 
#replace "None" values per 0 (area not calculated)
vac_1_table[c(vac_1_table$Area == "None"), 9] <- 0 

#Compute relative abundance
total_area_vac1 <- sum(as.numeric(vac_1_table$Area))
vac_1_table$relative_abundance <- (as.numeric(vac_1_table$Area) / total_area_vac1)


#associate with the stringtie gene id 
vac_1_table$stringTie_id <- gene_vector[as.character(vac_1_table$Accession)]
vac_1_table$cathacyc_annot <- annotation_vector[as.character(vac_1_table$Accession)]

vac_1_table <- vac_1_table[,c(14,1,15, 3:13)]
colnames(vac_1_table) <- c("Gene Id", "Caros acession", "cathacyc annotation", "score", "coverage", "proteins", "Unique peptides", "Peptides", "PSM", "Area", "AA", "kDA", "pI", "relative abundance")

write.table(vac_1_table, file = "protein_data_per_sample/vac_1_table.tsv", sep = "\t", quote = F, row.names = F)



###For vacule 2 ----
##calculate relative abundance 
#replace "None" values per 0 (area not calculated)
vac_2_table[c(vac_2_table$Area == "None"), 9] <- 0 

#Compute relative abundance
total_area_vac2 <- sum(as.numeric(vac_2_table$Area))
vac_2_table$relative_abundance <- (as.numeric(vac_2_table$Area) / total_area_vac2)


#associate with the stringtie gene id 
vac_2_table$stringTie_id <- gene_vector[as.character(vac_2_table$Accession)]
vac_2_table$cathacyc_annot <- annotation_vector[as.character(vac_2_table$Accession)]

vac_2_table <- vac_2_table[,c(14,1,15, 3:13)]
colnames(vac_2_table) <- c("Gene Id", "Caros acession", "cathacyc annotation", "score", "coverage", "proteins", "Unique peptides", "Peptides", "PSM", "Area", "AA", "kDA", "pI", "relative abundance")

write.table(vac_2_table, file = "protein_data_per_sample/vac_2_table.tsv", sep = "\t", quote = F, row.names = F)



###For tonoplast 1 ----
##calculate relative abundance 
#replace "None" values per 0 (area not calculated)
ton_1_table[c(ton_1_table$Area == "None"), 9] <- 0 

#Compute relative abundance
total_area_ton1 <- sum(as.numeric(ton_1_table$Area))
ton_1_table$relative_abundance <- (as.numeric(ton_1_table$Area) / total_area_ton1)


#associate with the stringtie gene id 
ton_1_table$stringTie_id <- gene_vector[as.character(ton_1_table$Accession)]
ton_1_table$cathacyc_annot <- annotation_vector[as.character(ton_1_table$Accession)]

ton_1_table <- ton_1_table[,c(14,1,15, 3:13)]
colnames(ton_1_table) <- c("Gene Id", "Caros acession", "cathacyc annotation", "score", "coverage", "proteins", "Unique peptides", "Peptides", "PSM", "Area", "AA", "kDA", "pI", "relative abundance")

write.table(ton_1_table, file = "protein_data_per_sample/ton_1_table.tsv", sep = "\t", quote = F, row.names = F)



###For tonoplast 2 ----
##calculate relative abundance 
#replace "None" values per 0 (area not calculated)
ton_2_table[c(ton_2_table$Area == "None"), 9] <- 0 

#Compute relative abundance
total_area_ton1 <- sum(as.numeric(ton_2_table$Area))
ton_2_table$relative_abundance <- (as.numeric(ton_2_table$Area) / total_area_ton1)

#associate with the stringtie gene id 
ton_2_table$stringTie_id <- gene_vector[as.character(ton_2_table$Accession)]
ton_2_table$cathacyc_annot <- annotation_vector[as.character(ton_2_table$Accession)]

ton_2_table <- ton_2_table[,c(14,1,15, 3:13)]
colnames(ton_2_table) <- c("Gene Id", "Caros acession", "cathacyc annotation", "score", "coverage", "proteins", "Unique peptides", "Peptides", "PSM", "Area", "AA", "kDA", "pI", "relative abundance")

write.table(ton_2_table, file = "protein_data_per_sample/ton_2_table.tsv", sep = "\t", quote = F, row.names = F)



###make veen diagram
myCol <- c("#000000", "#E69F00", "#56B4E9", "#009E73")


venn.diagram(x = list(vac_1_table$`Caros acession`, vac_2_table$`Caros acession`, ton_1_table$`Caros acession`, ton_2_table$`Caros acession`), 
             category.names = c("vac#1", "Vac#2", "Ton#1", "Ton#2"),
             filename = "Venn_diagram.png", 
             output =T, fill = myCol)

