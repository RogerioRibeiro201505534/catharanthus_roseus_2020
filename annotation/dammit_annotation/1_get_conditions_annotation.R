#use this script to treat individual DESeq2_tables 

#Set working directory 
setwd("C:/Users/Faculdade/Desktop/Dissertação/results 2020/mapping_to_genome/annotation/03_dammit")
rm(list = ls())
options(java.parameters = "-Xmx2048m")  

#load libraries 
library("xlsx")


##output dir  
results_dir <- "./results"
if (!(file.exists(results_dir))) {
  dir.create(results_dir)
}


#Load final dammit annotation
dammit_annot <- read.delim("StepsToFinalAnnotation/dammit_annotation_final.tsv", sep = "\t",  header= TRUE, row.names = 1, stringsAsFactors = FALSE)
swissprot_annot <- dammit_annot$Swissprot.name
names(swissprot_annot) <- row.names(dammit_annot)

#Load individual annotations 
idio_kallisto <- read.delim(file = "DESeq2_tables/meso_vs_idio_kallisto_Res.csv", sep = ",")
f1_kallisto <- read.delim(file = "DESeq2_tables/F1in_vs_F1out_kallisto_Res.csv", sep = ",")
f4_kallisto <- read.delim(file = "DESeq2_tables/f4in_vs_f4out_kallisto_Res.csv", sep = ",")

#Load expression levels and build a average_data_set
kallisto_expression_levels <- read.delim(file = "DeSeq2_tables/kallisto_normalizedCountsTable.txt")
kallisto_expression_levels <- data.frame("avgF1in" = rowMeans(kallisto_expression_levels[1:3]), "avgF1out" = rowMeans(kallisto_expression_levels[4:6]), "avgF4in" = rowMeans(kallisto_expression_levels[7:9]), "avgF4out" = rowMeans(kallisto_expression_levels[10:12]), "avgIdio" = rowMeans(kallisto_expression_levels[13:15]), "avgMeso" = rowMeans(kallisto_expression_levels[18:20]))
kallisto_expression_levels <- kallisto_expression_levels[,c(1,2,3,4,6,5)]


###For meso Vs idio annotation---- 
results_dir_idio <- paste(results_dir, "meso_vs_idio", sep = "/")
if (!(file.exists(results_dir_idio))) {
  dir.create(results_dir_idio)
}

#Add annotation
idio_kallisto$gene_name <- c()
idio_kallisto$proteome_location <- c()
idio_kallisto$cathacyc_annotation <- c()

for (i in 1:length(idio_kallisto$X)){
  idio_kallisto$gene_name[i] <- as.character(swissprot_annot[as.character(idio_kallisto$X[i])])
  idio_kallisto$proteome_location[i] <- as.character(dammit_annot[as.character(idio_kallisto$X[i]),15])
  idio_kallisto$cathacyc_annotation[i] <- as.character(dammit_annot[as.character(idio_kallisto$X[i]),16])
}

idio_kallisto_ordered <- idio_kallisto[order(abs(idio_kallisto$log2FoldChange), decreasing = T),]
row.names(idio_kallisto_ordered) <- idio_kallisto_ordered$X


#Get upregulated and downregulated genes row names to subset 
idio_kallisto_sig_up <- read.delim("DESeq2_tables/meso_vs_idio_kallisto_sigUP.csv", sep = ",")[,2]
idio_kallisto_sig_down <- read.delim("DESeq2_tables/meso_vs_idio_kallisto_sigDown.csv", sep = ",")[,2]

#Get tables for the idio up and idio down genes 
idio_kallisto_sig_up_annotated <- subset(idio_kallisto_ordered, idio_kallisto_ordered$X %in% idio_kallisto_sig_up)
idio_kallisto_sig_down_annotated <- subset(idio_kallisto_ordered, idio_kallisto_ordered$X %in% idio_kallisto_sig_down)

#Write the resulting tables in xlsx files 
write.xlsx(idio_kallisto_ordered, file = paste(results_dir_idio, "meso_vs_idio_annotated.xlsx", sep = "/"), sheetName = "all_ordered", row.names = F)
gc()
write.xlsx(idio_kallisto_sig_up_annotated, file = paste(results_dir_idio, "meso_vs_idio_annotated.xlsx", sep = "/"), sheetName = "up_ordered", append = T, row.names = F)
gc()
write.xlsx(idio_kallisto_sig_down_annotated, file = paste(results_dir_idio, "meso_vs_idio_annotated.xlsx", sep = "/"), sheetName = "down_ordered", append = T, row.names = F)
gc()

#Make a table for upregulated genes
#Row names, avg expression in 6 conditions, foldChange for each, log2foldchange, 2 collumssaying upregulated, geneName, 
f1_kallisto_sig_up <- read.delim("DESeq2_tables/f1in_vs_f1out_kallisto_sigUP.csv", sep = ",")[,2]
f4_kallisto_sig_up <- read.delim("DESeq2_tables/f4in_vs_f4out_kallisto_sigUP.csv", sep = ",")[,2]


idio_kallisto_sig_up_v2 <- kallisto_expression_levels[as.character(idio_kallisto_sig_up),]
idio_kallisto_sig_up_v2$Foldchange_meso_idio <- round((idio_kallisto_sig_up_v2$avgIdio / idio_kallisto_sig_up_v2$avgMeso), 2)
idio_kallisto_sig_up_v2$Foldchange_f1in_f1out <- round((idio_kallisto_sig_up_v2$avgF1out / idio_kallisto_sig_up_v2$avgF1in), 2)
idio_kallisto_sig_up_v2$Foldchange_f4in_f4out <- round((idio_kallisto_sig_up_v2$avgF4out / idio_kallisto_sig_up_v2$avgF4in), 2) 
idio_kallisto_sig_up_v2$sigF1 <- ifelse(rownames(idio_kallisto_sig_up_v2) %in% f1_kallisto_sig_up, "YES", NA)
idio_kallisto_sig_up_v2$sigF4 <- ifelse(rownames(idio_kallisto_sig_up_v2) %in% f4_kallisto_sig_up, "YES", NA)
idio_kallisto_sig_up_v2$gene_name <- dammit_annot[as.character(row.names(idio_kallisto_sig_up_v2)), 3]
idio_kallisto_sig_up_v2$proteome <- dammit_annot[as.character(row.names(idio_kallisto_sig_up_v2)), 15]

#add goslim terms
goslim_annotation <- read.delim(file = "StepsToFinalAnnotation/dammit_annotation_goslim.tsv", row.names = 1)[,c(10,12,14)]
idio_kallisto_sig_up_v2$goslim.BP <- goslim_annotation[as.character(row.names(idio_kallisto_sig_up_v2)), 1]
idio_kallisto_sig_up_v2$goslim.MF <- goslim_annotation[as.character(row.names(idio_kallisto_sig_up_v2)), 2]
idio_kallisto_sig_up_v2$goslim.CC <- goslim_annotation[as.character(row.names(idio_kallisto_sig_up_v2)), 3]

##classification of each protein in transporter, enzyme, TF, other or unkown can be made here
idio_kallisto_sig_up_v2$category <- rep(NA, length(idio_kallisto_sig_up_v2$avgF1in))


for (i in 1:length(idio_kallisto_sig_up_v2$avgF1in)){
  if (grepl(";", idio_kallisto_sig_up_v2$gene_name[i], fixed = T)){
    idio_kallisto_sig_up_v2$category[i] <- "ambigous"
  }
  if (grepl("secondary metabolic process", idio_kallisto_sig_up_v2$goslim.BP[i], fixed = T) && is.na(idio_kallisto_sig_up_v2$category[i])){
    idio_kallisto_sig_up_v2$category[i] <- "secondary metabolic process"
  }
  if (grepl("transporter", idio_kallisto_sig_up_v2$goslim.MF[i], fixed = T) && is.na(idio_kallisto_sig_up_v2$category[i])){
    idio_kallisto_sig_up_v2$category[i] <- "Transporter"
  }
  if (grepl("EC", idio_kallisto_sig_up_v2$gene_name[i], fixed = T) && is.na(idio_kallisto_sig_up_v2$category[i])){
    idio_kallisto_sig_up_v2$category[i] <- "Enzyme"
  }
  if (idio_kallisto_sig_up_v2$gene_name[i] == "" && is.na(idio_kallisto_sig_up_v2$category[i])){
    idio_kallisto_sig_up_v2$category[i] <- "Unknown"
  }
  if (is.na(idio_kallisto_sig_up_v2$category[i])){
    idio_kallisto_sig_up_v2$category[i] <- "Other"
  }
}

#solve some ambigous. If the two ambigous are enzymes annotation can be made with "enzyme"
for (i in 1:length(idio_kallisto_sig_up_v2$avgF1in)){
  if (idio_kallisto_sig_up_v2$category[i] == "ambigous"){
    individual_annotation <- strsplit(idio_kallisto_sig_up_v2$gene_name[i], ";")[[1]]
    if (all(grepl("EC", individual_annotation, fixed =T))){
      idio_kallisto_sig_up_v2$category[i] <- "Enzyme"
    }
  }
}

#flag possible transcription factors in the annotation 
for (i in 1:length(idio_kallisto_sig_up_v2$avgF1in)){
  if (grepl("transcription factor", idio_kallisto_sig_up_v2$gene_name[i], fixed = T) | grepl("DNA-binding transcription factor activity", idio_kallisto_sig_up_v2$goslim.MF[i], fixed = TRUE)){
    idio_kallisto_sig_up_v2$category[i] <- paste("(possible TF)", idio_kallisto_sig_up_v2$category[i], sep = " ")
  }
}


#add the realtive protein expression in each sample 
proteome_abundances <- read.xlsx(file = "Proteome_final_tables/proteome_sup_tableV2.xlsx", sheetIndex = 1)


abundances_vac_1 <- proteome_abundances$vac_1_relative_abundance
names(abundances_vac_1) <- proteome_abundances$Caros.acession

abundances_vac_2 <- proteome_abundances$vac_2_relative_abundance
names(abundances_vac_2) <- proteome_abundances$Caros.acession

abundances_ton_1 <- proteome_abundances$ton_1_relative_abundance
names(abundances_ton_1) <- proteome_abundances$Caros.acession

abundances_ton_2 <- proteome_abundances$ton_2_relative_abundance
names(abundances_ton_2) <- proteome_abundances$Caros.acession

abundances_vac_mean <- proteome_abundances$vac_mean
names(abundances_vac_mean) <- proteome_abundances$Caros.acession

abundances_ton_mean <- proteome_abundances$ton_mean
names(abundances_ton_mean) <- proteome_abundances$Caros.acession

for (i in 1:length(idio_kallisto_sig_up_v2$proteome)){
  if (!is.na(idio_kallisto_sig_up_v2$proteome[i])){
    caros <- strsplit(idio_kallisto_sig_up_v2$proteome[i], "; ", fixed = T)[[1]]
    for (j in 1:length(caros)){
      caro_acession_numer <- strsplit(caros[j], " ")[[1]][1]
      if (j == 1){
        vac_1 <- abundances_vac_1[as.character(caro_acession_numer)]
        vac_2 <- abundances_vac_2[as.character(caro_acession_numer)]
        ton_1 <- abundances_ton_1[as.character(caro_acession_numer)]
        ton_2 <- abundances_ton_2[as.character(caro_acession_numer)]
        vac_mean <- abundances_vac_mean[as.character(caro_acession_numer)]
        ton_mean <- abundances_ton_mean[as.character(caro_acession_numer)]
      }
      if (j >= 2){
        vac_1 <- paste(vac_1, abundances_vac_1[as.character(caro_acession_numer)], sep = ";")
        vac_2 <- paste(vac_2, abundances_vac_2[as.character(caro_acession_numer)], sep = ";")
        ton_1 <- paste(ton_1, abundances_ton_1[as.character(caro_acession_numer)], sep = ";")
        ton_2 <- paste(ton_2, abundances_ton_2[as.character(caro_acession_numer)], sep = ";")
        vac_mean <- paste(vac_mean, abundances_vac_mean[as.character(caro_acession_numer)], sep = ";")
        ton_mean <- paste(ton_mean, abundances_ton_mean[as.character(caro_acession_numer)], sep = ";")
      }
    }
    idio_kallisto_sig_up_v2$abudance_vac1[i] <- vac_1
    idio_kallisto_sig_up_v2$abudance_vac2[i] <- vac_2
    idio_kallisto_sig_up_v2$abudance_ton1[i] <- ton_1
    idio_kallisto_sig_up_v2$abudance_ton2[i] <- ton_2
    idio_kallisto_sig_up_v2$abudance_vac_mean[i] <- vac_mean
    idio_kallisto_sig_up_v2$abudance_ton_mean[i] <- ton_mean
  }
}

idio_kallisto_sig_up_v2 <- idio_kallisto_sig_up_v2[,c(1:12,17,14:16,13,18:23)]
write.xlsx(idio_kallisto_sig_up_v2, file = paste(results_dir_idio, "meso_vs_idio_annotated.xlsx", sep = "/"), sheetName = "Up.V2._2raw", append = T)

###For f1in Vs f1out annotation---- 
results_dir_f1 <- paste(results_dir, "f1in_vs_f1out", sep = "/")
if (!(file.exists(results_dir_f1))) {
  dir.create(results_dir_f1)
}

#Add annotation
f1_kallisto$gene_name <- c()
f1_kallisto$proteome_location <- c()
f1_kallisto$cathacyc_annotation <- c()

for (i in 1:length(f1_kallisto$X)){
  f1_kallisto$gene_name[i] <- as.character(swissprot_annot[as.character(f1_kallisto$X[i])])
  f1_kallisto$proteome_location[i] <- as.character(dammit_annot[as.character(f1_kallisto$X[i]),15])
  f1_kallisto$cathacyc_annotation[i] <- as.character(dammit_annot[as.character(f1_kallisto$X[i]),16])
}

f1_kallisto_ordered <- f1_kallisto[order(abs(f1_kallisto$log2FoldChange), decreasing = T),]
row.names(f1_kallisto_ordered) <- f1_kallisto_ordered$X


#Get upregulated and downregulated genes row names to subset 
f1_kallisto_sig_up <- read.delim("DESeq2_tables/F1in_vs_F1out_kallisto_sigUP.csv", sep = ",")[,2]
f1_kallisto_sig_down <- read.delim("DESeq2_tables/F1in_vs_F1out_kallisto_sigDOWN.csv", sep = ",")[,2]

#Get tables for the idio up and idio down genes 
f1_kallisto_sig_up_annotated <- subset(f1_kallisto_ordered, f1_kallisto_ordered$X %in% f1_kallisto_sig_up)
f1_kallisto_sig_down_annotated <- subset(f1_kallisto_ordered, f1_kallisto_ordered$X %in% f1_kallisto_sig_down)

#Write the resulting tables in xlsx files 
write.xlsx(f1_kallisto_ordered, file = paste(results_dir_f1, "f1in_vs_f1out_annotated.xlsx", sep = "/"), sheetName = "all_ordered", row.names = F)
gc()
write.xlsx(f1_kallisto_sig_up_annotated, file = paste(results_dir_f1, "f1in_vs_f1out_annotated.xlsx", sep = "/"), sheetName = "up_ordered", append = T, row.names = F)
gc()
write.xlsx(f1_kallisto_sig_down_annotated, file = paste(results_dir_f1, "f1in_vs_f1out_annotated.xlsx", sep = "/"), sheetName = "down_ordered", append = T, row.names = F)
gc()

#Make a table for upregulated genes
#Row names, avg expression in 6 conditions, foldChange for each, log2foldchange, 2 collumssaying upregulated, geneName, 
idio_kallisto_sig_up <- read.delim("DESeq2_tables/meso_vs_idio_kallisto_sigUP.csv", sep = ",")[,2]
f4_kallisto_sig_up <- read.delim("DESeq2_tables/f4in_vs_f4out_kallisto_sigUP.csv", sep = ",")[,2]


f1_kallisto_sig_up_v2 <- kallisto_expression_levels[as.character(f1_kallisto_sig_up),]
f1_kallisto_sig_up_v2$Foldchange_meso_idio <- round((f1_kallisto_sig_up_v2$avgIdio / f1_kallisto_sig_up_v2$avgMeso), 2)
f1_kallisto_sig_up_v2$Foldchange_f1in_f1out <- round((f1_kallisto_sig_up_v2$avgF1out / f1_kallisto_sig_up_v2$avgF1in), 2)
f1_kallisto_sig_up_v2$Foldchange_f4in_f4out <- round((f1_kallisto_sig_up_v2$avgF4out / f1_kallisto_sig_up_v2$avgF4in), 2) 
f1_kallisto_sig_up_v2$idio <- ifelse(rownames(f1_kallisto_sig_up_v2) %in% idio_kallisto_sig_up, "YES", NA)
f1_kallisto_sig_up_v2$sigF4 <- ifelse(rownames(f1_kallisto_sig_up_v2) %in% f4_kallisto_sig_up, "YES", NA)
f1_kallisto_sig_up_v2$gene_name <- dammit_annot[as.character(row.names(f1_kallisto_sig_up_v2)), 3]
f1_kallisto_sig_up_v2$proteome <- dammit_annot[as.character(row.names(f1_kallisto_sig_up_v2)), 15]

#add goslim terms
goslim_annotation <- read.delim(file = "StepsToFinalAnnotation/dammit_annotation_goslim.tsv", row.names = 1)[,c(10,12,14)]
f1_kallisto_sig_up_v2$goslim.BP <- goslim_annotation[as.character(row.names(f1_kallisto_sig_up_v2)), 1]
f1_kallisto_sig_up_v2$goslim.MF <- goslim_annotation[as.character(row.names(f1_kallisto_sig_up_v2)), 2]
f1_kallisto_sig_up_v2$goslim.CC <- goslim_annotation[as.character(row.names(f1_kallisto_sig_up_v2)), 3]

##classification of each protein in transporter, enzyme, TF, other or unkown can be made here
f1_kallisto_sig_up_v2$category <- rep(NA, length(f1_kallisto_sig_up_v2$gene_name))


for (i in 1:length(f1_kallisto_sig_up_v2$gene_name)){
  if (grepl(";", f1_kallisto_sig_up_v2$gene_name[i], fixed = T)){
    f1_kallisto_sig_up_v2$category[i] <- "ambigous"
  }
  if (grepl("secondary metabolic process", f1_kallisto_sig_up_v2$goslim.BP[i], fixed = T) && is.na(f1_kallisto_sig_up_v2$category[i])){
    f1_kallisto_sig_up_v2$category[i] <- "secondary metabolic process"
  }
  if (grepl("transporter", f1_kallisto_sig_up_v2$goslim.MF[i], fixed = T) && is.na(f1_kallisto_sig_up_v2$category[i])){
    f1_kallisto_sig_up_v2$category[i] <- "Transporter"
  }
  if (grepl("EC", f1_kallisto_sig_up_v2$gene_name[i], fixed = T) && is.na(f1_kallisto_sig_up_v2$category[i])){
    f1_kallisto_sig_up_v2$category[i] <- "Enzyme"
  }
  if (f1_kallisto_sig_up_v2$gene_name[i] == "" && is.na(f1_kallisto_sig_up_v2$category[i])){
    f1_kallisto_sig_up_v2$category[i] <- "Unknown"
  }
  if (is.na(f1_kallisto_sig_up_v2$category[i])){
    f1_kallisto_sig_up_v2$category[i] <- "Other"
  }
}

#solve some ambigous. If the two ambigous are enzymes annotation can be made with "enzyme"
for (i in 1:length(f1_kallisto_sig_up_v2$gene_name)){
  if (f1_kallisto_sig_up_v2$category[i] == "ambigous"){
    individual_annotation <- strsplit(f1_kallisto_sig_up_v2$gene_name[i], ";")[[1]]
    if (all(grepl("EC", individual_annotation, fixed =T))){
      f1_kallisto_sig_up_v2$category[i] <- "Enzyme"
    }
  }
}

#flag possible transcription factors in the annotation 
for (i in 1:length(f1_kallisto_sig_up_v2$gene_name)){
  if (grepl("transcription factor", f1_kallisto_sig_up_v2$gene_name[i], fixed = T) | grepl("DNA-binding transcription factor activity", f1_kallisto_sig_up_v2$goslim.MF[i], fixed = TRUE)){
    f1_kallisto_sig_up_v2$category[i] <- paste("(possible TF)", f1_kallisto_sig_up_v2$category[i], sep = " ")
  }
}


#add the realtive protein expression in each sample 
proteome_abundances <- read.xlsx(file = "Proteome_final_tables/proteome_sup_tableV2.xlsx", sheetIndex = 1)


abundances_vac_1 <- proteome_abundances$vac_1_relative_abundance
names(abundances_vac_1) <- proteome_abundances$Caros.acession

abundances_vac_2 <- proteome_abundances$vac_2_relative_abundance
names(abundances_vac_2) <- proteome_abundances$Caros.acession

abundances_ton_1 <- proteome_abundances$ton_1_relative_abundance
names(abundances_ton_1) <- proteome_abundances$Caros.acession

abundances_ton_2 <- proteome_abundances$ton_2_relative_abundance
names(abundances_ton_2) <- proteome_abundances$Caros.acession

abundances_vac_mean <- proteome_abundances$vac_mean
names(abundances_vac_mean) <- proteome_abundances$Caros.acession

abundances_ton_mean <- proteome_abundances$ton_mean
names(abundances_ton_mean) <- proteome_abundances$Caros.acession

for (i in 1:length(f1_kallisto_sig_up_v2$proteome)){
  if (!is.na(f1_kallisto_sig_up_v2$proteome[i])){
    caros <- strsplit(f1_kallisto_sig_up_v2$proteome[i], "; ", fixed = T)[[1]]
    for (j in 1:length(caros)){
      caro_acession_numer <- strsplit(caros[j], " ")[[1]][1]
      if (j == 1){
        vac_1 <- abundances_vac_1[as.character(caro_acession_numer)]
        vac_2 <- abundances_vac_2[as.character(caro_acession_numer)]
        ton_1 <- abundances_ton_1[as.character(caro_acession_numer)]
        ton_2 <- abundances_ton_2[as.character(caro_acession_numer)]
        vac_mean <- abundances_vac_mean[as.character(caro_acession_numer)]
        ton_mean <- abundances_ton_mean[as.character(caro_acession_numer)]
      }
      if (j >= 2){
        vac_1 <- paste(vac_1, abundances_vac_1[as.character(caro_acession_numer)], sep = ";")
        vac_2 <- paste(vac_2, abundances_vac_2[as.character(caro_acession_numer)], sep = ";")
        ton_1 <- paste(ton_1, abundances_ton_1[as.character(caro_acession_numer)], sep = ";")
        ton_2 <- paste(ton_2, abundances_ton_2[as.character(caro_acession_numer)], sep = ";")
        vac_mean <- paste(vac_mean, abundances_vac_mean[as.character(caro_acession_numer)], sep = ";")
        ton_mean <- paste(ton_mean, abundances_ton_mean[as.character(caro_acession_numer)], sep = ";")
      }
    }
    f1_kallisto_sig_up_v2$abudance_vac1[i] <- vac_1
    f1_kallisto_sig_up_v2$abudance_vac2[i] <- vac_2
    f1_kallisto_sig_up_v2$abudance_ton1[i] <- ton_1
    f1_kallisto_sig_up_v2$abudance_ton2[i] <- ton_2
    f1_kallisto_sig_up_v2$abudance_vac_mean[i] <- vac_mean
    f1_kallisto_sig_up_v2$abudance_ton_mean[i] <- ton_mean
  }
}

f1_kallisto_sig_up_v2 <- f1_kallisto_sig_up_v2[,c(1:12,17,14:16,13,18:23)]
write.xlsx(f1_kallisto_sig_up_v2, file = paste(results_dir_f1, "f1in_vs_f1out_annotated.xlsx", sep = "/"), sheetName = "Up.V2.raw", append = T)


###For f4in Vs f4out annotation---- 

results_dir_f4 <- paste(results_dir, "f4in_vs_f4out", sep = "/")
if (!(file.exists(results_dir_f4))) {
  dir.create(results_dir_f4)
}


#Add annotation
f4_kallisto$gene_name <- c()
f4_kallisto$proteome_location <- c()
f4_kallisto$cathacyc_annotation <- c()

for (i in 1:length(f4_kallisto$X)){
  f4_kallisto$gene_name[i] <- as.character(swissprot_annot[as.character(f4_kallisto$X[i])])
  f4_kallisto$proteome_location[i] <- as.character(dammit_annot[as.character(f4_kallisto$X[i]),15])
  f4_kallisto$cathacyc_annotation[i] <- as.character(dammit_annot[as.character(f4_kallisto$X[i]),16])
}

f4_kallisto_ordered <- f4_kallisto[order(abs(f4_kallisto$log2FoldChange), decreasing = T),]
row.names(f4_kallisto_ordered) <- f4_kallisto_ordered$X


#Get upregulated and downregulated genes row names to subset 
f4_kallisto_sig_up <- read.delim("DESeq2_tables/f4in_vs_f4out_kallisto_sigUP.csv", sep = ",")[,2]
f4_kallisto_sig_down <- read.delim("DESeq2_tables/f4in_vs_f4out_kallisto_sigDOWN.csv", sep = ",")[,2]

#Get tables for the idio up and idio down genes 
f4_kallisto_sig_up_annotated <- subset(f4_kallisto_ordered, f4_kallisto_ordered$X %in% f4_kallisto_sig_up)
f4_kallisto_sig_down_annotated <- subset(f4_kallisto_ordered, f4_kallisto_ordered$X %in% f4_kallisto_sig_down)

#Write the resulting tables in xlsx files 
write.xlsx(f4_kallisto_ordered, file = paste(results_dir_f4, "f4in_vs_f4out_annotated.xlsx", sep = "/"), sheetName = "all_ordered", row.names = F)
gc()
write.xlsx(f4_kallisto_sig_up_annotated, file = paste(results_dir_f4, "f4in_vs_f4out_annotated.xlsx", sep = "/"), sheetName = "up_ordered", append = T, row.names = F)
gc()
write.xlsx(f4_kallisto_sig_down_annotated, file = paste(results_dir_f4, "f4in_vs_f4out_annotated.xlsx", sep = "/"), sheetName = "down_ordered", append = T, row.names = F)
gc()

#Make a table for upregulated genes
#Row names, avg expression in 6 conditions, foldChange for each, log2foldchange, 2 collumssaying upregulated, geneName, 
idio_kallisto_sig_up <- read.delim("DESeq2_tables/meso_vs_idio_kallisto_sigUP.csv", sep = ",")[,2]
f4_kallisto_sig_up <- read.delim("DESeq2_tables/f4in_vs_f4out_kallisto_sigUP.csv", sep = ",")[,2]


f4_kallisto_sig_up_v2 <- kallisto_expression_levels[as.character(f4_kallisto_sig_up),]
f4_kallisto_sig_up_v2$Foldchange_meso_idio <- round((f4_kallisto_sig_up_v2$avgIdio / f4_kallisto_sig_up_v2$avgMeso), 2)
f4_kallisto_sig_up_v2$Foldchange_f1in_f1out <- round((f4_kallisto_sig_up_v2$avgF1out / f4_kallisto_sig_up_v2$avgF1in), 2)
f4_kallisto_sig_up_v2$Foldchange_f4in_f4out <- round((f4_kallisto_sig_up_v2$avgF4out / f4_kallisto_sig_up_v2$avgF4in), 2) 
f4_kallisto_sig_up_v2$sigidio <- ifelse(rownames(f4_kallisto_sig_up_v2) %in% idio_kallisto_sig_up, "YES", NA)
f4_kallisto_sig_up_v2$sigF1 <- ifelse(rownames(f4_kallisto_sig_up_v2) %in% f1_kallisto_sig_up, "YES", NA)
f4_kallisto_sig_up_v2$gene_name <- dammit_annot[as.character(row.names(f4_kallisto_sig_up_v2)), 3]
f4_kallisto_sig_up_v2$proteome <- dammit_annot[as.character(row.names(f4_kallisto_sig_up_v2)), 15]

#add goslim terms
goslim_annotation <- read.delim(file = "StepsToFinalAnnotation/dammit_annotation_goslim.tsv", row.names = 1)[,c(10,12,14)]
f4_kallisto_sig_up_v2$goslim.BP <- goslim_annotation[as.character(row.names(f4_kallisto_sig_up_v2)), 1]
f4_kallisto_sig_up_v2$goslim.MF <- goslim_annotation[as.character(row.names(f4_kallisto_sig_up_v2)), 2]
f4_kallisto_sig_up_v2$goslim.CC <- goslim_annotation[as.character(row.names(f4_kallisto_sig_up_v2)), 3]

##classification of each protein in transporter, enzyme, TF, other or unkown can be made here
f4_kallisto_sig_up_v2$category <- rep(NA, length(f4_kallisto_sig_up_v2$gene_name))


for (i in 1:length(f4_kallisto_sig_up_v2$gene_name)){
  if (grepl(";", f4_kallisto_sig_up_v2$gene_name[i], fixed = T)){
    f4_kallisto_sig_up_v2$category[i] <- "ambigous"
  }
  if (grepl("secondary metabolic process", f4_kallisto_sig_up_v2$goslim.BP[i], fixed = T) && is.na(f4_kallisto_sig_up_v2$category[i])){
    f4_kallisto_sig_up_v2$category[i] <- "secondary metabolic process"
  }
  if (grepl("transporter", f4_kallisto_sig_up_v2$goslim.MF[i], fixed = T) && is.na(f4_kallisto_sig_up_v2$category[i])){
    f4_kallisto_sig_up_v2$category[i] <- "Transporter"
  }
  if (grepl("EC", f4_kallisto_sig_up_v2$gene_name[i], fixed = T) && is.na(f4_kallisto_sig_up_v2$category[i])){
    f4_kallisto_sig_up_v2$category[i] <- "Enzyme"
  }
  if (f4_kallisto_sig_up_v2$gene_name[i] == "" && is.na(f4_kallisto_sig_up_v2$category[i])){
    f4_kallisto_sig_up_v2$category[i] <- "Unknown"
  }
  if (is.na(f4_kallisto_sig_up_v2$category[i])){
    f4_kallisto_sig_up_v2$category[i] <- "Other"
  }
}

#solve some ambigous. If the two ambigous are enzymes annotation can be made with "enzyme"
for (i in 1:length(f4_kallisto_sig_up_v2$gene_name)){
  if (f4_kallisto_sig_up_v2$category[i] == "ambigous"){
    individual_annotation <- strsplit(f4_kallisto_sig_up_v2$gene_name[i], ";")[[1]]
    if (all(grepl("EC", individual_annotation, fixed =T))){
      f4_kallisto_sig_up_v2$category[i] <- "Enzyme"
    }
  }
}

#flag possible transcription factors in the annotation 
for (i in 1:length(f4_kallisto_sig_up_v2$gene_name)){
  if (grepl("transcription factor", f4_kallisto_sig_up_v2$gene_name[i], fixed = T) | grepl("DNA-binding transcription factor activity", f4_kallisto_sig_up_v2$goslim.MF[i], fixed = TRUE)){
    f4_kallisto_sig_up_v2$category[i] <- paste("(possible TF)", f4_kallisto_sig_up_v2$category[i], sep = " ")
  }
}


#add the realtive protein expression in each sample 
proteome_abundances <- read.xlsx(file = "Proteome_final_tables/proteome_sup_tableV2.xlsx", sheetIndex = 1)


abundances_vac_1 <- proteome_abundances$vac_1_relative_abundance
names(abundances_vac_1) <- proteome_abundances$Caros.acession

abundances_vac_2 <- proteome_abundances$vac_2_relative_abundance
names(abundances_vac_2) <- proteome_abundances$Caros.acession

abundances_ton_1 <- proteome_abundances$ton_1_relative_abundance
names(abundances_ton_1) <- proteome_abundances$Caros.acession

abundances_ton_2 <- proteome_abundances$ton_2_relative_abundance
names(abundances_ton_2) <- proteome_abundances$Caros.acession

abundances_vac_mean <- proteome_abundances$vac_mean
names(abundances_vac_mean) <- proteome_abundances$Caros.acession

abundances_ton_mean <- proteome_abundances$ton_mean
names(abundances_ton_mean) <- proteome_abundances$Caros.acession

f4_kallisto_sig_up_v2$abudance_vac1 <- rep(NA, length(f4_kallisto_sig_up_v2$avgF1in))
f4_kallisto_sig_up_v2$abudance_vac2 <- rep(NA, length(f4_kallisto_sig_up_v2$avgF1in))
f4_kallisto_sig_up_v2$abudance_ton1 <- rep(NA, length(f4_kallisto_sig_up_v2$avgF1in))
f4_kallisto_sig_up_v2$abudance_ton2 <- rep(NA, length(f4_kallisto_sig_up_v2$avgF1in))
f4_kallisto_sig_up_v2$abudance_vac_mean <- rep(NA, length(f4_kallisto_sig_up_v2$avgF1in))
f4_kallisto_sig_up_v2$abudance_ton_mean <- rep(NA, length(f4_kallisto_sig_up_v2$avgF1in))



for (i in 1:length(f4_kallisto_sig_up_v2$avgF1in)){
  if (!is.na(f4_kallisto_sig_up_v2$proteome[i])){
    caros <- strsplit(f4_kallisto_sig_up_v2$proteome[i], "; ", fixed = T)[[1]]
    for (j in 1:length(caros)){
      caro_acession_numer <- strsplit(caros[j], " ")[[1]][1]
      if (j == 1){
        vac_1 <- abundances_vac_1[as.character(caro_acession_numer)]
        vac_2 <- abundances_vac_2[as.character(caro_acession_numer)]
        ton_1 <- abundances_ton_1[as.character(caro_acession_numer)]
        ton_2 <- abundances_ton_2[as.character(caro_acession_numer)]
        vac_mean <- abundances_vac_mean[as.character(caro_acession_numer)]
        ton_mean <- abundances_ton_mean[as.character(caro_acession_numer)]
      }
      if (j >= 2){
        vac_1 <- paste(vac_1, abundances_vac_1[as.character(caro_acession_numer)], sep = ";")
        vac_2 <- paste(vac_2, abundances_vac_2[as.character(caro_acession_numer)], sep = ";")
        ton_1 <- paste(ton_1, abundances_ton_1[as.character(caro_acession_numer)], sep = ";")
        ton_2 <- paste(ton_2, abundances_ton_2[as.character(caro_acession_numer)], sep = ";")
        vac_mean <- paste(vac_mean, abundances_vac_mean[as.character(caro_acession_numer)], sep = ";")
        ton_mean <- paste(ton_mean, abundances_ton_mean[as.character(caro_acession_numer)], sep = ";")
      }
    }
    f4_kallisto_sig_up_v2$abudance_vac1[i] <- vac_1
    f4_kallisto_sig_up_v2$abudance_vac2[i] <- vac_2
    f4_kallisto_sig_up_v2$abudance_ton1[i] <- ton_1
    f4_kallisto_sig_up_v2$abudance_ton2[i] <- ton_2
    f4_kallisto_sig_up_v2$abudance_vac_mean[i] <- vac_mean
    f4_kallisto_sig_up_v2$abudance_ton_mean[i] <- ton_mean
  }
}

f4_kallisto_sig_up_v2 <- f4_kallisto_sig_up_v2[,c(1:12,17,14:16,13,18:23)]
write.xlsx(f4_kallisto_sig_up_v2, file = paste(results_dir_f4, "f4in_vs_f4out_annotated.xlsx", sep = "/"), sheetName = "Up.V2.raw", append = T)

