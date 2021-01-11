#Some variables in this script may have to be loaded from run_WCGNA.R

#misc 
setwd("C:/Users/Faculdade/Desktop/Dissertação/results 2020/mapping_to_genome/11_WCGNA")


#library
library(xlsx)


#load the per module genes
modules_genes <- read.delim(file = "results/genes_per_modules.tsv")
table(modules_genes$module)

#Identify where are DEG genes 

#Load DEG genes
meso_idio_DEG_up <- read.delim(file = "DEG_genes/meso_vs_idio_sigUP.csv", sep = ",")
meso_idio_DEG_down <- read.delim(file = "DEG_genes/meso_vs_idio_sigDown.csv", sep = ",")
colnames(meso_idio_DEG_up) <- c("x", "geneSig")
colnames(meso_idio_DEG_down) <- c("x", "geneSig")
meso_idio_DEG <- rbind(meso_idio_DEG_up, meso_idio_DEG_down)


F1_DEG_up <- read.delim(file = "DEG_genes/F1in_vs_F1out_sigUP.csv", sep = ",")
F1_DEG_down <- read.delim(file = "DEG_genes/F1in_vs_F1out_sigDOWN.csv", sep = ",")
colnames(F1_DEG_up) <- c("x", "geneSig")
colnames(F1_DEG_down) <- c("x", "geneSig")
F1_DEG <- rbind(F1_DEG_up, F1_DEG_down)


F4_DEG_up <- read.delim(file = "DEG_genes/F4in_vs_F4out_sigUP.csv", sep = ",")
F4_DEG_down <- read.delim(file = "DEG_genes/F4in_vs_F4out_sigDOWN.csv", sep = ",")
colnames(F4_DEG_up) <- c("x", "geneSig")
colnames(F4_DEG_down) <- c("x", "geneSig")
F4_DEG <- rbind(F4_DEG_up, F4_DEG_down)


DEG_Up <- rbind(meso_idio_DEG_up, F1_DEG_up, F4_DEG_up)
DEG_UP <- unique(DEG_Up$geneSig)


DEG_Down <- rbind(meso_idio_DEG_down, F1_DEG_down, F4_DEG_down)
DEG_Down <- unique(DEG_Down$geneSig)


DEG <- rbind(meso_idio_DEG, F1_DEG, F4_DEG)
DEG <- unique(DEG$geneSig)


#Meso_vs_idio_DEG_UP---- 
meso_idio_up_modules <- modules_genes[as.character(meso_idio_DEG_up$geneSig), ]
table(meso_idio_up_modules$module)

#Meso_vs_idio_DEG_down----
meso_idio_down_modules <- modules_genes[as.character(meso_idio_DEG_down$geneSig), ]
table(meso_idio_down_modules$module)

#Meso_vs_idio_DEG_all----
meso_idio_modules <- modules_genes[as.character(meso_idio_DEG$geneSig), ]
table(meso_idio_modules$module)



#F1in_vs_F1out_DEG_UP---- 
F1_up_modules <- modules_genes[as.character(F1_DEG_up$geneSig), ]
table(F1_up_modules$module)

#F1in_vs_F1out_DEG_down----
F1_down_modules <- modules_genes[as.character(F1_DEG_down$geneSig), ]
table(F1_down_modules$module)

#F1in_vs_F1out_DEG_all----
F1_modules <- modules_genes[as.character(F1_DEG$geneSig), ]
table(F1_modules$module)


#F4in_vs_F4out_DEG_UP---- 
F4_up_modules <- modules_genes[as.character(F4_DEG_up$geneSig), ]
table(F4_up_modules$module)

#F4in_vs_F4out_DEG_down----
F4_down_modules <- modules_genes[as.character(F4_DEG_down$geneSig), ]
table(F4_down_modules$module)

#F4in_vs_F4out_DEG_all----
F4_modules <- modules_genes[as.character(F4_DEG$geneSig), ]
table(F4_modules$module)

#DEG_UP----
DEG_UP_module <- modules_genes[as.character(DEG_UP), ]
table(DEG_UP_module$module)

#DEG_DOWN---- 
DEG_Down_module <- modules_genes[as.character(DEG_Down), ]
table(DEG_Down_module$module)

#DEG
DEG_module <- modules_genes[as.character(DEG),]
table(DEG_module$module) 


#get modules with the pathway genes ---- 
manual_annotated_pathway <- read.xlsx(file = "annotation/manual_annotation_modified.xlsx", sheetIndex = 1)
manual_annotated_pathway <- manual_annotated_pathway[,c(2,4,3)]

manual_annotated_pathway$module <- modules_genes[as.character(manual_annotated_pathway$Gene_id),3]
manual_annotated_pathway$MM <- modules_genes[as.character(manual_annotated_pathway$Gene_id),4]
manual_annotated_pathway$MMPadj <- modules_genes[as.character(manual_annotated_pathway$Gene_id),5]
manual_annotated_pathway$GSvindoline <- modules_genes[as.character(manual_annotated_pathway$Gene_id),6]
manual_annotated_pathway$GSvindolinePadj <- modules_genes[as.character(manual_annotated_pathway$Gene_id),10]
manual_annotated_pathway$GSCatharanthina <- modules_genes[as.character(manual_annotated_pathway$Gene_id),7]
manual_annotated_pathway$GSCatharanthinaPadj <- modules_genes[as.character(manual_annotated_pathway$Gene_id),11]
manual_annotated_pathway$vinblastine <- modules_genes[as.character(manual_annotated_pathway$Gene_id),8]
manual_annotated_pathway$vinblastinePadj <- modules_genes[as.character(manual_annotated_pathway$Gene_id),12]
manual_annotated_pathway$AVLB <- modules_genes[as.character(manual_annotated_pathway$Gene_id),9]
manual_annotated_pathway$AVLBPadj <- modules_genes[as.character(manual_annotated_pathway$Gene_id),13]


  
write.xlsx(manual_annotated_pathway, file = "results/manual_annotation.xlsx")


#Get intersect 3 modules ---- 
intersect3 <- read.xlsx(file = "annotation/dammit_annotation_geral_kallisto.xlsx", sheetIndex = 1)

intersect3$module <- modules_genes[as.character(intersect3$NA.),3]
intersect3$module <- modules_genes[as.character(intersect3$NA.),3]
intersect3$MM <- modules_genes[as.character(intersect3$NA.),4]
intersect3$MMPadj <- modules_genes[as.character(intersect3$NA.),5]
intersect3$GSvindoline <- modules_genes[as.character(intersect3$NA.),6]
intersect3$GSvindolinePadj <- modules_genes[as.character(intersect3$NA.),10]
intersect3$GSCatharanthina <- modules_genes[as.character(intersect3$NA.),7]
intersect3$GSCatharanthinaPadj <- modules_genes[as.character(intersect3$NA.),11]
intersect3$vinblastine <- modules_genes[as.character(intersect3$NA.),8]
intersect3$vinblastinePadj <- modules_genes[as.character(intersect3$NA.),12]
intersect3$AVLB <- modules_genes[as.character(intersect3$NA.),9]
intersect3$AVLBPadj <- modules_genes[as.character(intersect3$NA.),13]


intersect3 <- intersect3[,c(1,3,4,24:34,10:17)]

write.xlsx(intersect3, file = "results/intersect_3_modules.xlsx")


#Make a table with the number of genes in each module and the number in of DEG (for each category) in each module ----

n_gene_mod <- table(modules_genes$module)
mod_no_dataframe <- data.frame(n_gene_mod)
colnames(mod_no_dataframe)[2] <- "Geral" 

#add DEG
temp_data_frame <- data.frame(table(DEG_module$module))
colnames(temp_data_frame)[2] <- "DEG"
mod_no_dataframe <- merge(mod_no_dataframe, temp_data_frame, by = 1, all.x = T)

mod_no_dataframe$DEG_percent <- round(100 * (mod_no_dataframe$DEG / mod_no_dataframe$Geral),2)


#add DEG up
temp_data_frame <- data.frame(table(DEG_UP_module$module))
colnames(temp_data_frame)[2] <- "DEG_up"
mod_no_dataframe <- merge(mod_no_dataframe, temp_data_frame, by = 1, all.x = T)

mod_no_dataframe$DEG_up_percent <- round(100 * (mod_no_dataframe$DEG_up / mod_no_dataframe$Geral),2)


#add DEG down
temp_data_frame <- data.frame(table(DEG_Down_module$module))
colnames(temp_data_frame)[2] <- "DEG_down"
mod_no_dataframe <- merge(mod_no_dataframe, temp_data_frame, by = 1, all.x = T)

mod_no_dataframe$DEG_down_percent <- round(100 * (mod_no_dataframe$DEG_down / mod_no_dataframe$Geral),2)


#add meso vs idio DEG
temp_data_frame <- data.frame(table(meso_idio_modules$module))
colnames(temp_data_frame)[2] <- "meso_vs_idio_DEG"
mod_no_dataframe <- merge(mod_no_dataframe, temp_data_frame, by = 1, all.x = T)

mod_no_dataframe$idio_DEG_percent <- round(100 * (mod_no_dataframe$meso_vs_idio_DEG / mod_no_dataframe$Geral),2)


#add meso vs idio up DEG
temp_data_frame <- data.frame(table(meso_idio_up_modules$module))
colnames(temp_data_frame)[2] <- "meso_vs_idio_up_DEG"
mod_no_dataframe <- merge(mod_no_dataframe, temp_data_frame, by = 1, all.x = T)

mod_no_dataframe$idio_up_DEG_percent <- round(100 * (mod_no_dataframe$meso_vs_idio_up_DEG / mod_no_dataframe$Geral),2)


#add meso vs idio down DEG
temp_data_frame <- data.frame(table(meso_idio_down_modules$module))
colnames(temp_data_frame)[2] <- "meso_vs_idio_down_DEG"
mod_no_dataframe <- merge(mod_no_dataframe, temp_data_frame, by = 1, all.x = T)

mod_no_dataframe$idio_down_DEG_percent <- round(100 * (mod_no_dataframe$meso_vs_idio_down_DEG / mod_no_dataframe$Geral),2)

  
#add F1in vs F1out DEG
temp_data_frame <- data.frame(table(F1_modules$module))
colnames(temp_data_frame)[2] <- "F1_DEG"
mod_no_dataframe <- merge(mod_no_dataframe, temp_data_frame, by = 1, all.x = T)

mod_no_dataframe$F1_DEG_percent <- round(100 * (mod_no_dataframe$F1_DEG / mod_no_dataframe$Geral),2)


#add F1in vs F1out up DEG
temp_data_frame <- data.frame(table(F1_up_modules$module))
colnames(temp_data_frame)[2] <- "F1_up_DEG"
mod_no_dataframe <- merge(mod_no_dataframe, temp_data_frame, by = 1, all.x = T)

mod_no_dataframe$F1_up_DEG_percent <- round(100 * (mod_no_dataframe$F1_up_DEG / mod_no_dataframe$Geral),2)


#add F1in vs F1out down DEG
temp_data_frame <- data.frame(table(F1_down_modules$module))
colnames(temp_data_frame)[2] <- "F1_down_DEG"
mod_no_dataframe <- merge(mod_no_dataframe, temp_data_frame, by = 1, all.x = T)

mod_no_dataframe$F1_down_DEG_percent <- round(100 * (mod_no_dataframe$F1_down_DEG / mod_no_dataframe$Geral),2)



#add F4in vs F4out DEG
temp_data_frame <- data.frame(table(F4_modules$module))
colnames(temp_data_frame)[2] <- "F4_DEG"
mod_no_dataframe <- merge(mod_no_dataframe, temp_data_frame, by = 1, all.x = T)

mod_no_dataframe$F4_DEG_percent <- round(100 * (mod_no_dataframe$F4_DEG / mod_no_dataframe$Geral),2)


#add F4in vs F4out up DEG
temp_data_frame <- data.frame(table(F4_up_modules$module))
colnames(temp_data_frame)[2] <- "F4_up_DEG"
mod_no_dataframe <- merge(mod_no_dataframe, temp_data_frame, by = 1, all.x = T)

mod_no_dataframe$F4_up_DEG_percent <- round(100 * (mod_no_dataframe$F4_up_DEG / mod_no_dataframe$Geral),2)


#add F4in vs F4out down DEG
temp_data_frame <- data.frame(table(F4_down_modules$module))
colnames(temp_data_frame)[2] <- "F4_down_DEG"
mod_no_dataframe <- merge(mod_no_dataframe, temp_data_frame, by = 1, all.x = T)

mod_no_dataframe$F4_down_DEG_percent <- round(100 * (mod_no_dataframe$F4_down_DEG / mod_no_dataframe$Geral),2)

mod_no_dataframe[is.na(mod_no_dataframe)] <- 0 

#save table 
write.xlsx(mod_no_dataframe, file = "results/number_genes_DEG_per_module.xlsx", row.names = F)

#make a smaller version of this table

mod_no_dataframe.abv <- mod_no_dataframe[,c(1,2)]
mod_no_dataframe.abv$DEG <- paste(mod_no_dataframe$DEG, " ", "(", mod_no_dataframe$DEG_percent, "%", ")", sep = "" )
mod_no_dataframe.abv$DEG_idio <- paste(mod_no_dataframe$meso_vs_idio_DEG, " ", "(", mod_no_dataframe$idio_DEG_percent, "%", ")", sep = "" )
mod_no_dataframe.abv$DEG_F1 <- paste(mod_no_dataframe$F1_DEG, " ", "(", mod_no_dataframe$F1_DEG_percent, "%", ")", sep = "" )
mod_no_dataframe.abv$DEG_F4 <- paste(mod_no_dataframe$F4_DEG, " ", "(", mod_no_dataframe$F4_DEG_percent, "%", ")", sep = "" )

write.xlsx(mod_no_dataframe.abv, file = "results/number_genes_DEG_per_module.xlsx", sheetName = "Abv table", append = T)


####Make an excel file with each sheat containing genes of one module---- 


allmodules <- names(table(modules_genes$module))

flag = T
for (i in allmodules){
  temp_table <- subset(modules_genes, modules_genes$module == i)
  temp_table$up_idio <- ifelse(row.names(temp_table) %in% meso_idio_DEG_up$geneSig, T, F)
  temp_table$up_F1 <- ifelse(row.names(temp_table) %in% F1_DEG_up$geneSig, T, F)
  temp_table$up_F4 <- ifelse(row.names(temp_table) %in% F4_DEG_up$geneSig, T, F)
  if (flag == T){
    write.xlsx(temp_table, file = "results/genesPermodule.xlsx", sheetName = i)
    flag = F
  }
  else{
    write.xlsx(temp_table, file = "results/genesPermodule.xlsx", append = T, sheetName = i)
    gc()
  }
}



###make a third table (to go to thesis)

table_geral <- read.xlsx(file = "results/number_genes_DEG_per_module.xlsx", sheetIndex = 4)

table_final <- table_geral[,c(1:5)]

table_final$DEG.F1.up <- paste0(table_geral$F1_up_DEG, " (" ,table_geral$F1_up_DEG_percent, "%)")
table_final$DEG.F1.down <- paste0(table_geral$F1_down_DEG, " (" ,table_geral$F1_down_DEG_percent, "%)")
table_final$DEG.F4.up <- paste0(table_geral$F4_up_DEG, " (" ,table_geral$F4_up_DEG_percent, "%)")
table_final$DEG.F4.down <- paste0(table_geral$F4_down_DEG, " (" ,table_geral$F4_down_DEG_percent, "%)")

write.xlsx(table_final, file = "results/number_genes_DEG_per_module.xlsx", sheetName = "Thesis table",append = TRUE)
