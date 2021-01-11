#Note before running: Take care as in one the of the steps bellow there is some comented lines
#these were comented due to errors I could not solve, and it does not make much difference
#If running with other data (not the WGCNA modules produces in 30/09, take care in this section)
#i.e uncoment lines


#misc 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

#load libraries
library(topGO)
library(ggplot2)
library(stringr)
library(xlsx)

#Load the datatable
dammit_annot.modules <- read.delim(file = "results/genes_per_modules.tsv")
dammit_annot.modules$gene_id <- row.names(dammit_annot.modules)


#result dir 
results.dir <- "results/topGO_enrichement"
if (!dir.exists(results.dir)){
  dir.create(results.dir)
}

#supplementary tables dir 
supp.dir <- "results/topGO_enrichement/supplementary"
if (!dir.exists(supp.dir)){
  dir.create(supp.dir)
}




#make background geneset ----
#using as background the 15204 genes 

##Go(P) <- biological process
background_Go.P <- data.frame("gene" = dammit_annot.modules$gene_id, "GO.BP" = dammit_annot.modules$GO.BP.)
background_Go.P$'GO.BP' <- str_replace_all(background_Go.P$'GO.BP', ";", ",")
write.table(background_Go.P, file = "background_GoP.txt", row.names = FALSE, col.names = F,quote = F, sep = "\t")


##Go(F) <- molecular function
background_Go.F <- data.frame("gene" = dammit_annot.modules$gene_id, "GO.MF" = dammit_annot.modules$GO.MF.)
background_Go.F$'GO.MF' <- str_replace_all(background_Go.F$'GO.MF', ";", ",")
write.table(background_Go.F, file = "background_GoF.txt", row.names = FALSE, col.names = F,quote = F, sep = "\t")

##Go(C) <- Cellular component
background_Go.C <- data.frame("gene" = dammit_annot.modules$gene_id, "GO.CC" = dammit_annot.modules$GO.CC.)
background_Go.C$'GO.CC' <- str_replace_all(background_Go.C$'GO.CC', ";", ",")
write.table(background_Go.C, file = "background_GoC.txt", row.names = FALSE, col.names = F,quote = F, sep = "\t")



#Load blackground annotation 
background.P <- readMappings(file = "background_GoP.txt")
bg_genes.P <- names(background.P)

background.F <- readMappings(file = "background_GoF.txt")
bg_genes.F <- names(background.F)

background.C <- readMappings(file = "background_GoC.txt")
bg_genes.C <- names(background.C)


table(dammit_annot.modules$module)


#darkgreen module ---- 
darkgreen_module <- subset(dammit_annot.modules, dammit_annot.modules$module == "darkgreen")

#GO.P
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.P %in% darkgreen_module$gene_id))
names(compared_genes) <- bg_genes.P

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_darkgreen.P <- new("topGOdata", ontology = "BP", 
                             description = "darkgreen module enrichement analysis",
                             allGenes = compared_genes, 
                             annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                             gene2GO = background.P,
                             nodeSize = 5) 

#Run fisher test
resultFisher_darkgreen.P <- runTest(GoData_darkgreen.P, 
                                  algorithm = "weight01", 
                                  statistic = "fisher")
resultFisher_darkgreen.P

allRes_darkgreen.ByP <- GenTable(GoData_darkgreen.P, 
                                       classicFisher =resultFisher_darkgreen.P,
                                       orderBy= "classicFisher", 
                                       numChar = 100, 
                                       topNodes = length(resultFisher_darkgreen.P@score))


allRes_darkgreen.ByP.sig <- allRes_darkgreen.ByP[which(as.numeric(allRes_darkgreen.ByP$classicFisher) < 0.01),]
allRes_darkgreen.ByP.sig <- allRes_darkgreen.ByP.sig[order(allRes_darkgreen.ByP.sig$Significant, decreasing = T), ]


#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% darkgreen_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_darkgreen.F <- new("topGOdata", ontology = "MF", 
                          description = "Black module enrichement analysis",
                          allGenes = compared_genes, 
                          annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                          gene2GO = background.F,
                          nodeSize = 5) 

#Run fisher test
resultFisher_darkgreen.F <- runTest(GoData_darkgreen.F, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_darkgreen.F

allRes_darkgreen.MoF <- GenTable(GoData_darkgreen.F, 
                                 classicFisher =resultFisher_darkgreen.F,
                                 orderBy= "classicFisher", 
                                 numChar = 100, 
                                 topNodes = length(resultFisher_darkgreen.F@score))


allRes_darkgreen.MoF.sig <- allRes_darkgreen.MoF[which(as.numeric(allRes_darkgreen.MoF$classicFisher) < 0.01),]
allRes_darkgreen.MoF.sig <- allRes_darkgreen.MoF.sig[order(allRes_darkgreen.MoF.sig$Significant, decreasing = T), ]



#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% darkgreen_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_darkgreen.F <- new("topGOdata", ontology = "MF", 
                          description = "darkgreen module enrichement analysis",
                          allGenes = compared_genes, 
                          annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                          gene2GO = background.F,
                          nodeSize = 5) 

#Run fisher test
resultFisher_darkgreen.F <- runTest(GoData_darkgreen.F, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_darkgreen.F

allRes_darkgreen.MoF <- GenTable(GoData_darkgreen.F, 
                                 classicFisher =resultFisher_darkgreen.F,
                                 orderBy= "classicFisher", 
                                 numChar = 100, 
                                 topNodes = length(resultFisher_darkgreen.F@score))


allRes_darkgreen.MoF.sig <- allRes_darkgreen.MoF[which(as.numeric(allRes_darkgreen.MoF$classicFisher) < 0.01),]
allRes_darkgreen.MoF.sig <- allRes_darkgreen.MoF.sig[order(allRes_darkgreen.MoF.sig$Significant, decreasing = T), ]



#GO.C
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.C %in% darkgreen_module$gene_id))
names(compared_genes) <- bg_genes.C

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_darkgreen.C <- new("topGOdata", ontology = "CC", 
                          description = "darkgreen module enrichement analysis",
                          allGenes = compared_genes, 
                          annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                          gene2GO = background.C,
                          nodeSize = 5) 

#Run fisher test
resultFisher_darkgreen.C <- runTest(GoData_darkgreen.C, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_darkgreen.C

allRes_darkgreen.CeC <- GenTable(GoData_darkgreen.C, 
                                 classicFisher =resultFisher_darkgreen.C,
                                 orderBy= "classicFisher", 
                                 numChar = 100, 
                                 topNodes = length(resultFisher_darkgreen.C@score))


allRes_darkgreen.CeC.sig <- allRes_darkgreen.CeC[which(as.numeric(allRes_darkgreen.CeC$classicFisher) < 0.01),]
allRes_darkgreen.CeC.sig <- allRes_darkgreen.CeC.sig[order(allRes_darkgreen.CeC.sig$Significant, decreasing = T), ]


allRes_darkgreen.combined <- rbind(allRes_darkgreen.ByP.sig, allRes_darkgreen.MoF.sig, allRes_darkgreen.CeC.sig)
allRes_darkgreen.combined$order <- seq(from = length(allRes_darkgreen.combined$GO.ID), to = 1, by = -1)

darkgreen_plot <- ggplot(allRes_darkgreen.combined, aes(x = reorder(Term, order, decreasing = F), y = Significant)) + 
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "darkgreen") + 
  coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.2) + 
  labs(title = "Module 14", 
       subtitle = "Enriched GO terms", 
       x = "", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))




png("results/topGO_enrichement/darkgreenModule.png", res = 300, width = 4500, height = 3000)
darkgreen_plot
dev.off()


#Make tables with the significant GO terms


allRes_darkgreen.ByP.sig$genes <- sapply(allRes_darkgreen.ByP.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_darkgreen.P, x)
  genes[[1]][genes[[1]] %in% row.names(darkgreen_module)]
})

allRes_darkgreen.MoF.sig$genes <- sapply(allRes_darkgreen.MoF.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_darkgreen.F, x)
  genes[[1]][genes[[1]] %in% row.names(darkgreen_module)]
})

allRes_darkgreen.CeC.sig$genes <- sapply(allRes_darkgreen.CeC.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_darkgreen.C, x)
  genes[[1]][genes[[1]] %in% row.names(darkgreen_module)]
})

allRes_darkgreen.genes.combined <- rbind(allRes_darkgreen.ByP.sig, allRes_darkgreen.MoF.sig, allRes_darkgreen.CeC.sig)

GoTerm_genes_darkgreen <- "IDs_darkgreen.xlsx"


for (i in 1:length(allRes_darkgreen.genes.combined$GO.ID)){
  genes_id <- allRes_darkgreen.genes.combined$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_darkgreen.genes.combined$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- dammit_annot.modules[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(24,25,1:23)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/topGO_enrichement/", GoTerm_genes_darkgreen), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/topGO_enrichement/", GoTerm_genes_darkgreen), sheetName = str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}

write.csv(allRes_darkgreen.combined[,-7], file = paste0(supp.dir, "/module_14.tsv"), row.names = F)



#yellowgreen module ---- 
yellowgreen_module <- subset(dammit_annot.modules, dammit_annot.modules$module == "yellowgreen")

#GO.P
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.P %in% yellowgreen_module$gene_id))
names(compared_genes) <- bg_genes.P

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_yellowgreen.P <- new("topGOdata", ontology = "BP", 
                          description = "yellowgreen module enrichement analysis",
                          allGenes = compared_genes, 
                          annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                          gene2GO = background.P,
                          nodeSize = 5) 

#Run fisher test
resultFisher_yellowgreen.P <- runTest(GoData_yellowgreen.P, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_yellowgreen.P

allRes_yellowgreen.ByP <- GenTable(GoData_yellowgreen.P, 
                                 classicFisher =resultFisher_yellowgreen.P,
                                 orderBy= "classicFisher", 
                                 numChar = 100, 
                                 topNodes = length(resultFisher_yellowgreen.P@score))


allRes_yellowgreen.ByP.sig <- allRes_yellowgreen.ByP[which(as.numeric(allRes_yellowgreen.ByP$classicFisher) < 0.01),]
allRes_yellowgreen.ByP.sig <- allRes_yellowgreen.ByP.sig[order(allRes_yellowgreen.ByP.sig$Significant, decreasing = T), ]


#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% yellowgreen_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_yellowgreen.F <- new("topGOdata", ontology = "MF", 
                          description = "Black module enrichement analysis",
                          allGenes = compared_genes, 
                          annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                          gene2GO = background.F,
                          nodeSize = 5) 

#Run fisher test
resultFisher_yellowgreen.F <- runTest(GoData_yellowgreen.F, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_yellowgreen.F

allRes_yellowgreen.MoF <- GenTable(GoData_yellowgreen.F, 
                                 classicFisher =resultFisher_yellowgreen.F,
                                 orderBy= "classicFisher", 
                                 numChar = 100, 
                                 topNodes = length(resultFisher_yellowgreen.F@score))


allRes_yellowgreen.MoF.sig <- allRes_yellowgreen.MoF[which(as.numeric(allRes_yellowgreen.MoF$classicFisher) < 0.01),]
allRes_yellowgreen.MoF.sig <- allRes_yellowgreen.MoF.sig[order(allRes_yellowgreen.MoF.sig$Significant, decreasing = T), ]



#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% yellowgreen_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_yellowgreen.F <- new("topGOdata", ontology = "MF", 
                          description = "yellowgreen module enrichement analysis",
                          allGenes = compared_genes, 
                          annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                          gene2GO = background.F,
                          nodeSize = 5) 

#Run fisher test
resultFisher_yellowgreen.F <- runTest(GoData_yellowgreen.F, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_yellowgreen.F

allRes_yellowgreen.MoF <- GenTable(GoData_yellowgreen.F, 
                                 classicFisher =resultFisher_yellowgreen.F,
                                 orderBy= "classicFisher", 
                                 numChar = 100, 
                                 topNodes = length(resultFisher_yellowgreen.F@score))


allRes_yellowgreen.MoF.sig <- allRes_yellowgreen.MoF[which(as.numeric(allRes_yellowgreen.MoF$classicFisher) < 0.01),]
allRes_yellowgreen.MoF.sig <- allRes_yellowgreen.MoF.sig[order(allRes_yellowgreen.MoF.sig$Significant, decreasing = T), ]



#GO.C
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.C %in% yellowgreen_module$gene_id))
names(compared_genes) <- bg_genes.C

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_yellowgreen.C <- new("topGOdata", ontology = "CC", 
                          description = "yellowgreen module enrichement analysis",
                          allGenes = compared_genes, 
                          annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                          gene2GO = background.C,
                          nodeSize = 5) 

#Run fisher test
resultFisher_yellowgreen.C <- runTest(GoData_yellowgreen.C, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_yellowgreen.C

allRes_yellowgreen.CeC <- GenTable(GoData_yellowgreen.C, 
                                 classicFisher =resultFisher_yellowgreen.C,
                                 orderBy= "classicFisher", 
                                 numChar = 100, 
                                 topNodes = length(resultFisher_yellowgreen.C@score))


allRes_yellowgreen.CeC.sig <- allRes_yellowgreen.CeC[which(as.numeric(allRes_yellowgreen.CeC$classicFisher) < 0.01),]
allRes_yellowgreen.CeC.sig <- allRes_yellowgreen.CeC.sig[order(allRes_yellowgreen.CeC.sig$Significant, decreasing = T), ]


allRes_yellowgreen.combined <- rbind(allRes_yellowgreen.ByP.sig, allRes_yellowgreen.MoF.sig, allRes_yellowgreen.CeC.sig)
allRes_yellowgreen.combined$order <- seq(from = length(allRes_yellowgreen.combined$GO.ID), to = 1, by = -1)

yellowgreen_plot <- ggplot(allRes_yellowgreen.combined, aes(x = reorder(Term, order, decreasing = F), y = Significant)) + 
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "yellowgreen") + 
  coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.2) + 
  labs(title = "Module 15", 
       subtitle = "Enriched GO terms", 
       x = "", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))




png("results/topGO_enrichement/yellowgreenModule.png", res = 300, width = 4500, height = 3000)
yellowgreen_plot
dev.off()


#Make tables with the significant GO terms


allRes_yellowgreen.ByP.sig$genes <- sapply(allRes_yellowgreen.ByP.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_yellowgreen.P, x)
  genes[[1]][genes[[1]] %in% row.names(yellowgreen_module)]
})

allRes_yellowgreen.MoF.sig$genes <- sapply(allRes_yellowgreen.MoF.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_yellowgreen.F, x)
  genes[[1]][genes[[1]] %in% row.names(yellowgreen_module)]
})

allRes_yellowgreen.CeC.sig$genes <- sapply(allRes_yellowgreen.CeC.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_yellowgreen.C, x)
  genes[[1]][genes[[1]] %in% row.names(yellowgreen_module)]
})

allRes_yellowgreen.genes.combined <- rbind(allRes_yellowgreen.ByP.sig, allRes_yellowgreen.MoF.sig, allRes_yellowgreen.CeC.sig)

GoTerm_genes_yellowgreen <- "IDs_yellowgreen.xlsx"


for (i in 1:length(allRes_yellowgreen.genes.combined$GO.ID)){
  genes_id <- allRes_yellowgreen.genes.combined$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_yellowgreen.genes.combined$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- dammit_annot.modules[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(24,25,1:23)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/topGO_enrichement/", GoTerm_genes_yellowgreen), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/topGO_enrichement/", GoTerm_genes_yellowgreen), sheetName = str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}

write.csv(allRes_yellowgreen.combined[,-7], file = paste0(supp.dir, "/module_14.tsv"), row.names = F)




#darkmagenta module ---- 
darkmagenta_module <- subset(dammit_annot.modules, dammit_annot.modules$module == "darkmagenta")

#GO.P
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.P %in% darkmagenta_module$gene_id))
names(compared_genes) <- bg_genes.P

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_darkmagenta.P <- new("topGOdata", ontology = "BP", 
                          description = "darkmagenta module enrichement analysis",
                          allGenes = compared_genes, 
                          annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                          gene2GO = background.P,
                          nodeSize = 5) 

#Run fisher test
resultFisher_darkmagenta.P <- runTest(GoData_darkmagenta.P, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_darkmagenta.P

allRes_darkmagenta.ByP <- GenTable(GoData_darkmagenta.P, 
                                 classicFisher =resultFisher_darkmagenta.P,
                                 orderBy= "classicFisher", 
                                 numChar = 100, 
                                 topNodes = length(resultFisher_darkmagenta.P@score))


allRes_darkmagenta.ByP.sig <- allRes_darkmagenta.ByP[which(as.numeric(allRes_darkmagenta.ByP$classicFisher) < 0.01),]
allRes_darkmagenta.ByP.sig <- allRes_darkmagenta.ByP.sig[order(allRes_darkmagenta.ByP.sig$Significant, decreasing = T), ]


#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% darkmagenta_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_darkmagenta.F <- new("topGOdata", ontology = "MF", 
                          description = "Black module enrichement analysis",
                          allGenes = compared_genes, 
                          annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                          gene2GO = background.F,
                          nodeSize = 5) 

#Run fisher test
resultFisher_darkmagenta.F <- runTest(GoData_darkmagenta.F, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_darkmagenta.F

allRes_darkmagenta.MoF <- GenTable(GoData_darkmagenta.F, 
                                 classicFisher =resultFisher_darkmagenta.F,
                                 orderBy= "classicFisher", 
                                 numChar = 100, 
                                 topNodes = length(resultFisher_darkmagenta.F@score))


allRes_darkmagenta.MoF.sig <- allRes_darkmagenta.MoF[which(as.numeric(allRes_darkmagenta.MoF$classicFisher) < 0.01),]
allRes_darkmagenta.MoF.sig <- allRes_darkmagenta.MoF.sig[order(allRes_darkmagenta.MoF.sig$Significant, decreasing = T), ]



#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% darkmagenta_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_darkmagenta.F <- new("topGOdata", ontology = "MF", 
                          description = "darkmagenta module enrichement analysis",
                          allGenes = compared_genes, 
                          annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                          gene2GO = background.F,
                          nodeSize = 5) 

#Run fisher test
resultFisher_darkmagenta.F <- runTest(GoData_darkmagenta.F, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_darkmagenta.F

allRes_darkmagenta.MoF <- GenTable(GoData_darkmagenta.F, 
                                 classicFisher =resultFisher_darkmagenta.F,
                                 orderBy= "classicFisher", 
                                 numChar = 100, 
                                 topNodes = length(resultFisher_darkmagenta.F@score))


allRes_darkmagenta.MoF.sig <- allRes_darkmagenta.MoF[which(as.numeric(allRes_darkmagenta.MoF$classicFisher) < 0.01),]
allRes_darkmagenta.MoF.sig <- allRes_darkmagenta.MoF.sig[order(allRes_darkmagenta.MoF.sig$Significant, decreasing = T), ]



#GO.C
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.C %in% darkmagenta_module$gene_id))
names(compared_genes) <- bg_genes.C

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_darkmagenta.C <- new("topGOdata", ontology = "CC", 
                          description = "darkmagenta module enrichement analysis",
                          allGenes = compared_genes, 
                          annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                          gene2GO = background.C,
                          nodeSize = 5) 

#Run fisher test
resultFisher_darkmagenta.C <- runTest(GoData_darkmagenta.C, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_darkmagenta.C

allRes_darkmagenta.CeC <- GenTable(GoData_darkmagenta.C, 
                                 classicFisher =resultFisher_darkmagenta.C,
                                 orderBy= "classicFisher", 
                                 numChar = 100, 
                                 topNodes = length(resultFisher_darkmagenta.C@score))


allRes_darkmagenta.CeC.sig <- allRes_darkmagenta.CeC[which(as.numeric(allRes_darkmagenta.CeC$classicFisher) < 0.01),]
allRes_darkmagenta.CeC.sig <- allRes_darkmagenta.CeC.sig[order(allRes_darkmagenta.CeC.sig$Significant, decreasing = T), ]


allRes_darkmagenta.combined <- rbind(allRes_darkmagenta.ByP.sig, allRes_darkmagenta.MoF.sig, allRes_darkmagenta.CeC.sig)
allRes_darkmagenta.combined$order <- seq(from = length(allRes_darkmagenta.combined$GO.ID), to = 1, by = -1)

darkmagenta_plot <- ggplot(allRes_darkmagenta.combined, aes(x = reorder(Term, order, decreasing = F), y = Significant)) + 
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "darkmagenta") + 
  coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.2) + 
  labs(title = "Module 1", 
       subtitle = "Enriched GO terms", 
       x = "", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))




png("results/topGO_enrichement/darkmagentaModule.png", res = 300, width = 4500, height = 3000)
darkmagenta_plot
dev.off()


#Make tables with the significant GO terms


allRes_darkmagenta.ByP.sig$genes <- sapply(allRes_darkmagenta.ByP.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_darkmagenta.P, x)
  genes[[1]][genes[[1]] %in% row.names(darkmagenta_module)]
  print(genes[[1]][genes[[1]] %in% row.names(darkmagenta_module)])
})

#For some reason the above code does not work

genes <- genesInTerm(GoData_darkmagenta.P, allRes_darkmagenta.ByP.sig$GO.ID[1])
temp <- genes[[1]][genes[[1]] %in% row.names(darkmagenta_module)]
allRes_darkmagenta.ByP.sig$genes[1] <- list(temp)

genes <- genesInTerm(GoData_darkmagenta.P, allRes_darkmagenta.ByP.sig$GO.ID[2])
temp <- genes[[1]][genes[[1]] %in% row.names(darkmagenta_module)]
allRes_darkmagenta.ByP.sig$genes[2] <- list(temp)

genes <- genesInTerm(GoData_darkmagenta.P, allRes_darkmagenta.ByP.sig$GO.ID[3])
temp <- genes[[1]][genes[[1]] %in% row.names(darkmagenta_module)]
allRes_darkmagenta.ByP.sig$genes[3] <- list(temp)
#also change some code bellow




allRes_darkmagenta.MoF.sig$genes <- sapply(allRes_darkmagenta.MoF.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_darkmagenta.F, x)
  genes[[1]][genes[[1]] %in% row.names(darkmagenta_module)]
})

allRes_darkmagenta.CeC.sig$genes <- sapply(allRes_darkmagenta.CeC.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_darkmagenta.C, x)
  genes[[1]][genes[[1]] %in% row.names(darkmagenta_module)]
})


allRes_darkmagenta.genes.combined <- rbind(allRes_darkmagenta.ByP.sig, allRes_darkmagenta.MoF.sig, allRes_darkmagenta.CeC.sig)

GoTerm_genes_darkmagenta <- "IDs_darkmagenta.xlsx"


for (i in 1:length(allRes_darkmagenta.genes.combined$GO.ID)){
  genes_id <- allRes_darkmagenta.genes.combined$genes[i] #get the genes from each term 
  Go_term <- allRes_darkmagenta.genes.combined$GO.ID[i] #Get the Go term code 
  Go_name <- allRes_darkmagenta.genes.combined$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- dammit_annot.modules[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(24,25,1:23)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/topGO_enrichement/", GoTerm_genes_darkmagenta), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/topGO_enrichement/", GoTerm_genes_darkmagenta), sheetName = str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}

View(allRes_darkmagenta.combined)

write.csv(allRes_darkmagenta.combined[,-7], file = paste0(supp.dir, "/module_1.tsv"), row.names = F)

#salmon module ---- 
salmon_module <- subset(dammit_annot.modules, dammit_annot.modules$module == "salmon")

#GO.P
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.P %in% salmon_module$gene_id))
names(compared_genes) <- bg_genes.P

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_salmon.P <- new("topGOdata", ontology = "BP", 
                            description = "salmon module enrichement analysis",
                            allGenes = compared_genes, 
                            annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                            gene2GO = background.P,
                            nodeSize = 5) 

#Run fisher test
resultFisher_salmon.P <- runTest(GoData_salmon.P, 
                                      algorithm = "weight01", 
                                      statistic = "fisher")
resultFisher_salmon.P

allRes_salmon.ByP <- GenTable(GoData_salmon.P, 
                                   classicFisher =resultFisher_salmon.P,
                                   orderBy= "classicFisher", 
                                   numChar = 100, 
                                   topNodes = length(resultFisher_salmon.P@score))


allRes_salmon.ByP.sig <- allRes_salmon.ByP[which(as.numeric(allRes_salmon.ByP$classicFisher) < 0.01),]
allRes_salmon.ByP.sig <- allRes_salmon.ByP.sig[order(allRes_salmon.ByP.sig$Significant, decreasing = T), ]


#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% salmon_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_salmon.F <- new("topGOdata", ontology = "MF", 
                            description = "Black module enrichement analysis",
                            allGenes = compared_genes, 
                            annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                            gene2GO = background.F,
                            nodeSize = 5) 

#Run fisher test
resultFisher_salmon.F <- runTest(GoData_salmon.F, 
                                      algorithm = "weight01", 
                                      statistic = "fisher")
resultFisher_salmon.F

allRes_salmon.MoF <- GenTable(GoData_salmon.F, 
                                   classicFisher =resultFisher_salmon.F,
                                   orderBy= "classicFisher", 
                                   numChar = 100, 
                                   topNodes = length(resultFisher_salmon.F@score))


allRes_salmon.MoF.sig <- allRes_salmon.MoF[which(as.numeric(allRes_salmon.MoF$classicFisher) < 0.01),]
allRes_salmon.MoF.sig <- allRes_salmon.MoF.sig[order(allRes_salmon.MoF.sig$Significant, decreasing = T), ]



#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% salmon_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_salmon.F <- new("topGOdata", ontology = "MF", 
                            description = "salmon module enrichement analysis",
                            allGenes = compared_genes, 
                            annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                            gene2GO = background.F,
                            nodeSize = 5) 

#Run fisher test
resultFisher_salmon.F <- runTest(GoData_salmon.F, 
                                      algorithm = "weight01", 
                                      statistic = "fisher")
resultFisher_salmon.F

allRes_salmon.MoF <- GenTable(GoData_salmon.F, 
                                   classicFisher =resultFisher_salmon.F,
                                   orderBy= "classicFisher", 
                                   numChar = 100, 
                                   topNodes = length(resultFisher_salmon.F@score))


allRes_salmon.MoF.sig <- allRes_salmon.MoF[which(as.numeric(allRes_salmon.MoF$classicFisher) < 0.01),]
allRes_salmon.MoF.sig <- allRes_salmon.MoF.sig[order(allRes_salmon.MoF.sig$Significant, decreasing = T), ]



#GO.C
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.C %in% salmon_module$gene_id))
names(compared_genes) <- bg_genes.C

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_salmon.C <- new("topGOdata", ontology = "CC", 
                            description = "salmon module enrichement analysis",
                            allGenes = compared_genes, 
                            annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                            gene2GO = background.C,
                            nodeSize = 5) 

#Run fisher test
resultFisher_salmon.C <- runTest(GoData_salmon.C, 
                                      algorithm = "weight01", 
                                      statistic = "fisher")
resultFisher_salmon.C

allRes_salmon.CeC <- GenTable(GoData_salmon.C, 
                                   classicFisher =resultFisher_salmon.C,
                                   orderBy= "classicFisher", 
                                   numChar = 100, 
                                   topNodes = length(resultFisher_salmon.C@score))


allRes_salmon.CeC.sig <- allRes_salmon.CeC[which(as.numeric(allRes_salmon.CeC$classicFisher) < 0.01),]
allRes_salmon.CeC.sig <- allRes_salmon.CeC.sig[order(allRes_salmon.CeC.sig$Significant, decreasing = T), ]


allRes_salmon.combined <- rbind(allRes_salmon.ByP.sig, allRes_salmon.MoF.sig, allRes_salmon.CeC.sig)
allRes_salmon.combined$order <- seq(from = length(allRes_salmon.combined$GO.ID), to = 1, by = -1)

salmon_plot <- ggplot(allRes_salmon.combined, aes(x = reorder(Term, order, decreasing = F), y = Significant)) + 
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "salmon") + 
  coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.2) + 
  labs(title = "Module 2", 
       subtitle = "Enriched GO terms", 
       x = "", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))




png("results/topGO_enrichement/salmonModule.png", res = 300, width = 4500, height = 3000)
salmon_plot
dev.off()



#Make tables with the significant GO terms


allRes_salmon.ByP.sig$genes <- sapply(allRes_salmon.ByP.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_salmon.P, x)
  genes[[1]][genes[[1]] %in% row.names(salmon_module)]
})


allRes_salmon.MoF.sig$genes <- sapply(allRes_salmon.MoF.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_salmon.F, x)
  genes[[1]][genes[[1]] %in% row.names(salmon_module)]
})



#allRes_salmon.CeC.sig$genes <- sapply(allRes_salmon.CeC.sig$GO.ID, function(x){
#  genes <- genesInTerm(GoData_salmon.C, x)
#  genes[[1]][genes[[1]] %in% row.names(salmon_module)]
#})
#the above code does not work

genes <- genesInTerm(GoData_salmon.C, allRes_salmon.CeC.sig$GO.ID[1])
temp <- genes[[1]][genes[[1]] %in% row.names(salmon_module)]
allRes_salmon.CeC.sig$genes[1] <- list(temp)


allRes_salmon.genes.combined <- rbind(allRes_salmon.ByP.sig, allRes_salmon.MoF.sig, allRes_salmon.CeC.sig)

GoTerm_genes_salmon <- "IDs_salmon.xlsx"


for (i in 1:length(allRes_salmon.genes.combined$GO.ID)){
  genes_id <- allRes_salmon.genes.combined$genes[i] #get the genes from each term 
  Go_term <- allRes_salmon.genes.combined$GO.ID[i] #Get the Go term code 
  Go_name <- allRes_salmon.genes.combined$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- dammit_annot.modules[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(24,25,1:23)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/topGO_enrichement/", GoTerm_genes_salmon), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/topGO_enrichement/", GoTerm_genes_salmon), sheetName = str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}


write.csv(allRes_salmon.combined[,-7], file = paste0(supp.dir, "/module_2.tsv"), row.names = F)



#tan module ---- 
tan_module <- subset(dammit_annot.modules, dammit_annot.modules$module == "tan")

#GO.P
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.P %in% tan_module$gene_id))
names(compared_genes) <- bg_genes.P

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_tan.P <- new("topGOdata", ontology = "BP", 
                            description = "tan module enrichement analysis",
                            allGenes = compared_genes, 
                            annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                            gene2GO = background.P,
                            nodeSize = 5) 

#Run fisher test
resultFisher_tan.P <- runTest(GoData_tan.P, 
                                      algorithm = "weight01", 
                                      statistic = "fisher")
resultFisher_tan.P

allRes_tan.ByP <- GenTable(GoData_tan.P, 
                                   classicFisher =resultFisher_tan.P,
                                   orderBy= "classicFisher", 
                                   numChar = 100, 
                                   topNodes = length(resultFisher_tan.P@score))


allRes_tan.ByP.sig <- allRes_tan.ByP[which(as.numeric(allRes_tan.ByP$classicFisher) < 0.01),]
allRes_tan.ByP.sig <- allRes_tan.ByP.sig[order(allRes_tan.ByP.sig$Significant, decreasing = T), ]


#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% tan_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_tan.F <- new("topGOdata", ontology = "MF", 
                            description = "Black module enrichement analysis",
                            allGenes = compared_genes, 
                            annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                            gene2GO = background.F,
                            nodeSize = 5) 

#Run fisher test
resultFisher_tan.F <- runTest(GoData_tan.F, 
                                      algorithm = "weight01", 
                                      statistic = "fisher")
resultFisher_tan.F

allRes_tan.MoF <- GenTable(GoData_tan.F, 
                                   classicFisher =resultFisher_tan.F,
                                   orderBy= "classicFisher", 
                                   numChar = 100, 
                                   topNodes = length(resultFisher_tan.F@score))


allRes_tan.MoF.sig <- allRes_tan.MoF[which(as.numeric(allRes_tan.MoF$classicFisher) < 0.01),]
allRes_tan.MoF.sig <- allRes_tan.MoF.sig[order(allRes_tan.MoF.sig$Significant, decreasing = T), ]



#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% tan_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_tan.F <- new("topGOdata", ontology = "MF", 
                            description = "tan module enrichement analysis",
                            allGenes = compared_genes, 
                            annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                            gene2GO = background.F,
                            nodeSize = 5) 

#Run fisher test
resultFisher_tan.F <- runTest(GoData_tan.F, 
                                      algorithm = "weight01", 
                                      statistic = "fisher")
resultFisher_tan.F

allRes_tan.MoF <- GenTable(GoData_tan.F, 
                                   classicFisher =resultFisher_tan.F,
                                   orderBy= "classicFisher", 
                                   numChar = 100, 
                                   topNodes = length(resultFisher_tan.F@score))


allRes_tan.MoF.sig <- allRes_tan.MoF[which(as.numeric(allRes_tan.MoF$classicFisher) < 0.01),]
allRes_tan.MoF.sig <- allRes_tan.MoF.sig[order(allRes_tan.MoF.sig$Significant, decreasing = T), ]



#GO.C
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.C %in% tan_module$gene_id))
names(compared_genes) <- bg_genes.C

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_tan.C <- new("topGOdata", ontology = "CC", 
                            description = "tan module enrichement analysis",
                            allGenes = compared_genes, 
                            annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                            gene2GO = background.C,
                            nodeSize = 5) 

#Run fisher test
resultFisher_tan.C <- runTest(GoData_tan.C, 
                                      algorithm = "weight01", 
                                      statistic = "fisher")
resultFisher_tan.C

allRes_tan.CeC <- GenTable(GoData_tan.C, 
                                   classicFisher =resultFisher_tan.C,
                                   orderBy= "classicFisher", 
                                   numChar = 100, 
                                   topNodes = length(resultFisher_tan.C@score))


allRes_tan.CeC.sig <- allRes_tan.CeC[which(as.numeric(allRes_tan.CeC$classicFisher) < 0.01),]
allRes_tan.CeC.sig <- allRes_tan.CeC.sig[order(allRes_tan.CeC.sig$Significant, decreasing = T), ]


allRes_tan.combined <- rbind(allRes_tan.ByP.sig, allRes_tan.MoF.sig, allRes_tan.CeC.sig)
allRes_tan.combined$order <- seq(from = length(allRes_tan.combined$GO.ID), to = 1, by = -1)

tan_plot <- ggplot(allRes_tan.combined, aes(x = reorder(Term, order, decreasing = F), y = Significant)) + 
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "tan") + 
  coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.2) + 
  labs(title = "Module 7", 
       subtitle = "Enriched GO terms", 
       x = "", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))




png("results/topGO_enrichement/tanModule.png", res = 300, width = 4500, height = 3000)
tan_plot
dev.off()


#Make tables with the significant GO terms


allRes_tan.ByP.sig$genes <- sapply(allRes_tan.ByP.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_tan.P, x)
  genes[[1]][genes[[1]] %in% row.names(tan_module)]
})

allRes_tan.MoF.sig$genes <- sapply(allRes_tan.MoF.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_tan.F, x)
  genes[[1]][genes[[1]] %in% row.names(tan_module)]
})

allRes_tan.CeC.sig$genes <- sapply(allRes_tan.CeC.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_tan.C, x)
  genes[[1]][genes[[1]] %in% row.names(tan_module)]
})

allRes_tan.genes.combined <- rbind(allRes_tan.ByP.sig, allRes_tan.MoF.sig, allRes_tan.CeC.sig)

GoTerm_genes_tan <- "IDs_tan.xlsx"


for (i in 1:length(allRes_tan.genes.combined$GO.ID)){
  genes_id <- allRes_tan.genes.combined$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_tan.genes.combined$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- dammit_annot.modules[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(24,25,1:23)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/topGO_enrichement/", GoTerm_genes_tan), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/topGO_enrichement/", GoTerm_genes_tan), sheetName = str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}

write.csv(allRes_tan.combined[,-7], file = paste0(supp.dir, "/module_7.tsv"), row.names = F)



#skyblue module ---- 
skyblue_module <- subset(dammit_annot.modules, dammit_annot.modules$module == "skyblue")

#GO.P
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.P %in% skyblue_module$gene_id))
names(compared_genes) <- bg_genes.P

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_skyblue.P <- new("topGOdata", ontology = "BP", 
                    description = "skyblue module enrichement analysis",
                    allGenes = compared_genes, 
                    annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                    gene2GO = background.P,
                    nodeSize = 5) 

#Run fisher test
resultFisher_skyblue.P <- runTest(GoData_skyblue.P, 
                              algorithm = "weight01", 
                              statistic = "fisher")
resultFisher_skyblue.P

allRes_skyblue.ByP <- GenTable(GoData_skyblue.P, 
                           classicFisher =resultFisher_skyblue.P,
                           orderBy= "classicFisher", 
                           numChar = 100, 
                           topNodes = length(resultFisher_skyblue.P@score))


allRes_skyblue.ByP.sig <- allRes_skyblue.ByP[which(as.numeric(allRes_skyblue.ByP$classicFisher) < 0.01),]
allRes_skyblue.ByP.sig <- allRes_skyblue.ByP.sig[order(allRes_skyblue.ByP.sig$Significant, decreasing = T), ]


#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% skyblue_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_skyblue.F <- new("topGOdata", ontology = "MF", 
                    description = "Black module enrichement analysis",
                    allGenes = compared_genes, 
                    annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                    gene2GO = background.F,
                    nodeSize = 5) 

#Run fisher test
resultFisher_skyblue.F <- runTest(GoData_skyblue.F, 
                              algorithm = "weight01", 
                              statistic = "fisher")
resultFisher_skyblue.F

allRes_skyblue.MoF <- GenTable(GoData_skyblue.F, 
                           classicFisher =resultFisher_skyblue.F,
                           orderBy= "classicFisher", 
                           numChar = 100, 
                           topNodes = length(resultFisher_skyblue.F@score))


allRes_skyblue.MoF.sig <- allRes_skyblue.MoF[which(as.numeric(allRes_skyblue.MoF$classicFisher) < 0.01),]
allRes_skyblue.MoF.sig <- allRes_skyblue.MoF.sig[order(allRes_skyblue.MoF.sig$Significant, decreasing = T), ]



#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% skyblue_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_skyblue.F <- new("topGOdata", ontology = "MF", 
                    description = "skyblue module enrichement analysis",
                    allGenes = compared_genes, 
                    annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                    gene2GO = background.F,
                    nodeSize = 5) 

#Run fisher test
resultFisher_skyblue.F <- runTest(GoData_skyblue.F, 
                              algorithm = "weight01", 
                              statistic = "fisher")
resultFisher_skyblue.F

allRes_skyblue.MoF <- GenTable(GoData_skyblue.F, 
                           classicFisher =resultFisher_skyblue.F,
                           orderBy= "classicFisher", 
                           numChar = 100, 
                           topNodes = length(resultFisher_skyblue.F@score))


allRes_skyblue.MoF.sig <- allRes_skyblue.MoF[which(as.numeric(allRes_skyblue.MoF$classicFisher) < 0.01),]
allRes_skyblue.MoF.sig <- allRes_skyblue.MoF.sig[order(allRes_skyblue.MoF.sig$Significant, decreasing = T), ]



#GO.C
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.C %in% skyblue_module$gene_id))
names(compared_genes) <- bg_genes.C

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_skyblue.C <- new("topGOdata", ontology = "CC", 
                    description = "skyblue module enrichement analysis",
                    allGenes = compared_genes, 
                    annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                    gene2GO = background.C,
                    nodeSize = 5) 

#Run fisher test
resultFisher_skyblue.C <- runTest(GoData_skyblue.C, 
                              algorithm = "weight01", 
                              statistic = "fisher")
resultFisher_skyblue.C

allRes_skyblue.CeC <- GenTable(GoData_skyblue.C, 
                           classicFisher =resultFisher_skyblue.C,
                           orderBy= "classicFisher", 
                           numChar = 100, 
                           topNodes = length(resultFisher_skyblue.C@score))


allRes_skyblue.CeC.sig <- allRes_skyblue.CeC[which(as.numeric(allRes_skyblue.CeC$classicFisher) < 0.01),]
allRes_skyblue.CeC.sig <- allRes_skyblue.CeC.sig[order(allRes_skyblue.CeC.sig$Significant, decreasing = T), ]


allRes_skyblue.combined <- rbind(allRes_skyblue.ByP.sig, allRes_skyblue.MoF.sig, allRes_skyblue.CeC.sig)
allRes_skyblue.combined$order <- seq(from = length(allRes_skyblue.combined$GO.ID), to = 1, by = -1)

skyblue_plot <- ggplot(allRes_skyblue.combined, aes(x = reorder(Term, order, decreasing = F), y = Significant)) + 
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "skyblue") + 
  coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.2) + 
  labs(title = "Module 10", 
       subtitle = "Enriched GO terms", 
       x = "", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))




png("results/topGO_enrichement/skyblueModule.png", res = 300, width = 4500, height = 3000)
skyblue_plot
dev.off()


#Make tables with the significant GO terms


allRes_skyblue.ByP.sig$genes <- sapply(allRes_skyblue.ByP.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_skyblue.P, x)
  genes[[1]][genes[[1]] %in% row.names(skyblue_module)]
})

#ignoring the molecular function (does not seem to relevant in this case ...)
#allRes_skyblue.MoF.sig$genes <- sapply(allRes_skyblue.MoF.sig$GO.ID, function(x){
#  genes <- genesInTerm(GoData_skyblue.F, x)
#  genes[[1]][genes[[1]] %in% row.names(skyblue_module)]
#})

allRes_skyblue.CeC.sig$genes <- sapply(allRes_skyblue.CeC.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_skyblue.C, x)
  genes[[1]][genes[[1]] %in% row.names(skyblue_module)]
})

allRes_skyblue.genes.combined <- rbind(allRes_skyblue.ByP.sig, allRes_skyblue.CeC.sig)

GoTerm_genes_skyblue <- "IDs_skyblue.xlsx"


for (i in 1:length(allRes_skyblue.genes.combined$GO.ID)){
  genes_id <- allRes_skyblue.genes.combined$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_skyblue.genes.combined$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- dammit_annot.modules[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(24,25,1:23)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/topGO_enrichement/", GoTerm_genes_skyblue), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/topGO_enrichement/", GoTerm_genes_skyblue), sheetName = str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}

write.csv(allRes_skyblue.combined[,-7], file = paste0(supp.dir, "/module_10.tsv"), row.names = F)




#darkolivegreen module ---- 
darkolivegreen_module <- subset(dammit_annot.modules, dammit_annot.modules$module == "darkolivegreen")

#GO.P
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.P %in% darkolivegreen_module$gene_id))
names(compared_genes) <- bg_genes.P

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_darkolivegreen.P <- new("topGOdata", ontology = "BP", 
                          description = "darkolivegreen module enrichement analysis",
                          allGenes = compared_genes, 
                          annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                          gene2GO = background.P,
                          nodeSize = 5) 

#Run fisher test
resultFisher_darkolivegreen.P <- runTest(GoData_darkolivegreen.P, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_darkolivegreen.P

allRes_darkolivegreen.ByP <- GenTable(GoData_darkolivegreen.P, 
                                 classicFisher =resultFisher_darkolivegreen.P,
                                 orderBy= "classicFisher", 
                                 numChar = 100, 
                                 topNodes = length(resultFisher_darkolivegreen.P@score))


allRes_darkolivegreen.ByP.sig <- allRes_darkolivegreen.ByP[which(as.numeric(allRes_darkolivegreen.ByP$classicFisher) < 0.01),]
allRes_darkolivegreen.ByP.sig <- allRes_darkolivegreen.ByP.sig[order(allRes_darkolivegreen.ByP.sig$Significant, decreasing = T), ]


#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% darkolivegreen_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_darkolivegreen.F <- new("topGOdata", ontology = "MF", 
                          description = "Black module enrichement analysis",
                          allGenes = compared_genes, 
                          annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                          gene2GO = background.F,
                          nodeSize = 5) 

#Run fisher test
resultFisher_darkolivegreen.F <- runTest(GoData_darkolivegreen.F, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_darkolivegreen.F

allRes_darkolivegreen.MoF <- GenTable(GoData_darkolivegreen.F, 
                                 classicFisher =resultFisher_darkolivegreen.F,
                                 orderBy= "classicFisher", 
                                 numChar = 100, 
                                 topNodes = length(resultFisher_darkolivegreen.F@score))


allRes_darkolivegreen.MoF.sig <- allRes_darkolivegreen.MoF[which(as.numeric(allRes_darkolivegreen.MoF$classicFisher) < 0.01),]
allRes_darkolivegreen.MoF.sig <- allRes_darkolivegreen.MoF.sig[order(allRes_darkolivegreen.MoF.sig$Significant, decreasing = T), ]



#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% darkolivegreen_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_darkolivegreen.F <- new("topGOdata", ontology = "MF", 
                          description = "darkolivegreen module enrichement analysis",
                          allGenes = compared_genes, 
                          annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                          gene2GO = background.F,
                          nodeSize = 5) 

#Run fisher test
resultFisher_darkolivegreen.F <- runTest(GoData_darkolivegreen.F, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_darkolivegreen.F

allRes_darkolivegreen.MoF <- GenTable(GoData_darkolivegreen.F, 
                                 classicFisher =resultFisher_darkolivegreen.F,
                                 orderBy= "classicFisher", 
                                 numChar = 100, 
                                 topNodes = length(resultFisher_darkolivegreen.F@score))


allRes_darkolivegreen.MoF.sig <- allRes_darkolivegreen.MoF[which(as.numeric(allRes_darkolivegreen.MoF$classicFisher) < 0.01),]
allRes_darkolivegreen.MoF.sig <- allRes_darkolivegreen.MoF.sig[order(allRes_darkolivegreen.MoF.sig$Significant, decreasing = T), ]



#GO.C
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.C %in% darkolivegreen_module$gene_id))
names(compared_genes) <- bg_genes.C

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_darkolivegreen.C <- new("topGOdata", ontology = "CC", 
                          description = "darkolivegreen module enrichement analysis",
                          allGenes = compared_genes, 
                          annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                          gene2GO = background.C,
                          nodeSize = 5) 

#Run fisher test
resultFisher_darkolivegreen.C <- runTest(GoData_darkolivegreen.C, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_darkolivegreen.C

allRes_darkolivegreen.CeC <- GenTable(GoData_darkolivegreen.C, 
                                 classicFisher =resultFisher_darkolivegreen.C,
                                 orderBy= "classicFisher", 
                                 numChar = 100, 
                                 topNodes = length(resultFisher_darkolivegreen.C@score))


allRes_darkolivegreen.CeC.sig <- allRes_darkolivegreen.CeC[which(as.numeric(allRes_darkolivegreen.CeC$classicFisher) < 0.01),]
allRes_darkolivegreen.CeC.sig <- allRes_darkolivegreen.CeC.sig[order(allRes_darkolivegreen.CeC.sig$Significant, decreasing = T), ]


allRes_darkolivegreen.combined <- rbind(allRes_darkolivegreen.ByP.sig, allRes_darkolivegreen.MoF.sig, allRes_darkolivegreen.CeC.sig)
allRes_darkolivegreen.combined$order <- seq(from = length(allRes_darkolivegreen.combined$GO.ID), to = 1, by = -1)

darkolivegreen_plot <- ggplot(allRes_darkolivegreen.combined, aes(x = reorder(Term, order, decreasing = F), y = Significant)) + 
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "darkolivegreen") + 
  coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.2) + 
  labs(title = "Module 12", 
       subtitle = "Enriched GO terms", 
       x = "", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))




png("results/topGO_enrichement/darkolivegreenModule.png", res = 300, width = 4500, height = 3000)
darkolivegreen_plot
dev.off()




allRes_darkolivegreen.ByP.sig$genes <- sapply(allRes_darkolivegreen.ByP.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_darkolivegreen.P, x)
  genes[[1]][genes[[1]] %in% row.names(darkolivegreen_module)]
})

allRes_darkolivegreen.MoF.sig$genes <- sapply(allRes_darkolivegreen.MoF.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_darkolivegreen.F, x)
  genes[[1]][genes[[1]] %in% row.names(darkolivegreen_module)]
})

allRes_darkolivegreen.CeC.sig$genes <- sapply(allRes_darkolivegreen.CeC.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_darkolivegreen.C, x)
  genes[[1]][genes[[1]] %in% row.names(darkolivegreen_module)]
})

allRes_darkolivegreen.genes.combined <- rbind(allRes_darkolivegreen.ByP.sig, allRes_darkolivegreen.MoF.sig, allRes_darkolivegreen.CeC.sig)

GoTerm_genes_darkolivegreen <- "IDs_darkolivegreen.xlsx"


for (i in 1:length(allRes_darkolivegreen.genes.combined$GO.ID)){
  genes_id <- allRes_darkolivegreen.genes.combined$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_darkolivegreen.genes.combined$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- dammit_annot.modules[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(24,25,1:23)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/topGO_enrichement/", GoTerm_genes_darkolivegreen), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/topGO_enrichement/", GoTerm_genes_darkolivegreen), sheetName = str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}

write.csv(allRes_darkolivegreen.combined[,-7], file = paste0(supp.dir, "/module_12.tsv"), row.names = F)


#darkorange module ---- 
darkorange_module <- subset(dammit_annot.modules, dammit_annot.modules$module == "darkorange")

#GO.P
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.P %in% darkorange_module$gene_id))
names(compared_genes) <- bg_genes.P

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_darkorange.P <- new("topGOdata", ontology = "BP", 
                               description = "darkorange module enrichement analysis",
                               allGenes = compared_genes, 
                               annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                               gene2GO = background.P,
                               nodeSize = 5) 

#Run fisher test
resultFisher_darkorange.P <- runTest(GoData_darkorange.P, 
                                         algorithm = "weight01", 
                                         statistic = "fisher")
resultFisher_darkorange.P

allRes_darkorange.ByP <- GenTable(GoData_darkorange.P, 
                                      classicFisher =resultFisher_darkorange.P,
                                      orderBy= "classicFisher", 
                                      numChar = 100, 
                                      topNodes = length(resultFisher_darkorange.P@score))


allRes_darkorange.ByP.sig <- allRes_darkorange.ByP[which(as.numeric(allRes_darkorange.ByP$classicFisher) < 0.01),]
allRes_darkorange.ByP.sig <- allRes_darkorange.ByP.sig[order(allRes_darkorange.ByP.sig$Significant, decreasing = T), ]


#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% darkorange_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_darkorange.F <- new("topGOdata", ontology = "MF", 
                               description = "Black module enrichement analysis",
                               allGenes = compared_genes, 
                               annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                               gene2GO = background.F,
                               nodeSize = 5) 

#Run fisher test
resultFisher_darkorange.F <- runTest(GoData_darkorange.F, 
                                         algorithm = "weight01", 
                                         statistic = "fisher")
resultFisher_darkorange.F

allRes_darkorange.MoF <- GenTable(GoData_darkorange.F, 
                                      classicFisher =resultFisher_darkorange.F,
                                      orderBy= "classicFisher", 
                                      numChar = 100, 
                                      topNodes = length(resultFisher_darkorange.F@score))


allRes_darkorange.MoF.sig <- allRes_darkorange.MoF[which(as.numeric(allRes_darkorange.MoF$classicFisher) < 0.01),]
allRes_darkorange.MoF.sig <- allRes_darkorange.MoF.sig[order(allRes_darkorange.MoF.sig$Significant, decreasing = T), ]



#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% darkorange_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_darkorange.F <- new("topGOdata", ontology = "MF", 
                               description = "darkorange module enrichement analysis",
                               allGenes = compared_genes, 
                               annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                               gene2GO = background.F,
                               nodeSize = 5) 

#Run fisher test
resultFisher_darkorange.F <- runTest(GoData_darkorange.F, 
                                         algorithm = "weight01", 
                                         statistic = "fisher")
resultFisher_darkorange.F

allRes_darkorange.MoF <- GenTable(GoData_darkorange.F, 
                                      classicFisher =resultFisher_darkorange.F,
                                      orderBy= "classicFisher", 
                                      numChar = 100, 
                                      topNodes = length(resultFisher_darkorange.F@score))


allRes_darkorange.MoF.sig <- allRes_darkorange.MoF[which(as.numeric(allRes_darkorange.MoF$classicFisher) < 0.01),]
allRes_darkorange.MoF.sig <- allRes_darkorange.MoF.sig[order(allRes_darkorange.MoF.sig$Significant, decreasing = T), ]



#GO.C
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.C %in% darkorange_module$gene_id))
names(compared_genes) <- bg_genes.C

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_darkorange.C <- new("topGOdata", ontology = "CC", 
                               description = "darkorange module enrichement analysis",
                               allGenes = compared_genes, 
                               annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                               gene2GO = background.C,
                               nodeSize = 5) 

#Run fisher test
resultFisher_darkorange.C <- runTest(GoData_darkorange.C, 
                                         algorithm = "weight01", 
                                         statistic = "fisher")
resultFisher_darkorange.C

allRes_darkorange.CeC <- GenTable(GoData_darkorange.C, 
                                      classicFisher =resultFisher_darkorange.C,
                                      orderBy= "classicFisher", 
                                      numChar = 100, 
                                      topNodes = length(resultFisher_darkorange.C@score))


allRes_darkorange.CeC.sig <- allRes_darkorange.CeC[which(as.numeric(allRes_darkorange.CeC$classicFisher) < 0.01),]
allRes_darkorange.CeC.sig <- allRes_darkorange.CeC.sig[order(allRes_darkorange.CeC.sig$Significant, decreasing = T), ]


allRes_darkorange.combined <- rbind(allRes_darkorange.ByP.sig, allRes_darkorange.MoF.sig, allRes_darkorange.CeC.sig)
allRes_darkorange.combined$order <- seq(from = length(allRes_darkorange.combined$GO.ID), to = 1, by = -1)

darkorange_plot <- ggplot(allRes_darkorange.combined, aes(x = reorder(Term, order, decreasing = F), y = Significant)) + 
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "darkorange") + 
  coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.2) + 
  labs(title = "darkorangeModule", 
       subtitle = "Enriched GO terms", 
       x = "", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))




png("results/topGO_enrichement/darkorangeModule.png", res = 300, width = 4500, height = 3000)
darkorange_plot
dev.off()




#grey60 module ---- 
grey60_module <- subset(dammit_annot.modules, dammit_annot.modules$module == "grey60")

#GO.P
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.P %in% grey60_module$gene_id))
names(compared_genes) <- bg_genes.P

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_grey60.P <- new("topGOdata", ontology = "BP", 
                           description = "grey60 module enrichement analysis",
                           allGenes = compared_genes, 
                           annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                           gene2GO = background.P,
                           nodeSize = 5) 

#Run fisher test
resultFisher_grey60.P <- runTest(GoData_grey60.P, 
                                     algorithm = "weight01", 
                                     statistic = "fisher")
resultFisher_grey60.P

allRes_grey60.ByP <- GenTable(GoData_grey60.P, 
                                  classicFisher =resultFisher_grey60.P,
                                  orderBy= "classicFisher", 
                                  numChar = 100, 
                                  topNodes = length(resultFisher_grey60.P@score))


allRes_grey60.ByP.sig <- allRes_grey60.ByP[which(as.numeric(allRes_grey60.ByP$classicFisher) < 0.01),]
allRes_grey60.ByP.sig <- allRes_grey60.ByP.sig[order(allRes_grey60.ByP.sig$Significant, decreasing = T), ]


#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% grey60_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_grey60.F <- new("topGOdata", ontology = "MF", 
                           description = "Black module enrichement analysis",
                           allGenes = compared_genes, 
                           annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                           gene2GO = background.F,
                           nodeSize = 5) 

#Run fisher test
resultFisher_grey60.F <- runTest(GoData_grey60.F, 
                                     algorithm = "weight01", 
                                     statistic = "fisher")
resultFisher_grey60.F

allRes_grey60.MoF <- GenTable(GoData_grey60.F, 
                                  classicFisher =resultFisher_grey60.F,
                                  orderBy= "classicFisher", 
                                  numChar = 100, 
                                  topNodes = length(resultFisher_grey60.F@score))


allRes_grey60.MoF.sig <- allRes_grey60.MoF[which(as.numeric(allRes_grey60.MoF$classicFisher) < 0.01),]
allRes_grey60.MoF.sig <- allRes_grey60.MoF.sig[order(allRes_grey60.MoF.sig$Significant, decreasing = T), ]



#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% grey60_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_grey60.F <- new("topGOdata", ontology = "MF", 
                           description = "grey60 module enrichement analysis",
                           allGenes = compared_genes, 
                           annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                           gene2GO = background.F,
                           nodeSize = 5) 

#Run fisher test
resultFisher_grey60.F <- runTest(GoData_grey60.F, 
                                     algorithm = "weight01", 
                                     statistic = "fisher")
resultFisher_grey60.F

allRes_grey60.MoF <- GenTable(GoData_grey60.F, 
                                  classicFisher =resultFisher_grey60.F,
                                  orderBy= "classicFisher", 
                                  numChar = 100, 
                                  topNodes = length(resultFisher_grey60.F@score))


allRes_grey60.MoF.sig <- allRes_grey60.MoF[which(as.numeric(allRes_grey60.MoF$classicFisher) < 0.01),]
allRes_grey60.MoF.sig <- allRes_grey60.MoF.sig[order(allRes_grey60.MoF.sig$Significant, decreasing = T), ]



#GO.C
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.C %in% grey60_module$gene_id))
names(compared_genes) <- bg_genes.C

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_grey60.C <- new("topGOdata", ontology = "CC", 
                           description = "grey60 module enrichement analysis",
                           allGenes = compared_genes, 
                           annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                           gene2GO = background.C,
                           nodeSize = 5) 

#Run fisher test
resultFisher_grey60.C <- runTest(GoData_grey60.C, 
                                     algorithm = "weight01", 
                                     statistic = "fisher")
resultFisher_grey60.C

allRes_grey60.CeC <- GenTable(GoData_grey60.C, 
                                  classicFisher =resultFisher_grey60.C,
                                  orderBy= "classicFisher", 
                                  numChar = 100, 
                                  topNodes = length(resultFisher_grey60.C@score))


allRes_grey60.CeC.sig <- allRes_grey60.CeC[which(as.numeric(allRes_grey60.CeC$classicFisher) < 0.01),]
allRes_grey60.CeC.sig <- allRes_grey60.CeC.sig[order(allRes_grey60.CeC.sig$Significant, decreasing = T), ]


allRes_grey60.combined <- rbind(allRes_grey60.ByP.sig, allRes_grey60.MoF.sig, allRes_grey60.CeC.sig)
allRes_grey60.combined$order <- seq(from = length(allRes_grey60.combined$GO.ID), to = 1, by = -1)

grey60_plot <- ggplot(allRes_grey60.combined, aes(x = reorder(Term, order, decreasing = F), y = Significant)) + 
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "grey60") + 
  coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.2) + 
  labs(title = "grey60Module", 
       subtitle = "Enriched GO terms", 
       x = "", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))




png("results/topGO_enrichement/grey60Module.png", res = 300, width = 4500, height = 3000)
grey60_plot
dev.off()




#magenta module ---- 
magenta_module <- subset(dammit_annot.modules, dammit_annot.modules$module == "magenta")

#GO.P
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.P %in% magenta_module$gene_id))
names(compared_genes) <- bg_genes.P

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_magenta.P <- new("topGOdata", ontology = "BP", 
                           description = "magenta module enrichement analysis",
                           allGenes = compared_genes, 
                           annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                           gene2GO = background.P,
                           nodeSize = 5) 

#Run fisher test
resultFisher_magenta.P <- runTest(GoData_magenta.P, 
                                     algorithm = "weight01", 
                                     statistic = "fisher")
resultFisher_magenta.P

allRes_magenta.ByP <- GenTable(GoData_magenta.P, 
                                  classicFisher =resultFisher_magenta.P,
                                  orderBy= "classicFisher", 
                                  numChar = 100, 
                                  topNodes = length(resultFisher_magenta.P@score))


allRes_magenta.ByP.sig <- allRes_magenta.ByP[which(as.numeric(allRes_magenta.ByP$classicFisher) < 0.01),]
allRes_magenta.ByP.sig <- allRes_magenta.ByP.sig[order(allRes_magenta.ByP.sig$Significant, decreasing = T), ]


#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% magenta_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_magenta.F <- new("topGOdata", ontology = "MF", 
                           description = "Black module enrichement analysis",
                           allGenes = compared_genes, 
                           annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                           gene2GO = background.F,
                           nodeSize = 5) 

#Run fisher test
resultFisher_magenta.F <- runTest(GoData_magenta.F, 
                                     algorithm = "weight01", 
                                     statistic = "fisher")
resultFisher_magenta.F

allRes_magenta.MoF <- GenTable(GoData_magenta.F, 
                                  classicFisher =resultFisher_magenta.F,
                                  orderBy= "classicFisher", 
                                  numChar = 100, 
                                  topNodes = length(resultFisher_magenta.F@score))


allRes_magenta.MoF.sig <- allRes_magenta.MoF[which(as.numeric(allRes_magenta.MoF$classicFisher) < 0.01),]
allRes_magenta.MoF.sig <- allRes_magenta.MoF.sig[order(allRes_magenta.MoF.sig$Significant, decreasing = T), ]



#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% magenta_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_magenta.F <- new("topGOdata", ontology = "MF", 
                           description = "magenta module enrichement analysis",
                           allGenes = compared_genes, 
                           annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                           gene2GO = background.F,
                           nodeSize = 5) 

#Run fisher test
resultFisher_magenta.F <- runTest(GoData_magenta.F, 
                                     algorithm = "weight01", 
                                     statistic = "fisher")
resultFisher_magenta.F

allRes_magenta.MoF <- GenTable(GoData_magenta.F, 
                                  classicFisher =resultFisher_magenta.F,
                                  orderBy= "classicFisher", 
                                  numChar = 100, 
                                  topNodes = length(resultFisher_magenta.F@score))


allRes_magenta.MoF.sig <- allRes_magenta.MoF[which(as.numeric(allRes_magenta.MoF$classicFisher) < 0.01),]
allRes_magenta.MoF.sig <- allRes_magenta.MoF.sig[order(allRes_magenta.MoF.sig$Significant, decreasing = T), ]



#GO.C
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.C %in% magenta_module$gene_id))
names(compared_genes) <- bg_genes.C

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_magenta.C <- new("topGOdata", ontology = "CC", 
                           description = "magenta module enrichement analysis",
                           allGenes = compared_genes, 
                           annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                           gene2GO = background.C,
                           nodeSize = 5) 

#Run fisher test
resultFisher_magenta.C <- runTest(GoData_magenta.C, 
                                     algorithm = "weight01", 
                                     statistic = "fisher")
resultFisher_magenta.C

allRes_magenta.CeC <- GenTable(GoData_magenta.C, 
                                  classicFisher =resultFisher_magenta.C,
                                  orderBy= "classicFisher", 
                                  numChar = 100, 
                                  topNodes = length(resultFisher_magenta.C@score))


allRes_magenta.CeC.sig <- allRes_magenta.CeC[which(as.numeric(allRes_magenta.CeC$classicFisher) < 0.01),]
allRes_magenta.CeC.sig <- allRes_magenta.CeC.sig[order(allRes_magenta.CeC.sig$Significant, decreasing = T), ]


allRes_magenta.combined <- rbind(allRes_magenta.ByP.sig, allRes_magenta.MoF.sig, allRes_magenta.CeC.sig)
allRes_magenta.combined$order <- seq(from = length(allRes_magenta.combined$GO.ID), to = 1, by = -1)

magenta_plot <- ggplot(allRes_magenta.combined, aes(x = reorder(Term, order, decreasing = F), y = Significant)) + 
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "magenta") + 
  coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.2) + 
  labs(title = "magentaModule", 
       subtitle = "Enriched GO terms", 
       x = "", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))




png("results/topGO_enrichement/magentaModule.png", res = 300, width = 4500, height = 3000)
magenta_plot
dev.off()



#midnightblue module ---- 
midnightblue_module <- subset(dammit_annot.modules, dammit_annot.modules$module == "midnightblue")

#GO.P
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.P %in% midnightblue_module$gene_id))
names(compared_genes) <- bg_genes.P

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_midnightblue.P <- new("topGOdata", ontology = "BP", 
                           description = "midnightblue module enrichement analysis",
                           allGenes = compared_genes, 
                           annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                           gene2GO = background.P,
                           nodeSize = 5) 

#Run fisher test
resultFisher_midnightblue.P <- runTest(GoData_midnightblue.P, 
                                     algorithm = "weight01", 
                                     statistic = "fisher")
resultFisher_midnightblue.P

allRes_midnightblue.ByP <- GenTable(GoData_midnightblue.P, 
                                  classicFisher =resultFisher_midnightblue.P,
                                  orderBy= "classicFisher", 
                                  numChar = 100, 
                                  topNodes = length(resultFisher_midnightblue.P@score))


allRes_midnightblue.ByP.sig <- allRes_midnightblue.ByP[which(as.numeric(allRes_midnightblue.ByP$classicFisher) < 0.01),]
allRes_midnightblue.ByP.sig <- allRes_midnightblue.ByP.sig[order(allRes_midnightblue.ByP.sig$Significant, decreasing = T), ]


#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% midnightblue_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_midnightblue.F <- new("topGOdata", ontology = "MF", 
                           description = "Black module enrichement analysis",
                           allGenes = compared_genes, 
                           annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                           gene2GO = background.F,
                           nodeSize = 5) 

#Run fisher test
resultFisher_midnightblue.F <- runTest(GoData_midnightblue.F, 
                                     algorithm = "weight01", 
                                     statistic = "fisher")
resultFisher_midnightblue.F

allRes_midnightblue.MoF <- GenTable(GoData_midnightblue.F, 
                                  classicFisher =resultFisher_midnightblue.F,
                                  orderBy= "classicFisher", 
                                  numChar = 100, 
                                  topNodes = length(resultFisher_midnightblue.F@score))


allRes_midnightblue.MoF.sig <- allRes_midnightblue.MoF[which(as.numeric(allRes_midnightblue.MoF$classicFisher) < 0.01),]
allRes_midnightblue.MoF.sig <- allRes_midnightblue.MoF.sig[order(allRes_midnightblue.MoF.sig$Significant, decreasing = T), ]



#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% midnightblue_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_midnightblue.F <- new("topGOdata", ontology = "MF", 
                           description = "midnightblue module enrichement analysis",
                           allGenes = compared_genes, 
                           annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                           gene2GO = background.F,
                           nodeSize = 5) 

#Run fisher test
resultFisher_midnightblue.F <- runTest(GoData_midnightblue.F, 
                                     algorithm = "weight01", 
                                     statistic = "fisher")
resultFisher_midnightblue.F

allRes_midnightblue.MoF <- GenTable(GoData_midnightblue.F, 
                                  classicFisher =resultFisher_midnightblue.F,
                                  orderBy= "classicFisher", 
                                  numChar = 100, 
                                  topNodes = length(resultFisher_midnightblue.F@score))


allRes_midnightblue.MoF.sig <- allRes_midnightblue.MoF[which(as.numeric(allRes_midnightblue.MoF$classicFisher) < 0.01),]
allRes_midnightblue.MoF.sig <- allRes_midnightblue.MoF.sig[order(allRes_midnightblue.MoF.sig$Significant, decreasing = T), ]



#GO.C
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.C %in% midnightblue_module$gene_id))
names(compared_genes) <- bg_genes.C

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_midnightblue.C <- new("topGOdata", ontology = "CC", 
                           description = "midnightblue module enrichement analysis",
                           allGenes = compared_genes, 
                           annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                           gene2GO = background.C,
                           nodeSize = 5) 

#Run fisher test
resultFisher_midnightblue.C <- runTest(GoData_midnightblue.C, 
                                     algorithm = "weight01", 
                                     statistic = "fisher")
resultFisher_midnightblue.C

allRes_midnightblue.CeC <- GenTable(GoData_midnightblue.C, 
                                  classicFisher =resultFisher_midnightblue.C,
                                  orderBy= "classicFisher", 
                                  numChar = 100, 
                                  topNodes = length(resultFisher_midnightblue.C@score))


allRes_midnightblue.CeC.sig <- allRes_midnightblue.CeC[which(as.numeric(allRes_midnightblue.CeC$classicFisher) < 0.01),]
allRes_midnightblue.CeC.sig <- allRes_midnightblue.CeC.sig[order(allRes_midnightblue.CeC.sig$Significant, decreasing = T), ]


allRes_midnightblue.combined <- rbind(allRes_midnightblue.ByP.sig, allRes_midnightblue.MoF.sig, allRes_midnightblue.CeC.sig)
allRes_midnightblue.combined$order <- seq(from = length(allRes_midnightblue.combined$GO.ID), to = 1, by = -1)

midnightblue_plot <- ggplot(allRes_midnightblue.combined, aes(x = reorder(Term, order, decreasing = F), y = Significant)) + 
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "midnightblue") + 
  coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.2) + 
  labs(title = "midnightblueModule", 
       subtitle = "Enriched GO terms", 
       x = "", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))




png("results/topGO_enrichement/midnightblueModule.png", res = 300, width = 4500, height = 3000)
midnightblue_plot
dev.off()


#paleturquoise module ---- 
paleturquoise_module <- subset(dammit_annot.modules, dammit_annot.modules$module == "paleturquoise")

#GO.P
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.P %in% paleturquoise_module$gene_id))
names(compared_genes) <- bg_genes.P

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_paleturquoise.P <- new("topGOdata", ontology = "BP", 
                           description = "paleturquoise module enrichement analysis",
                           allGenes = compared_genes, 
                           annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                           gene2GO = background.P,
                           nodeSize = 5) 

#Run fisher test
resultFisher_paleturquoise.P <- runTest(GoData_paleturquoise.P, 
                                     algorithm = "weight01", 
                                     statistic = "fisher")
resultFisher_paleturquoise.P

allRes_paleturquoise.ByP <- GenTable(GoData_paleturquoise.P, 
                                  classicFisher =resultFisher_paleturquoise.P,
                                  orderBy= "classicFisher", 
                                  numChar = 100, 
                                  topNodes = length(resultFisher_paleturquoise.P@score))


allRes_paleturquoise.ByP.sig <- allRes_paleturquoise.ByP[which(as.numeric(allRes_paleturquoise.ByP$classicFisher) < 0.01),]
allRes_paleturquoise.ByP.sig <- allRes_paleturquoise.ByP.sig[order(allRes_paleturquoise.ByP.sig$Significant, decreasing = T), ]


#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% paleturquoise_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_paleturquoise.F <- new("topGOdata", ontology = "MF", 
                           description = "Black module enrichement analysis",
                           allGenes = compared_genes, 
                           annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                           gene2GO = background.F,
                           nodeSize = 5) 

#Run fisher test
resultFisher_paleturquoise.F <- runTest(GoData_paleturquoise.F, 
                                     algorithm = "weight01", 
                                     statistic = "fisher")
resultFisher_paleturquoise.F

allRes_paleturquoise.MoF <- GenTable(GoData_paleturquoise.F, 
                                  classicFisher =resultFisher_paleturquoise.F,
                                  orderBy= "classicFisher", 
                                  numChar = 100, 
                                  topNodes = length(resultFisher_paleturquoise.F@score))


allRes_paleturquoise.MoF.sig <- allRes_paleturquoise.MoF[which(as.numeric(allRes_paleturquoise.MoF$classicFisher) < 0.01),]
allRes_paleturquoise.MoF.sig <- allRes_paleturquoise.MoF.sig[order(allRes_paleturquoise.MoF.sig$Significant, decreasing = T), ]



#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% paleturquoise_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_paleturquoise.F <- new("topGOdata", ontology = "MF", 
                           description = "paleturquoise module enrichement analysis",
                           allGenes = compared_genes, 
                           annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                           gene2GO = background.F,
                           nodeSize = 5) 

#Run fisher test
resultFisher_paleturquoise.F <- runTest(GoData_paleturquoise.F, 
                                     algorithm = "weight01", 
                                     statistic = "fisher")
resultFisher_paleturquoise.F

allRes_paleturquoise.MoF <- GenTable(GoData_paleturquoise.F, 
                                  classicFisher =resultFisher_paleturquoise.F,
                                  orderBy= "classicFisher", 
                                  numChar = 100, 
                                  topNodes = length(resultFisher_paleturquoise.F@score))


allRes_paleturquoise.MoF.sig <- allRes_paleturquoise.MoF[which(as.numeric(allRes_paleturquoise.MoF$classicFisher) < 0.01),]
allRes_paleturquoise.MoF.sig <- allRes_paleturquoise.MoF.sig[order(allRes_paleturquoise.MoF.sig$Significant, decreasing = T), ]



#GO.C
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.C %in% paleturquoise_module$gene_id))
names(compared_genes) <- bg_genes.C

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_paleturquoise.C <- new("topGOdata", ontology = "CC", 
                           description = "paleturquoise module enrichement analysis",
                           allGenes = compared_genes, 
                           annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                           gene2GO = background.C,
                           nodeSize = 5) 

#Run fisher test
resultFisher_paleturquoise.C <- runTest(GoData_paleturquoise.C, 
                                     algorithm = "weight01", 
                                     statistic = "fisher")
resultFisher_paleturquoise.C

allRes_paleturquoise.CeC <- GenTable(GoData_paleturquoise.C, 
                                  classicFisher =resultFisher_paleturquoise.C,
                                  orderBy= "classicFisher", 
                                  numChar = 100, 
                                  topNodes = length(resultFisher_paleturquoise.C@score))


allRes_paleturquoise.CeC.sig <- allRes_paleturquoise.CeC[which(as.numeric(allRes_paleturquoise.CeC$classicFisher) < 0.01),]
allRes_paleturquoise.CeC.sig <- allRes_paleturquoise.CeC.sig[order(allRes_paleturquoise.CeC.sig$Significant, decreasing = T), ]


allRes_paleturquoise.combined <- rbind(allRes_paleturquoise.ByP.sig, allRes_paleturquoise.MoF.sig, allRes_paleturquoise.CeC.sig)
allRes_paleturquoise.combined$order <- seq(from = length(allRes_paleturquoise.combined$GO.ID), to = 1, by = -1)

paleturquoise_plot <- ggplot(allRes_paleturquoise.combined, aes(x = reorder(Term, order, decreasing = F), y = Significant)) + 
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "paleturquoise") + 
  coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.2) + 
  labs(title = "paleturquoiseModule", 
       subtitle = "Enriched GO terms", 
       x = "", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))




png("results/topGO_enrichement/paleturquoiseModule.png", res = 300, width = 4500, height = 3000)
paleturquoise_plot
dev.off()



#darkred module ---- 
darkred_module <- subset(dammit_annot.modules, dammit_annot.modules$module == "darkred")

#GO.P
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.P %in% darkred_module$gene_id))
names(compared_genes) <- bg_genes.P

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_darkred.P <- new("topGOdata", ontology = "BP", 
                              description = "darkred module enrichement analysis",
                              allGenes = compared_genes, 
                              annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                              gene2GO = background.P,
                              nodeSize = 5) 

#Run fisher test
resultFisher_darkred.P <- runTest(GoData_darkred.P, 
                                        algorithm = "weight01", 
                                        statistic = "fisher")
resultFisher_darkred.P

allRes_darkred.ByP <- GenTable(GoData_darkred.P, 
                                     classicFisher =resultFisher_darkred.P,
                                     orderBy= "classicFisher", 
                                     numChar = 100, 
                                     topNodes = length(resultFisher_darkred.P@score))


allRes_darkred.ByP.sig <- allRes_darkred.ByP[which(as.numeric(allRes_darkred.ByP$classicFisher) < 0.01),]
allRes_darkred.ByP.sig <- allRes_darkred.ByP.sig[order(allRes_darkred.ByP.sig$Significant, decreasing = T), ]


#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% darkred_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_darkred.F <- new("topGOdata", ontology = "MF", 
                              description = "Black module enrichement analysis",
                              allGenes = compared_genes, 
                              annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                              gene2GO = background.F,
                              nodeSize = 5) 

#Run fisher test
resultFisher_darkred.F <- runTest(GoData_darkred.F, 
                                        algorithm = "weight01", 
                                        statistic = "fisher")
resultFisher_darkred.F

allRes_darkred.MoF <- GenTable(GoData_darkred.F, 
                                     classicFisher =resultFisher_darkred.F,
                                     orderBy= "classicFisher", 
                                     numChar = 100, 
                                     topNodes = length(resultFisher_darkred.F@score))


allRes_darkred.MoF.sig <- allRes_darkred.MoF[which(as.numeric(allRes_darkred.MoF$classicFisher) < 0.01),]
allRes_darkred.MoF.sig <- allRes_darkred.MoF.sig[order(allRes_darkred.MoF.sig$Significant, decreasing = T), ]



#GO.F
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.F %in% darkred_module$gene_id))
names(compared_genes) <- bg_genes.F

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_darkred.F <- new("topGOdata", ontology = "MF", 
                              description = "darkred module enrichement analysis",
                              allGenes = compared_genes, 
                              annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                              gene2GO = background.F,
                              nodeSize = 5) 

#Run fisher test
resultFisher_darkred.F <- runTest(GoData_darkred.F, 
                                        algorithm = "weight01", 
                                        statistic = "fisher")
resultFisher_darkred.F

allRes_darkred.MoF <- GenTable(GoData_darkred.F, 
                                     classicFisher =resultFisher_darkred.F,
                                     orderBy= "classicFisher", 
                                     numChar = 100, 
                                     topNodes = length(resultFisher_darkred.F@score))


allRes_darkred.MoF.sig <- allRes_darkred.MoF[which(as.numeric(allRes_darkred.MoF$classicFisher) < 0.01),]
allRes_darkred.MoF.sig <- allRes_darkred.MoF.sig[order(allRes_darkred.MoF.sig$Significant, decreasing = T), ]



#GO.C
#get interesting genes 
compared_genes <- factor(as.integer(bg_genes.C %in% darkred_module$gene_id))
names(compared_genes) <- bg_genes.C

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_darkred.C <- new("topGOdata", ontology = "CC", 
                              description = "darkred module enrichement analysis",
                              allGenes = compared_genes, 
                              annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                              gene2GO = background.C,
                              nodeSize = 5) 

#Run fisher test
resultFisher_darkred.C <- runTest(GoData_darkred.C, 
                                        algorithm = "weight01", 
                                        statistic = "fisher")
resultFisher_darkred.C

allRes_darkred.CeC <- GenTable(GoData_darkred.C, 
                                     classicFisher =resultFisher_darkred.C,
                                     orderBy= "classicFisher", 
                                     numChar = 100, 
                                     topNodes = length(resultFisher_darkred.C@score))


allRes_darkred.CeC.sig <- allRes_darkred.CeC[which(as.numeric(allRes_darkred.CeC$classicFisher) < 0.01),]
allRes_darkred.CeC.sig <- allRes_darkred.CeC.sig[order(allRes_darkred.CeC.sig$Significant, decreasing = T), ]


allRes_darkred.combined <- rbind(allRes_darkred.ByP.sig, allRes_darkred.MoF.sig, allRes_darkred.CeC.sig)
allRes_darkred.combined$order <- seq(from = length(allRes_darkred.combined$GO.ID), to = 1, by = -1)

darkred_plot <- ggplot(allRes_darkred.combined, aes(x = reorder(Term, order, decreasing = F), y = Significant)) + 
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "darkred") + 
  coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.2) + 
  labs(title = "Module 9", 
       subtitle = "Enriched GO terms", 
       x = "", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))




png("results/topGO_enrichement/darkredModule.png", res = 300, width = 4500, height = 3000)
darkred_plot
dev.off()




allRes_darkred.ByP.sig$genes <- sapply(allRes_darkred.ByP.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_darkred.P, x)
  genes[[1]][genes[[1]] %in% row.names(darkred_module)]
})

allRes_darkred.MoF.sig$genes <- sapply(allRes_darkred.MoF.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_darkred.F, x)
  genes[[1]][genes[[1]] %in% row.names(darkred_module)]
})

allRes_darkred.CeC.sig$genes <- sapply(allRes_darkred.CeC.sig$GO.ID, function(x){
  genes <- genesInTerm(GoData_darkred.C, x)
  genes[[1]][genes[[1]] %in% row.names(darkred_module)]
})

allRes_darkred.genes.combined <- rbind(allRes_darkred.ByP.sig, allRes_darkred.MoF.sig, allRes_darkred.CeC.sig)

GoTerm_genes_darkred <- "IDs_darkred.xlsx"


for (i in 1:length(allRes_darkred.genes.combined$GO.ID)){
  genes_id <- allRes_darkred.genes.combined$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_darkred.genes.combined$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- dammit_annot.modules[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(24,25,1:23)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/topGO_enrichement/", GoTerm_genes_darkred), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/topGO_enrichement/", GoTerm_genes_darkred), sheetName = str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}

write.csv(allRes_darkred.combined[,-7], file = paste0(supp.dir, "/module_9.tsv"), row.names = F)
