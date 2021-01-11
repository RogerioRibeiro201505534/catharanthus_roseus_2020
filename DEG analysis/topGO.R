#Misc
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())
options(stringsAsFactors = F)

#Load library
library(topGO)
library(ggplot2)
library(stringr)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)
library(DBI)
library(xlsx)

#Load annotations files 
annotation_file <- read.csv("dammit_annotation_final.tsv", sep = "\t", row.names = 1)

#create results dir 
results.dir <- "results"
if (!dir.exists(results.dir)){
  dir.create(results.dir)
}

###########################################
####meso vs idio enrichement analysis #####
###########################################

#make results directory
results.dir_2 <- paste(results.dir, "meso_vs_idio", sep = "/")
if (!dir.exists(results.dir_2)){
  dir.create(results.dir_2)
}

results.dir_2.1 <- paste(results.dir_2, "FDR_0.1", sep = "/")
if (!dir.exists(results.dir_2.1)){
  dir.create(results.dir_2.1)
}

results.dir_2.05 <- paste(results.dir_2, "FDR_0.05", sep = "/")
if (!dir.exists(results.dir_2.05)){
  dir.create(results.dir_2.05)
}


results.dir_noFDR <- paste(results.dir_2, "Before_muliple_testing", sep = "/")
if (!dir.exists(results.dir_noFDR)){
  dir.create(results.dir_noFDR)
}

meso_vs_idio_bg_genes <- read.csv(file = "DESeq2_tables/meso_vs_idio_Res.csv")$X
meso_vs_idio_sub_annot <- annotation_file[meso_vs_idio_bg_genes,]


###Make backgroud annotation ---- 
#For each compariton, the background annotation is different as we only consider the gene with
#gene counts >0 for the compariton. For this filtering, the genes in the res table from DesEq will be used, 
#as genes with rowcounts = 0 were filtered out before 


##Subset annotation files for backgroud annotation 
##Go(P) <- biological process
idio_Go.P <- data.frame("gene" = row.names(meso_vs_idio_sub_annot), "GO.BP" = meso_vs_idio_sub_annot[,9])
idio_Go.P$'GO.BP' <- str_replace_all(idio_Go.P$'GO.BP', ";", ",")
write.table(idio_Go.P, file = "idio_GoP.txt", row.names = FALSE, col.names = F,quote = F, sep = "\t")

##Go(F) <- molecular function
idio_Go.F <- data.frame("gene" = row.names(meso_vs_idio_sub_annot), "GO.MF" = meso_vs_idio_sub_annot[,11])
idio_Go.F$'GO.MF' <- str_replace_all(idio_Go.F$'GO.MF', ";", ",")
write.table(idio_Go.F, file = "idio_GoF.txt", row.names = FALSE, col.names = F,quote = F, sep = "\t")

##Go(C) <- celular compartement
idio_Go.C <- data.frame("gene" = row.names(meso_vs_idio_sub_annot), "GO.CC" = meso_vs_idio_sub_annot[,13])
idio_Go.C$'GO.CC' <- str_replace_all(idio_Go.C$'GO.CC', ";", ",")
write.table(idio_Go.C, file = "idio_GoC.txt", row.names = FALSE, col.names = F,quote = F, sep = "\t")


##By Biological process ----
#Get background annotation 
idio_GOByP <- readMappings(file = "idio_GoP.txt")
bg_genes <- names(idio_GOByP)


###Upregulated genes ----
#get gene IDs for the enrichement 
genes_UP_meso_vs_idio <- read.csv("DESeq2_tables/meso_vs_idio_sigUP.csv", header = TRUE)[,2]

#Get annotation for upregulated genes 
UP.annotated <- meso_vs_idio_sub_annot[genes_UP_meso_vs_idio,]

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_UP_meso_vs_idio))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_meso_idio_up.P <- new("topGOdata", ontology = "BP", 
                  description = "meso_vs_idio enrichment using kallisto quantification and dammit annotation",
                  allGenes = compared_genes, 
                  annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                  gene2GO = idio_GOByP,
                  nodeSize = 5) 

#Run fisher test
resultFisher_idio_up.P <- runTest(GoData_meso_idio_up.P, 
                        algorithm = "weight01", 
                        statistic = "fisher")
resultFisher_idio_up.P

#4467 GO terms scored, 20 with p < 0.01
#(801 out of 1328 significant genes, and 13406 feasible genes, or gene universe)

allRes_meso_vs_idio_up.ByP <- GenTable(GoData_meso_idio_up.P, 
                                    classicFisher = resultFisher_idio_up.P,
                                    orderBy= "classicFisher", 
                                    numChar = 100, 
                                    topNodes = length(resultFisher_idio_up.P@score))

#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the upregulated genes 
#were removed 
allRes_meso_vs_idio_up.ByP_filtered <- allRes_meso_vs_idio_up.ByP[which(allRes_meso_vs_idio_up.ByP$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_meso_vs_idio_up.ByP_filtered$classicFisher, method = "fdr")

allRes_meso_vs_idio_up.ByP_filtered$padj <- adj.pvalues

write.xlsx(allRes_meso_vs_idio_up.ByP_filtered, file = paste(results.dir_2,"byP_UP_resTable.xlsx", sep = "/"), row.names = F) 

#Make plot for 0.1 and 0.05 FDR levels 
n_term_sig <- sum(allRes_meso_vs_idio_up.ByP_filtered$padj <= 0.1)

#Plot the sig Go terms 
idio_BP.up <- ggplot(allRes_meso_vs_idio_up.ByP_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "Meso vs idio upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.1, "BP_UP_FDR_0.1.png", sep = "/"), width = 1000, height = 600)
idio_BP.up
dev.off()


n_term_sig_2 <- sum(allRes_meso_vs_idio_up.ByP_filtered$padj <= 0.05)

#Plot the sig Go terms 
idio_BP.up <- ggplot(allRes_meso_vs_idio_up.ByP_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "Meso vs idio upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.05, "BP_UP_FDR_0.05.png", sep = "/"), width = 1000, height = 600)
idio_BP.up
dev.off()


allRes_meso_vs_idio.ByP_ordered <- allRes_meso_vs_idio_up.ByP_filtered[which(as.numeric(allRes_meso_vs_idio_up.ByP_filtered$classicFisher) < 0.01),]
allRes_meso_vs_idio.ByP_ordered <- allRes_meso_vs_idio.ByP_ordered[order(allRes_meso_vs_idio.ByP_ordered$Significant, decreasing = T),]

allRes_meso_vs_idio.ByP_ordered[15,2] <- "regulation of Cys-type endopeptidase activity involved in apoptotic process"

idio_BP_2.up <- ggplot(allRes_meso_vs_idio.ByP_ordered[c(1:20),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Mesophyll vs idioblast upregulated", 
       subtitle = "Biological process", y = "Significant genes number", x = "") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_noFDR, "BP_UP_top20.png", sep = "/"), width = 4500, height = 3000, res = 300)
idio_BP_2.up
dev.off()


###Downregulated genes ----

#get gene IDs for the enrichement 
genes_Down_meso_vs_idio <- read.csv("DESeq2_tables/meso_vs_idio_sigDOWN.csv", header = TRUE)[,2]

#Get annotation for upregulated genes 
Down.annotated <- meso_vs_idio_sub_annot[genes_Down_meso_vs_idio,]

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_Down_meso_vs_idio))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_meso_idio_Down.P <- new("topGOdata", ontology = "BP", 
                             description = "meso_vs_idio enrichment using kallisto quantification and dammit annotation",
                             allGenes = compared_genes, 
                             annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                             gene2GO = idio_GOByP,
                             nodeSize = 5) 

#Run fisher test
resultFisher_idio_Down.P <- runTest(GoData_meso_idio_Down.P, 
                                  algorithm = "weight01", 
                                  statistic = "fisher")
resultFisher_idio_Down.P

#4498 GO terms scored, 41 with p < 0.01
#(788 out of 1308 significant genes, and 13406 feasible genes, or gene universe)

allRes_meso_vs_idio_Down.ByP <- GenTable(GoData_meso_idio_Down.P, 
                                    classicFisher = resultFisher_idio_Down.P,
                                    orderBy= "classicFisher", 
                                    numChar = 100, 
                                    topNodes = length(resultFisher_idio_Down.P@score))


#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the upregulated genes 
#were removed 
allRes_meso_vs_idio_Down.ByP_filtered <- allRes_meso_vs_idio_Down.ByP[which(allRes_meso_vs_idio_Down.ByP$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_meso_vs_idio_Down.ByP_filtered$classicFisher, method = "fdr")

allRes_meso_vs_idio_Down.ByP_filtered$padj <- adj.pvalues

write.xlsx(allRes_meso_vs_idio_Down.ByP_filtered, file = paste(results.dir_2,"byP_Down_resTable.xlsx", sep = "/"), row.names = F) 

#plot significant terms at 0.1 and 0.05 significance
n_term_sig <- sum(allRes_meso_vs_idio_Down.ByP_filtered$padj <= 0.1)

#Plot the sig Go terms 
idio_BP.down <- ggplot(allRes_meso_vs_idio_Down.ByP_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "Meso vs idio downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text=element_text(size=12),
                            axis.title=element_text(size=14,face="bold"), 
                            title = element_text(size=16,face="bold"))
  


png(paste(results.dir_2.1, "BP_Down_FDR_0.1.png", sep = "/"), width = 1100, height = 600)
idio_BP.down
dev.off()


n_term_sig <- sum(allRes_meso_vs_idio_Down.ByP_filtered$padj <= 0.05)

#Plot the sig Go terms 
idio_BP.down <- ggplot(allRes_meso_vs_idio_Down.ByP_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "Meso vs idio downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.05, "BP_Down_FDR_0.05.png", sep = "/"), width = 1100, height = 600)
idio_BP.down
dev.off()

allRes_meso_vs_idio.ByP_ordered <- allRes_meso_vs_idio_Down.ByP_filtered[which(as.numeric(allRes_meso_vs_idio_Down.ByP_filtered$classicFisher) < 0.01),]
allRes_meso_vs_idio.ByP_ordered <- allRes_meso_vs_idio.ByP_ordered[order(allRes_meso_vs_idio.ByP_ordered$Significant, decreasing = T),]


idio_BP_2.down <- ggplot(allRes_meso_vs_idio.ByP_ordered[c(1:20),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Mesophyll vs idioblast downregulated", 
       subtitle = "Biological process", 
       x = "", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_noFDR, "BP_Down_top20.png", sep = "/"), width = 4500, height = 3000, res = 300)
idio_BP_2.down
dev.off()


###Diff expressed ---- 
#(over representation of terms that are in both under and over represented)

#concatenate datasets from up and down annotated 
diff_idio.BP <- rbind(UP.annotated, Down.annotated)
genes_diff_expressed <- row.names(diff_idio.BP)

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_diff_expressed))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_meso_idio_diff_expressed.P <- new("topGOdata", ontology = "BP", 
                               description = "meso_vs_idio enrichment using kallisto quantification and dammit annotation",
                               allGenes = compared_genes, 
                               annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                               gene2GO = idio_GOByP,
                               nodeSize = 5) 

#Run fisher test
resultFisher_idio_diff_expressed.P <- runTest(GoData_meso_idio_diff_expressed.P, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_idio_diff_expressed.P

#4498 GO terms scored, 30 with p < 0.01
#(1682 out of 2562 significant genes, and 13406 feasible genes, or gene universe)

allRes_meso_vs_idio_diff_expressed.ByP <- GenTable(GoData_meso_idio_diff_expressed.P, 
                                         classicFisher = resultFisher_idio_diff_expressed.P,
                                         orderBy= "classicFisher", 
                                         numChar = 100, 
                                         topNodes = length(resultFisher_idio_diff_expressed.P@score))


#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the diff expressed genes 
#were removed 
allRes_meso_vs_idio_diff_expressed.ByP_filtered <- allRes_meso_vs_idio_diff_expressed.ByP[which(allRes_meso_vs_idio_diff_expressed.ByP$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_meso_vs_idio_diff_expressed.ByP_filtered$classicFisher, method = "fdr")

allRes_meso_vs_idio_diff_expressed.ByP_filtered$padj <- adj.pvalues

write.xlsx(allRes_meso_vs_idio_diff_expressed.ByP_filtered, file = paste(results.dir_2,"byP_diff_expressed_resTable.xlsx", sep = "/"), row.names = F) 

n_term_sig <- sum(allRes_meso_vs_idio_diff_expressed.ByP_filtered$padj <= 0.1)

#Plot the sig Go terms 
idio_BP.diff_expressed <- ggplot(allRes_meso_vs_idio_diff_expressed.ByP_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "Meso vs idio differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.1, "BP_diff_expressed_FDR_0.1.png", sep = "/"), width = 1100, height = 600)
idio_BP.diff_expressed
dev.off()

n_term_sig <- sum(allRes_meso_vs_idio_diff_expressed.ByP_filtered$padj <= 0.05)

#Plot the sig Go terms 
idio_BP.diff_expressed <- ggplot(allRes_meso_vs_idio_diff_expressed.ByP_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "Meso vs idio differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.05, "BP_diff_expressed_FDR_0.05.png", sep = "/"), width = 1100, height = 600)
idio_BP.diff_expressed
dev.off()



allRes_meso_vs_idio.ByP_ordered <- allRes_meso_vs_idio_diff_expressed.ByP_filtered[which(as.numeric(allRes_meso_vs_idio_diff_expressed.ByP_filtered$classicFisher) < 0.01),]
allRes_meso_vs_idio.ByP_ordered <- allRes_meso_vs_idio.ByP_ordered[order(allRes_meso_vs_idio.ByP_ordered$Significant, decreasing = T),]


idio_BP_2.diff_expressed <- ggplot(allRes_meso_vs_idio.ByP_ordered[c(1:20),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Mesophyll vs idioblast differential expressed", 
       subtitle = "Biological process", y = "Significant genes number", x = "") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_noFDR, "BP_diff_expressed_top20.png", sep = "/"), width = 4500, height = 3000, res = 300)
idio_BP_2.diff_expressed
dev.off()


##By Molecular function ----
#Get background annotation 
idio_GOMoF <- readMappings(file = "idio_GoF.txt")
bg_genes <- names(idio_GOMoF)


###Upregulated genes ----
#get gene IDs for the enrichement 
genes_UP_meso_vs_idio <- read.csv("DESeq2_tables/meso_vs_idio_sigUP.csv", header = TRUE)[,2]

#Get annotation for upregulated genes 
UP.annotated <- meso_vs_idio_sub_annot[genes_UP_meso_vs_idio,]

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_UP_meso_vs_idio))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_meso_idio_up.F <- new("topGOdata", ontology = "MF", 
                             description = "meso_vs_idio enrichment using kallisto quantification and dammit annotation",
                             allGenes = compared_genes, 
                             annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                             gene2GO = idio_GOMoF,
                             nodeSize = 5) 

#Run fisher test
resultFisher_idio_up.F <- runTest(GoData_meso_idio_up.F, 
                                  algorithm = "weight01", 
                                  statistic = "fisher")
resultFisher_idio_up.F

#4498 GO terms scored, 23 with p < 0.01
#(788 out of 1308 significant genes, and 13406 feasible genes, or gene universe)

allRes_meso_vs_idio_up.MoF <- GenTable(GoData_meso_idio_up.F, 
                                       classicFisher = resultFisher_idio_up.F,
                                       orderBy= "classicFisher", 
                                       numChar = 100, 
                                       topNodes = length(resultFisher_idio_up.F@score))

#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the upregulated genes 
#were removed 
allRes_meso_vs_idio_up.MoF_filtered <- allRes_meso_vs_idio_up.MoF[which(allRes_meso_vs_idio_up.MoF$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_meso_vs_idio_up.MoF_filtered$classicFisher, method = "fdr")

allRes_meso_vs_idio_up.MoF_filtered$padj <- adj.pvalues

write.xlsx(allRes_meso_vs_idio_up.MoF_filtered, file = paste(results.dir_2,"MoF_UP_resTable.xlsx", sep = "/"), row.names = F) 

n_term_sig <- sum(allRes_meso_vs_idio_up.MoF_filtered$padj <= 0.05)

#Plot the sig Go terms 
idio_MF.up <- ggplot(allRes_meso_vs_idio_up.MoF_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular Function", 
       subtitle = "Meso vs idio upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.05, "MF_UP_FDR_0.05.png", sep = "/"), width = 1000, height = 600)
idio_MF.up
dev.off()


n_term_sig <- sum(allRes_meso_vs_idio_up.MoF_filtered$padj <= 0.1)

#Plot the sig Go terms 
idio_MF.up <- ggplot(allRes_meso_vs_idio_up.MoF_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular Function", 
       subtitle = "Meso vs idio upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.1, "MF_UP_FDR_0.1.png", sep = "/"), width = 1000, height = 600)
idio_MF.up
dev.off()

allRes_meso_vs_idio.MoF_ordered <- allRes_meso_vs_idio_up.MoF_filtered[which(as.numeric(allRes_meso_vs_idio_up.MoF_filtered$classicFisher) < 0.01),]
allRes_meso_vs_idio.MoF_ordered <- allRes_meso_vs_idio.MoF_ordered[order(allRes_meso_vs_idio.MoF_ordered$Significant, decreasing = T),]


#Change the goterm name for plot reasons 
allRes_meso_vs_idio.MoF_ordered$Term[1] <- "oxidoreductase activity, acting on paired donors,\nwith incorporation or reduction of molecular oxygen"


idio_MF_2.up <- ggplot(allRes_meso_vs_idio.MoF_ordered[c(1:12),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Mesophyll vs idioblast upregulated", 
       subtitle = "Molecular function", y = "Significant genes number", x = " ") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))

png(paste(results.dir_noFDR, "MF_UP_top12.png", sep = "/"), width = 4500, height = 3000, res = 300)
idio_MF_2.up
dev.off()



###Downregulated genes ----

#get gene IDs for the enrichement 
genes_Down_meso_vs_idio <- read.csv("DESeq2_tables/meso_vs_idio_sigDOWN.csv", header = TRUE)[,2]

#Get annotation for upregulated genes 
Down.annotated <- meso_vs_idio_sub_annot[genes_Down_meso_vs_idio,]

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_Down_meso_vs_idio))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_meso_idio_Down.F <- new("topGOdata", ontology = "MF", 
                               description = "meso_vs_idio enrichment using kallisto quantification and dammit annotation",
                               allGenes = compared_genes, 
                               annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                               gene2GO = idio_GOMoF,
                               nodeSize = 5) 

#Run fisher test
resultFisher_idio_Down.F <- runTest(GoData_meso_idio_Down.F, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_idio_Down.F

#4498 GO terms scored, 28 with p < 0.01
#(788 out of 1308 significant genes, and 13406 feasible genes, or gene universe)

allRes_meso_vs_idio_Down.MoF <- GenTable(GoData_meso_idio_Down.F, 
                                         classicFisher = resultFisher_idio_Down.F,
                                         orderBy= "classicFisher", 
                                         numChar = 100, 
                                         topNodes = length(resultFisher_idio_Down.F@score))


#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the upregulated genes 
#were removed 
allRes_meso_vs_idio_Down.MoF_filtered <- allRes_meso_vs_idio_Down.MoF[which(allRes_meso_vs_idio_Down.MoF$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_meso_vs_idio_Down.MoF_filtered$classicFisher, method = "fdr")

allRes_meso_vs_idio_Down.MoF_filtered$padj <- adj.pvalues

write.xlsx(allRes_meso_vs_idio_Down.MoF_filtered, file = paste(results.dir_2,"MoF_Down_resTable.xlsx", sep = "/"), row.names = F) 

n_term_sig <- sum(allRes_meso_vs_idio_Down.MoF_filtered$padj <= 0.05)

#Plot the sig Go terms 
idio_MF.down <- ggplot(allRes_meso_vs_idio_Down.MoF_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular function", 
       subtitle = "Meso vs idio downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.05, "MF_Down_FDR_0.05.png", sep = "/"), width = 1300, height = 600)
idio_MF.down
dev.off()


n_term_sig <- sum(allRes_meso_vs_idio_Down.MoF_filtered$padj <= 0.1)

#Plot the sig Go terms 
idio_MF.down <- ggplot(allRes_meso_vs_idio_Down.MoF_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular function", 
       subtitle = "Meso vs idio downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.1, "MF_Down_FDR_0.1.png", sep = "/"), width = 1400, height = 600)
idio_MF.down
dev.off()

allRes_meso_vs_idio.MoF_ordered <- allRes_meso_vs_idio_Down.MoF_filtered[which(as.numeric(allRes_meso_vs_idio_Down.MoF_filtered$classicFisher) < 0.01),]
allRes_meso_vs_idio.MoF_ordered <- allRes_meso_vs_idio.MoF_ordered[order(allRes_meso_vs_idio.MoF_ordered$Significant, decreasing = T),]


allRes_meso_vs_idio.MoF_ordered$Term[3] <- "oxidoreductase activity, acting on paired donors,\nwith incorporation or reduction of molecular oxygen"
allRes_meso_vs_idio.MoF_ordered$Term[17] <- "electron transporter, transferring electrons within the \n cyclic electron transport pathway of photosynthesis activity"

idio_MF_2.down <- ggplot(allRes_meso_vs_idio.MoF_ordered[c(1:20),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Mesophyll vs idioblast downregulated", 
       subtitle = "Molecular function", y = "Significant genes number", x = "") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_noFDR, "MF_Down_top20.png", sep = "/"), width = 5000, height = 3000, res = 300)
idio_MF_2.down
dev.off()

###Diff expressed ---- 
#(over representation of terms that are in both under and over represented)

#concatenate datasets from up and down annotated 
diff_idio.MF <- rbind(UP.annotated, Down.annotated)
genes_diff_expressed <- row.names(diff_idio.MF)

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_diff_expressed))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_meso_idio_diff_expressed.F <- new("topGOdata", ontology = "MF", 
                                         description = "meso_vs_idio enrichment using kallisto quantification and dammit annotation",
                                         allGenes = compared_genes, 
                                         annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                                         gene2GO = idio_GOMoF,
                                         nodeSize = 5) 

#Run fisher test
resultFisher_idio_diff_expressed.F <- runTest(GoData_meso_idio_diff_expressed.F, 
                                              algorithm = "weight01", 
                                              statistic = "fisher")
resultFisher_idio_diff_expressed.F

#1308 GO terms scored, 27  with p < 0.01
#(1747 out of 2562 significant genes, and 13406 feasible genes, or gene universe)

allRes_meso_vs_idio_diff_expressed.MoF <- GenTable(GoData_meso_idio_diff_expressed.F, 
                                                   classicFisher = resultFisher_idio_diff_expressed.F,
                                                   orderBy= "classicFisher", 
                                                   numChar = 100, 
                                                   topNodes = length(resultFisher_idio_diff_expressed.F@score))


#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the diff expressed genes 
#were removed 
allRes_meso_vs_idio_diff_expressed.MoF_filtered <- allRes_meso_vs_idio_diff_expressed.MoF[which(allRes_meso_vs_idio_diff_expressed.MoF$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_meso_vs_idio_diff_expressed.MoF_filtered$classicFisher, method = "fdr")

allRes_meso_vs_idio_diff_expressed.MoF_filtered$padj <- adj.pvalues

write.xlsx(allRes_meso_vs_idio_diff_expressed.MoF_filtered, file = paste(results.dir_2,"MoF_diff_expressed_resTable.xlsx", sep = "/"), row.names = F) 

n_term_sig <- sum(allRes_meso_vs_idio_diff_expressed.MoF_filtered$padj <= 0.05)

#Plot the sig Go terms 
idio_MF.diff_expressed <- ggplot(allRes_meso_vs_idio_diff_expressed.MoF_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular function", 
       subtitle = "Meso vs idio differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.05, "MF_diff_expressed_FDR_0.05.png", sep = "/"), width = 1200, height = 600)
idio_MF.diff_expressed
dev.off()


n_term_sig <- sum(allRes_meso_vs_idio_diff_expressed.MoF_filtered$padj <= 0.1)

#Plot the sig Go terms 
idio_MF.diff_expressed <- ggplot(allRes_meso_vs_idio_diff_expressed.MoF_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular function", 
       subtitle = "Meso vs idio differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.1, "MF_diff_expressed_FDR_0.01.png", sep = "/"), width = 1200, height = 600)
idio_MF.diff_expressed
dev.off()

allRes_meso_vs_idio.MoF_ordered <- allRes_meso_vs_idio_diff_expressed.MoF_filtered[which(as.numeric(allRes_meso_vs_idio_diff_expressed.MoF_filtered$classicFisher) < 0.01),]
allRes_meso_vs_idio.MoF_ordered <- allRes_meso_vs_idio.MoF_ordered[order(allRes_meso_vs_idio.MoF_ordered$Significant, decreasing = T),]

#Fix some names 
allRes_meso_vs_idio.MoF_ordered$Term[2] <- "oxidoreductase activity, acting on paired donors\nwith incorporation or reduction of molecular oxygen"
allRes_meso_vs_idio.MoF_ordered$Term[7] <- "oxidoreductase activity, acting on the CH-OH group of donors,\n NAD or NADP as acceptor"
allRes_meso_vs_idio.MoF_ordered$Term[10] <- "oxidoreductase activity, acting on NAD(P)H\nquinone or similar compound as acceptor"

idio_MF_2.diff_expressed <- ggplot(allRes_meso_vs_idio.MoF_ordered[c(1:20),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Mesophyll vs idioblast differential expressed", 
       subtitle = "Molecular function", y = "Significant genes number", x = "") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))

png(paste(results.dir_noFDR, "MF_diff_expressed_top20.png", sep = "/"), width = 1000, height = 600)
idio_MF_2.diff_expressed
dev.off()



##By cellular_component ----
#Get background annotation 
idio_GOCeC <- readMappings(file = "idio_GoC.txt")
bg_genes <- names(idio_GOCeC)


###Upregulated genes ----
#get gene IDs for the enrichement 
genes_UP_meso_vs_idio <- read.csv("DESeq2_tables/meso_vs_idio_sigUP.csv", header = TRUE)[,2]

#Get annotation for upregulated genes 
UP.annotated <- meso_vs_idio_sub_annot[genes_UP_meso_vs_idio,]

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_UP_meso_vs_idio))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_meso_idio_up.C <- new("topGOdata", ontology = "CC", 
                             description = "meso_vs_idio enrichment using kallisto quantification and dammit annotation",
                             allGenes = compared_genes, 
                             annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                             gene2GO = idio_GOCeC,
                             nodeSize = 5) 

#Run fisher test
resultFisher_idio_up.C <- runTest(GoData_meso_idio_up.C, 
                                  algorithm = "weight01", 
                                  statistic = "fisher")
resultFisher_idio_up.C

#825 GO terms scored, 4 with p < 0.01
#(788 out of 831 significant genes, and 13406 feasible genes, or gene universe)

allRes_meso_vs_idio_up.CeC <- GenTable(GoData_meso_idio_up.C, 
                                       classicFisher = resultFisher_idio_up.C,
                                       orderBy= "classicFisher", 
                                       numChar = 100, 
                                       topNodes = length(resultFisher_idio_up.C@score))

#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the upregulated genes 
#were removed 
allRes_meso_vs_idio_up.CeC_filtered <- allRes_meso_vs_idio_up.CeC[which(allRes_meso_vs_idio_up.CeC$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_meso_vs_idio_up.CeC_filtered$classicFisher, method = "fdr")

allRes_meso_vs_idio_up.CeC_filtered$padj <- adj.pvalues

write.xlsx(allRes_meso_vs_idio_up.CeC_filtered, file = paste(results.dir_2,"CeC_UP_resTable.xlsx", sep = "/"), row.names = F) 

n_term_sig <- sum(allRes_meso_vs_idio_up.CeC_filtered$padj <= 0.05)

#Plot the sig Go terms 
idio_CC.up <- ggplot(allRes_meso_vs_idio_up.CeC_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular compartement", 
       subtitle = "Meso vs idio upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.05, "CC_UP_FDR_0.05.png", sep = "/"), width = 1200, height = 600)
idio_CC.up
dev.off()

allRes_meso_vs_idio.CeC_ordered <- allRes_meso_vs_idio_up.CeC_filtered[which(as.numeric(allRes_meso_vs_idio_up.CeC_filtered$classicFisher) < 0.01),]
allRes_meso_vs_idio.CeC_ordered <- allRes_meso_vs_idio.CeC_ordered[order(allRes_meso_vs_idio.CeC_ordered$Significant, decreasing = T),]


n_term_sig <- sum(allRes_meso_vs_idio_up.CeC_filtered$padj <= 0.1)

#Plot the sig Go terms 
idio_CC.up <- ggplot(allRes_meso_vs_idio_up.CeC_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular compartement", 
       subtitle = "Meso vs idio upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.1, "CC_UP_FDR_0.1.png", sep = "/"), width = 1200, height = 600)
idio_CC.up
dev.off()

allRes_meso_vs_idio.CeC_ordered <- allRes_meso_vs_idio_up.CeC_filtered[which(as.numeric(allRes_meso_vs_idio_up.CeC_filtered$classicFisher) < 0.01),]
allRes_meso_vs_idio.CeC_ordered <- allRes_meso_vs_idio.CeC_ordered[order(allRes_meso_vs_idio.CeC_ordered$Significant, decreasing = T),]



idio_CC_2.up <- ggplot(allRes_meso_vs_idio.CeC_ordered[c(1:3),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Mesophyll vs idioblast upregulated", 
       subtitle = "Cellular component", y = "Significant genes number", x = "") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_noFDR, "CC_UP_top3.png", sep = "/"), width = 4500, height = 3000, res = 300)
idio_CC_2.up
dev.off()


 
###Downregulated genes ----

#get gene IDs for the enrichement 
genes_Down_meso_vs_idio <- read.csv("DESeq2_tables/meso_vs_idio_sigDOWN.csv", header = TRUE)[,2]

#Get annotation for upregulated genes 
Down.annotated <- meso_vs_idio_sub_annot[genes_Down_meso_vs_idio,]

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_Down_meso_vs_idio))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_meso_idio_Down.C <- new("topGOdata", ontology = "CC", 
                               description = "meso_vs_idio enrichment using kallisto quantification and dammit annotation",
                               allGenes = compared_genes, 
                               annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                               gene2GO = idio_GOCeC,
                               nodeSize = 5) 

#Run fisher test
resultFisher_idio_Down.C <- runTest(GoData_meso_idio_Down.C, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_idio_Down.C

#4498 GO terms scored, 28 with p < 0.01
#(788 out of 1308 significant genes, and 13406 feasible genes, or gene universe)

allRes_meso_vs_idio_Down.CeC <- GenTable(GoData_meso_idio_Down.C, 
                                         classicFisher = resultFisher_idio_Down.C,
                                         orderBy= "classicFisher", 
                                         numChar = 100, 
                                         topNodes = length(resultFisher_idio_Down.C@score))

#one of the values has "<" 
allRes_meso_vs_idio_Down.CeC$classicFisher <- str_replace(allRes_meso_vs_idio_Down.CeC$classicFisher, "< ", "") 

#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the upregulated genes 
#were removed 
allRes_meso_vs_idio_Down.CeC_filtered <- allRes_meso_vs_idio_Down.CeC[which(allRes_meso_vs_idio_Down.CeC$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_meso_vs_idio_Down.CeC_filtered$classicFisher, method = "fdr")

allRes_meso_vs_idio_Down.CeC_filtered$padj <- adj.pvalues

write.xlsx(allRes_meso_vs_idio_Down.CeC_filtered, file = paste(results.dir_2,"CeC_Down_resTable.xlsx", sep = "/"), row.names = F) 

n_term_sig <- sum(allRes_meso_vs_idio_Down.CeC_filtered$padj <= 0.05)

#Plot the sig Go terms 
idio_CC.down <- ggplot(allRes_meso_vs_idio_Down.CeC_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular compartement", 
       subtitle = "Meso vs idio downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.05, "CC_Down_FDR_0.05.png", sep = "/"), width = 1300, height = 600)
idio_CC.down
dev.off()

n_term_sig <- sum(allRes_meso_vs_idio_Down.CeC_filtered$padj <= 0.1)

#Plot the sig Go terms 
idio_CC.down <- ggplot(allRes_meso_vs_idio_Down.CeC_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular compartement", 
       subtitle = "Meso vs idio downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.1, "CC_Down_FDR_0.1.png", sep = "/"), width = 1300, height = 600)
idio_CC.down
dev.off()


allRes_meso_vs_idio.CeC_ordered <- allRes_meso_vs_idio_Down.CeC_filtered[which(as.numeric(allRes_meso_vs_idio_Down.CeC_filtered$classicFisher) < 0.01),]
allRes_meso_vs_idio.CeC_ordered <- allRes_meso_vs_idio.CeC_ordered[order(allRes_meso_vs_idio.CeC_ordered$Significant, decreasing = T),]


idio_CC_2.down <- ggplot(allRes_meso_vs_idio.CeC_ordered[c(1:20),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Mesophyll vs idioblast downregulated", 
       subtitle = "Cellular compartement", y = "Significant genes number", x = "") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_noFDR, "CC_Down_top20.png", sep = "/"), width = 4500, height = 3000, res = 300)
idio_CC_2.down
dev.off()

###Dif expressed ---- 
#(over representation of terms that are in both under and over represented)

#concatenate datasets from up and down annotated 
diff_idio.CC <- rbind(UP.annotated, Down.annotated)
genes_diff_expressed <- row.names(diff_idio.CC)

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_diff_expressed))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_meso_idio_diff_expressed.C <- new("topGOdata", ontology = "CC", 
                                         description = "meso_vs_idio enrichment using kallisto quantification and dammit annotation",
                                         allGenes = compared_genes, 
                                         annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                                         gene2GO = idio_GOCeC,
                                         nodeSize = 5) 

#Run fisher test
resultFisher_idio_diff_expressed.C <- runTest(GoData_meso_idio_diff_expressed.C, 
                                              algorithm = "weight01", 
                                              statistic = "fisher")
resultFisher_idio_diff_expressed.C

#825 GO terms scored, 23 with p < 0.01
#(1815 out of 2562 significant genes, and 13406 feasible genes, or gene universe)

allRes_meso_vs_idio_diff_expressed.CeC <- GenTable(GoData_meso_idio_diff_expressed.C, 
                                                   classicFisher = resultFisher_idio_diff_expressed.C,
                                                   orderBy= "classicFisher", 
                                                   numChar = 100, 
                                                   topNodes = length(resultFisher_idio_diff_expressed.C@score))


#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the diff expressed genes 
#were removed 
allRes_meso_vs_idio_diff_expressed.CeC_filtered <- allRes_meso_vs_idio_diff_expressed.CeC[which(allRes_meso_vs_idio_diff_expressed.CeC$Significant > 0),]

allRes_meso_vs_idio_diff_expressed.CeC_filtered$classicFisher[1] <- 1^-30

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_meso_vs_idio_diff_expressed.CeC_filtered$classicFisher, method = "fdr")

allRes_meso_vs_idio_diff_expressed.CeC_filtered$padj <- adj.pvalues

write.xlsx(allRes_meso_vs_idio_diff_expressed.CeC_filtered, file = paste(results.dir_2,"CeC_diff_expressed_resTable.xlsx", sep = "/"), row.names = F) 

n_term_sig <- sum(allRes_meso_vs_idio_diff_expressed.CeC_filtered$padj <= 0.05)

#Plot the sig Go terms 
idio_CC.diff_expressed <- ggplot(allRes_meso_vs_idio_diff_expressed.CeC_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular function", 
       subtitle = "Meso vs idio differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.05, "CC_diff_expressed_FRD_0.05.png", sep = "/"), width = 1200, height = 600)
idio_CC.diff_expressed
dev.off()


n_term_sig <- sum(allRes_meso_vs_idio_diff_expressed.CeC_filtered$padj <= 0.1)

#Plot the sig Go terms 
idio_CC.diff_expressed <- ggplot(allRes_meso_vs_idio_diff_expressed.CeC_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular function", 
       subtitle = "Meso vs idio differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.1, "CC_diff_expressed_FDR_0.1.png", sep = "/"), width = 1200, height = 600)
idio_CC.diff_expressed
dev.off()



allRes_meso_vs_idio.CeC_ordered <- allRes_meso_vs_idio_diff_expressed.CeC_filtered[which(as.numeric(allRes_meso_vs_idio_diff_expressed.CeC_filtered$classicFisher) < 0.01),]
allRes_meso_vs_idio.CeC_ordered <- allRes_meso_vs_idio.CeC_ordered[order(allRes_meso_vs_idio.CeC_ordered$Significant, decreasing = T),]


idio_CC_2.diff_expressed <- ggplot(allRes_meso_vs_idio.CeC_ordered[c(1:20),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Mesophyll vs idioblast differential expressed", 
       subtitle = "Cellular component", y = "Significant genes number", x = "") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=16),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))

png(paste(results.dir_noFDR, "CC_diff_expressed_top20.png", sep = "/"),  width = 4500, height = 3000, res = 300)
idio_CC_2.diff_expressed
dev.off()


############################################
#### F1in vs F1out enrichement analysis ####
############################################

#make results directory
results.dir_2 <- paste(results.dir, "F1in_vs_f1out", sep = "/")
if (!dir.exists(results.dir_2)){
  dir.create(results.dir_2)
}


results.dir_2.1 <- paste(results.dir_2, "FDR_0.1", sep = "/")
if (!dir.exists(results.dir_2.1)){
  dir.create(results.dir_2.1)
}

results.dir_2.05 <- paste(results.dir_2, "FDR_0.05", sep = "/")
if (!dir.exists(results.dir_2.05)){
  dir.create(results.dir_2.05)
}


results.dir_noFDR <- paste(results.dir_2, "Before_muliple_testing", sep = "/")
if (!dir.exists(results.dir_noFDR)){
  dir.create(results.dir_noFDR)
}


f1_bg_genes <- read.csv(file = "DESeq2_tables/F1in_vs_F1out_Res.csv")$X
f1_sub_annot <- annotation_file[f1_bg_genes,]

###Make backgroud annotation ---- 
#For each compariton, the background annotation is different as we only consider the gene with
#gene counts > 0 for the compariton. For this filtering, the genes in the res table from DesEq will be used, 
#as genes with rowcounts = 0 were filtered out before 


##Subset annotation files for backgroud annotation 
##Go(P) <- biological process
f1_Go.P <- data.frame("gene" = row.names(f1_sub_annot), "GO.BP" = f1_sub_annot[,9])
f1_Go.P$'GO.BP' <- str_replace_all(f1_Go.P$'GO.BP', ";", ",")
write.table(f1_Go.P, file = "f1_GoP.txt", row.names = FALSE, col.names = F,quote = F, sep = "\t")

##Go(F) <- molecular function
f1_Go.F <- data.frame("gene" = row.names(f1_sub_annot), "GO.MF" = f1_sub_annot[,11])
f1_Go.F$'GO.MF' <- str_replace_all(f1_Go.F$'GO.MF', ";", ",")
write.table(f1_Go.F, file = "f1_GoF.txt", row.names = FALSE, col.names = F,quote = F, sep = "\t")

##Go(C) <- celular compartement
f1_Go.C <- data.frame("gene" = row.names(f1_sub_annot), "GO.CC" = f1_sub_annot[,13])
f1_Go.C$'GO.CC' <- str_replace_all(f1_Go.C$'GO.CC', ";", ",")
write.table(f1_Go.C, file = "f1_GoC.txt", row.names = FALSE, col.names = F,quote = F, sep = "\t")


##By Biological process ----
#Get background annotation 
f1_GOByP <- readMappings(file = "f1_GoP.txt")
bg_genes <- names(f1_GOByP)


###Upregulated genes ----
#get gene IDs for the enrichement 
genes_UP_f1 <- read.csv("DESeq2_tables/F1in_vs_F1out_sigUP.csv", header = TRUE)[,2]

#Get annotation for upregulated genes 
UP.annotated <- f1_sub_annot[genes_UP_f1,]

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_UP_f1))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_f1_up.P <- new("topGOdata", ontology = "BP", 
                             description = "f1in vs f1out enrichment using kallisto quantification and dammit annotation",
                             allGenes = compared_genes, 
                             annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                             gene2GO = f1_GOByP,
                             nodeSize = 5) 

#Run fisher test
resultFisher_f1_up.P <- runTest(GoData_f1_up.P, 
                                  algorithm = "weight01", 
                                  statistic = "fisher")
resultFisher_f1_up.P

#4393 GO terms scored, 21 with p < 0.01
#(417 out of 705 significant genes, and 13406 feasible genes, or gene universe)

allRes_f1_up.ByP <- GenTable(GoData_f1_up.P, 
                                       classicFisher = resultFisher_f1_up.P,
                                       orderBy= "classicFisher", 
                                       numChar = 100, 
                                       topNodes = length(resultFisher_f1_up.P@score))

#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the upregulated genes 
#were removed 
allRes_f1_up.ByP_filtered <- allRes_f1_up.ByP[which(allRes_f1_up.ByP$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_f1_up.ByP_filtered$classicFisher, method = "fdr")

allRes_f1_up.ByP_filtered$padj <- adj.pvalues

write.xlsx(allRes_f1_up.ByP_filtered, file = paste(results.dir_2,"byP_UP_resTable.xlsx", sep = "/"), row.names = F) 

n_term_sig <- sum(allRes_f1_up.ByP_filtered$padj <= 0.05)

#Plot the sig Go terms 
f1_BP.up <- ggplot(allRes_f1_up.ByP_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "f1in vs f1out upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.05, "BP_UP_FDR_0.05.png", sep = "/"), width = 1000, height = 600)
f1_BP.up
dev.off()


n_term_sig <- sum(allRes_f1_up.ByP_filtered$padj <= 0.1)

#Plot the sig Go terms 
f1_BP.up <- ggplot(allRes_f1_up.ByP_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "f1in vs f1out upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.1, "BP_UP_FDR_0.1.png", sep = "/"), width = 1000, height = 600)
f1_BP.up
dev.off()

allRes_f1.ByP_ordered <- allRes_f1_up.ByP_filtered[which(as.numeric(allRes_f1_up.ByP_filtered$classicFisher) < 0.01),]
allRes_f1.ByP_ordered <- allRes_f1.ByP_ordered[order(allRes_f1.ByP_ordered$Significant, decreasing = T),]


f1_BP_2.up <- ggplot(allRes_f1.ByP_ordered[c(1:20),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "F1in vs F1out upregulated", 
       subtitle = "Biological process", y = "Significant genes number", x = "") + 
  theme_classic() +   
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_noFDR, "BP_UP_top20.png", sep = "/"), width = 4500, height = 3000, res = 300)
f1_BP_2.up
dev.off()



###Downregulated genes ----

#get gene IDs for the enrichement 
genes_Down_f1 <- read.csv("DESeq2_tables/F1in_vs_F1out_sigDOWN.csv", header = TRUE)[,2]

#Get annotation for upregulated genes 
Down.annotated <- f1_sub_annot[genes_Down_f1,]

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_Down_f1))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_f1_Down.P <- new("topGOdata", ontology = "BP", 
                               description = "f1 enrichment using kallisto quantification and dammit annotation",
                               allGenes = compared_genes, 
                               annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                               gene2GO = f1_GOByP,
                               nodeSize = 5) 

#Run fisher test
resultFisher_f1_Down.P <- runTest(GoData_f1_Down.P, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_f1_Down.P

#4393 GO terms scored, 19 with p < 0.01
#(304 out of 509 significant genes, and 13406 feasible genes, or gene universe)

allRes_f1_Down.ByP <- GenTable(GoData_f1_Down.P, 
                                         classicFisher = resultFisher_f1_Down.P,
                                         orderBy= "classicFisher", 
                                         numChar = 100, 
                                         topNodes = length(resultFisher_f1_Down.P@score))


#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the upregulated genes 
#were removed 
allRes_f1_Down.ByP_filtered <- allRes_f1_Down.ByP[which(allRes_f1_Down.ByP$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_f1_Down.ByP_filtered$classicFisher, method = "fdr")

allRes_f1_Down.ByP_filtered$padj <- adj.pvalues

write.xlsx(allRes_f1_Down.ByP_filtered, file = paste(results.dir_2,"byP_Down_resTable.xlsx", sep = "/"), row.names = F) 

n_term_sig <- sum(allRes_f1_Down.ByP_filtered$padj <= 0.05)

#Plot the sig Go terms 
f1_BP.down <- ggplot(allRes_f1_Down.ByP_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "f1in vs f1out downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))




png(paste(results.dir_2.05, "BP_Down_FDR_0.05.png", sep = "/"), width = 1100, height = 600)
f1_BP.down
dev.off()


n_term_sig <- sum(allRes_f1_Down.ByP_filtered$padj <= 0.1)

#Plot the sig Go terms 
f1_BP.down <- ggplot(allRes_f1_Down.ByP_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "f1in vs f1out downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))




png(paste(results.dir_2.1, "BP_Down_FDR_0.1.png", sep = "/"), width = 1100, height = 600)
f1_BP.down
dev.off()


allRes_f1.ByP_ordered <- allRes_f1_Down.ByP_filtered[which(as.numeric(allRes_f1_Down.ByP_filtered$classicFisher) < 0.01),]
allRes_f1.ByP_ordered <- allRes_f1.ByP_ordered[order(allRes_f1.ByP_ordered$Significant, decreasing = T),]


f1_BP_2.down <- ggplot(allRes_f1.ByP_ordered[c(1:20),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "F1in vs F1out downregulated", 
       subtitle = "Biological process", y = "Significant genes number", x = "") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_noFDR, "BP_Down_top20.png", sep = "/"), width = 4500, height = 3000, res = 300)
f1_BP_2.down
dev.off()

###Dif expressed ---- 
#(over representation of terms that are in both under and over represented)

#concatenate datasets from up and down annotated 
diff_f1.BP <- rbind(UP.annotated, Down.annotated)
genes_diff_expressed <- row.names(diff_f1.BP)

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_diff_expressed))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_f1_diff_expressed.P <- new("topGOdata", ontology = "BP", 
                                         description = "f1in_vs_f1out enrichment using kallisto quantification and dammit annotation",
                                         allGenes = compared_genes, 
                                         annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                                         gene2GO = f1_GOByP,
                                         nodeSize = 5) 

#Run fisher test
resultFisher_f1_diff_expressed.P <- runTest(GoData_f1_diff_expressed.P, 
                                              algorithm = "weight01", 
                                              statistic = "fisher")
resultFisher_f1_diff_expressed.P

#4498 GO terms scored, 28 with p < 0.01
#(757 out of 1214 significant genes, and 13406 feasible genes, or gene universe)

allRes_f1in_diff_expressed.ByP <- GenTable(GoData_f1_diff_expressed.P, 
                                                   classicFisher = resultFisher_f1_diff_expressed.P,
                                                   orderBy= "classicFisher", 
                                                   numChar = 100, 
                                                   topNodes = length(resultFisher_f1_diff_expressed.P@score))


#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the diff expressed genes 
#were removed 
allRes_f1in_diff_expressed.ByP_filtered <- allRes_f1in_diff_expressed.ByP[which(allRes_f1in_diff_expressed.ByP$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_f1in_diff_expressed.ByP_filtered$classicFisher, method = "fdr")

allRes_f1in_diff_expressed.ByP_filtered$padj <- adj.pvalues

write.xlsx(allRes_f1in_diff_expressed.ByP_filtered, file = paste(results.dir_2,"byP_diff_expressed_resTable.xlsx", sep = "/"), row.names = F) 

n_term_sig <- sum(allRes_f1in_diff_expressed.ByP_filtered$padj <= 0.05)

#Plot the sig Go terms 
f1_BP.diff_expressed <- ggplot(allRes_f1in_diff_expressed.ByP_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "f1in vs f1out differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.05, "BP_diff_expressed_FDR_0.05.png", sep = "/"), width = 1100, height = 600)
f1_BP.diff_expressed
dev.off()


n_term_sig <- sum(allRes_f1in_diff_expressed.ByP_filtered$padj <= 0.1)

#Plot the sig Go terms 
f1_BP.diff_expressed <- ggplot(allRes_f1in_diff_expressed.ByP_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "f1in vs f1out differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.1, "BP_diff_expressed_FDR_0.1.png", sep = "/"), width = 1100, height = 600)
f1_BP.diff_expressed
dev.off()


allRes_f1in.ByP_ordered <- allRes_f1in_diff_expressed.ByP_filtered[which(as.numeric(allRes_f1in_diff_expressed.ByP_filtered$classicFisher) < 0.01),]
allRes_f1in.ByP_ordered <- allRes_f1in.ByP_ordered[order(allRes_f1in.ByP_ordered$Significant, decreasing = T),]


f1_BP_2.diff_expressed <- ggplot(allRes_f1in.ByP_ordered[c(1:20),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "F1in vs F1out differential expressed", 
       subtitle = "Biological process", y = "Significant genes number", x = "") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_noFDR, "BP_diff_expressed_top20.png", sep = "/"), width = 1000, height = 600)
f1_BP_2.diff_expressed
dev.off()



##By Molecular function ----
#Get background annotation 
f1_GOMoF <- readMappings(file = "f1_GoF.txt")
bg_genes <- names(f1_GOMoF)


###Upregulated genes ----
#get gene IDs for the enrichement 
genes_UP_f1 <- read.csv("DESeq2_tables/F1in_vs_F1out_sigUP.csv", header = TRUE)[,2]

#Get annotation for upregulated genes 
UP.annotated <- f1_sub_annot[genes_UP_f1,]

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_UP_f1))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_f1_up.F <- new("topGOdata", ontology = "MF", 
                             description = "f1 enrichment using kallisto quantification and dammit annotation",
                             allGenes = compared_genes, 
                             annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                             gene2GO = f1_GOMoF,
                             nodeSize = 5) 

#Run fisher test
resultFisher_f1_up.F <- runTest(GoData_f1_up.F, 
                                  algorithm = "weight01", 
                                  statistic = "fisher")
resultFisher_f1_up.F

#1317 GO terms scored, 16 with p < 0.01
#(454 out of 705 significant genes, and 13470 feasible genes, or gene universe)

allRes_f1_up.MoF <- GenTable(GoData_f1_up.F, 
                                       classicFisher = resultFisher_f1_up.F,
                                       orderBy= "classicFisher", 
                                       numChar = 100, 
                                       topNodes = length(resultFisher_f1_up.F@score))

#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the upregulated genes 
#were removed 
allRes_f1_up.MoF_filtered <- allRes_f1_up.MoF[which(allRes_f1_up.MoF$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_f1_up.MoF_filtered$classicFisher, method = "fdr")

allRes_f1_up.MoF_filtered$padj <- adj.pvalues

write.xlsx(allRes_f1_up.MoF_filtered, file = paste(results.dir_2,"moF_UP_resTable.xlsx", sep = "/"), row.names = F) 

n_term_sig <- sum(allRes_f1_up.MoF_filtered$padj <= 0.05)

#Plot the sig Go terms 
f1_MF.up <- ggplot(allRes_f1_up.MoF_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular Function", 
       subtitle = "f1 upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.05, "MF_UP_sig_FDR_0.05.png", sep = "/"), width = 1000, height = 600)
f1_MF.up
dev.off()



n_term_sig <- sum(allRes_f1_up.MoF_filtered$padj <= 0.1)

#Plot the sig Go terms 
f1_MF.up <- ggplot(allRes_f1_up.MoF_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular Function", 
       subtitle = "f1 upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.1, "MF_UP_sig_FDR_0.1.png", sep = "/"), width = 1000, height = 600)
f1_MF.up
dev.off()



allRes_f1.MoF_ordered <- allRes_f1_up.MoF_filtered[which(as.numeric(allRes_f1_up.MoF_filtered$classicFisher) < 0.01),]
allRes_f1.MoF_ordered <- allRes_f1.MoF_ordered[order(allRes_f1.MoF_ordered$Significant, decreasing = T),]


#fix names 
allRes_f1.MoF_ordered$Term[11] <- "oxidoreductase activity, acting on single donors with incorporation\nof molecular oxygen, incorporation of two atoms of oxygen"

f1_MF_2.up <- ggplot(allRes_f1.MoF_ordered[c(1:18),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "F1in vs F1out upregulated", 
       subtitle = "Molecular function", y = "Significant genes number", x = "") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_noFDR, "MF_UP_top18.png", sep = "/"), width = 4500, height = 3000, res = 300)
f1_MF_2.up
dev.off()



###Downregulated genes ----

#get gene IDs for the enrichement 
genes_Down_f1 <- read.csv("DESeq2_tables/F1in_vs_F1out_sigDOWN.csv", header = TRUE)[,2]

#Get annotation for upregulated genes 
Down.annotated <- f1_sub_annot[genes_Down_f1,]

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_Down_f1))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_f1_Down.F <- new("topGOdata", ontology = "MF", 
                               description = "f1 enrichment using kallisto quantification and dammit annotation",
                               allGenes = compared_genes, 
                               annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                               gene2GO = f1_GOMoF,
                               nodeSize = 5) 

#Run fisher test
resultFisher_f1_Down.F <- runTest(GoData_f1_Down.F, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_f1_Down.F

#1317 GO terms scored, 15 with p < 0.01
#(336  out of 509 significant genes, and 13406 feasible genes, or gene universe)

allRes_f1_Down.MoF <- GenTable(GoData_f1_Down.F, 
                                         classicFisher = resultFisher_f1_Down.F,
                                         orderBy= "classicFisher", 
                                         numChar = 100, 
                                         topNodes = length(resultFisher_f1_Down.F@score))


#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the upregulated genes 
#were removed 
allRes_f1_Down.MoF_filtered <- allRes_f1_Down.MoF[which(allRes_f1_Down.MoF$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_f1_Down.MoF_filtered$classicFisher, method = "fdr")

allRes_f1_Down.MoF_filtered$padj <- adj.pvalues

write.xlsx(allRes_f1_Down.MoF_filtered, file = paste(results.dir_2,"MoF_Down_resTable.xlsx", sep = "/"), row.names = F) 

n_term_sig <- sum(allRes_f1_Down.MoF_filtered$padj <= 0.05)

#Plot the sig Go terms 
f1_MF.down <- ggplot(allRes_f1_Down.MoF_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular function", 
       subtitle = "f1 downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.05, "MF_Down_sig_FDR_0.05.png", sep = "/"), width = 1100, height = 600)
f1_MF.down
dev.off()


n_term_sig <- sum(allRes_f1_Down.MoF_filtered$padj <= 0.1)

#Plot the sig Go terms 
f1_MF.down <- ggplot(allRes_f1_Down.MoF_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular function", 
       subtitle = "f1 downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.1, "MF_Down_sig_FDR_0.1.png", sep = "/"), width = 1100, height = 600)
f1_MF.down
dev.off()


allRes_f1.MoF_ordered <- allRes_f1_Down.MoF_filtered[which(as.numeric(allRes_f1_Down.MoF_filtered$classicFisher) < 0.01),]
allRes_f1.MoF_ordered <- allRes_f1.MoF_ordered[order(allRes_f1.MoF_ordered$Significant, decreasing = T),]

allRes_f1.MoF_ordered$Term[10] <- "oxidoreductase activity, acting on the CH-NH2\ngroup of donors, NAD or NADP as acceptor"

f1_MF_2.down <- ggplot(allRes_f1.MoF_ordered[c(1:11),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular function", 
       subtitle = "f1in vs f1out downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))

png(paste(results.dir_noFDR, "MF_Down_top11.png", sep = "/"), width = 1000, height = 600)
f1_MF_2.down
dev.off()

###Dif expressed ---- 
#(over representation of terms that are in both under and over represented)

#concatenate datasets from up and down annotated 
diff_f1.MF <- rbind(UP.annotated, Down.annotated)
genes_diff_expressed <- row.names(diff_f1.MF)

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_diff_expressed))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_f1_diff_expressed.M <- new("topGOdata", ontology = "MF", 
                                  description = "f1in_vs_f1out enrichment using kallisto quantification and dammit annotation",
                                  allGenes = compared_genes, 
                                  annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                                  gene2GO = f1_GOMoF,
                                  nodeSize = 5) 

#Run fisher test
resultFisher_f1_diff_expressed.M <- runTest(GoData_f1_diff_expressed.M, 
                                            algorithm = "weight01", 
                                            statistic = "fisher")
resultFisher_f1_diff_expressed.M

#1317 GO terms scored, 27 with p < 0.01
#(790 out of 1214 significant genes, and 13406 feasible genes, or gene universe)

allRes_f1in_diff_expressed.MoF <- GenTable(GoData_f1_diff_expressed.M, 
                                           classicFisher = resultFisher_f1_diff_expressed.M,
                                           orderBy= "classicFisher", 
                                           numChar = 100, 
                                           topNodes = length(resultFisher_f1_diff_expressed.M@score))


#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the diff expressed genes 
#were removed 
allRes_f1in_diff_expressed.MoF_filtered <- allRes_f1in_diff_expressed.MoF[which(allRes_f1in_diff_expressed.MoF$Significant > 0),]

#correct for multiple testing with FDR
adj.Mvalues <- p.adjust(allRes_f1in_diff_expressed.MoF_filtered$classicFisher, method = "fdr")

allRes_f1in_diff_expressed.MoF_filtered$padj <- adj.Mvalues

write.xlsx(allRes_f1in_diff_expressed.MoF_filtered, file = paste(results.dir_2,"MoF_diff_expressed_resTable.xlsx", sep = "/"), row.names = F) 

n_term_sig <- sum(allRes_f1in_diff_expressed.MoF_filtered$padj <= 0.05)

#Plot the sig Go terms 
f1_MF.diff_expressed <- ggplot(allRes_f1in_diff_expressed.MoF_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular function", 
       subtitle = "f1in vs f1out differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.05, "MF_diff_expressed_FDR_0.05.png", sep = "/"), width = 1100, height = 600)
f1_MF.diff_expressed
dev.off()


n_term_sig <- sum(allRes_f1in_diff_expressed.MoF_filtered$padj <= 0.1)

#Plot the sig Go terms 
f1_MF.diff_expressed <- ggplot(allRes_f1in_diff_expressed.MoF_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular function", 
       subtitle = "f1in vs f1out differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.1, "MF_diff_expressed_FDR_0.1.png", sep = "/"), width = 1100, height = 600)
f1_MF.diff_expressed
dev.off()

allRes_f1in.MoF_ordered <- allRes_f1in_diff_expressed.MoF_filtered[which(as.numeric(allRes_f1in_diff_expressed.MoF_filtered$classicFisher) < 0.01),]
allRes_f1in.MoF_ordered <- allRes_f1in.MoF_ordered[order(allRes_f1in.MoF_ordered$Significant, decreasing = T),]


f1_MF_2.diff_expressed <- ggplot(allRes_f1in.MoF_ordered[c(1:20),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular function", 
       subtitle = "f1in vs f1out differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))




png(paste(results.dir_noFDR, "MF_diff_expressed_top20.png", sep = "/"), width = 1000, height = 600)
f1_MF_2.diff_expressed
dev.off()


##By cellular_component ----
#Get background annotation 
f1_GOCeC <- readMappings(file = "f1_GoC.txt")
bg_genes <- names(f1_GOCeC)


###Upregulated genes ----
#get gene IDs for the enrichement 
genes_UP_f1 <- read.csv("DESeq2_tables/F1in_vs_F1out_sigUP.csv", header = TRUE)[,2]

#Get annotation for upregulated genes 
UP.annotated <- f1_sub_annot[genes_UP_f1,]

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_UP_f1))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_f1_up.C <- new("topGOdata", ontology = "CC", 
                             description = "f1 enrichment using kallisto quantification and dammit annotation",
                             allGenes = compared_genes, 
                             annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                             gene2GO = f1_GOCeC,
                             nodeSize = 5) 

#Run fisher test
resultFisher_f1_up.C <- runTest(GoData_f1_up.C, 
                                  algorithm = "weight01", 
                                  statistic = "fisher")
resultFisher_f1_up.C

#819 GO terms scored, 7 with p < 0.01
#(788 out of 831 significant genes, and 13406 feasible genes, or gene universe)

allRes_f1_up.CeC <- GenTable(GoData_f1_up.C, 
                                       classicFisher = resultFisher_f1_up.C,
                                       orderBy= "classicFisher", 
                                       numChar = 100, 
                                       topNodes = length(resultFisher_f1_up.C@score))

#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the upregulated genes 
#were removed 
allRes_f1_up.CeC_filtered <- allRes_f1_up.CeC[which(allRes_f1_up.CeC$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_f1_up.CeC_filtered$classicFisher, method = "fdr")

allRes_f1_up.CeC_filtered$padj <- adj.pvalues

write.xlsx(allRes_f1_up.CeC_filtered, file = paste(results.dir_2,"CeC_UP_resTable.xlsx", sep = "/"), row.names = F) 

n_term_sig <- sum(allRes_f1_up.CeC_filtered$padj <= 0.05)

#Plot the sig Go terms 
f1_CC.up <- ggplot(allRes_f1_up.CeC_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular compartement", 
       subtitle = "f1 upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.05, "CC_UP_sig_FDR_0.05.png", sep = "/"), width = 1000, height = 600)
f1_CC.up
dev.off()


n_term_sig <- sum(allRes_f1_up.CeC_filtered$padj <= 0.1)

#Plot the sig Go terms 
f1_CC.up <- ggplot(allRes_f1_up.CeC_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular compartement", 
       subtitle = "f1 upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))




png(paste(results.dir_2.1, "CC_UP_sig_FDR_0.1.png", sep = "/"), width = 1000, height = 600)
f1_CC.up
dev.off()

allRes_f1.CeC_ordered <- allRes_f1_up.CeC_filtered[which(as.numeric(allRes_f1_up.CeC_filtered$classicFisher) < 0.01),]
allRes_f1.CeC_ordered <- allRes_f1.CeC_ordered[order(allRes_f1.CeC_ordered$Significant, decreasing = T),]


f1_CC_2.up <- ggplot(allRes_f1.CeC_ordered[c(1:7),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular compartement", 
       subtitle = "F1in vs F1out upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold")) + 
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))




png(paste(results.dir_noFDR, "CC_UP_top7.png", sep = "/"), width = 1000, height = 600)
f1_CC_2.up
dev.off()



###Downregulated genes ----

#get gene IDs for the enrichement 
genes_Down_f1 <- read.csv("DESeq2_tables/F1in_vs_F1out_sigDOWN.csv", header = TRUE)[,2]

#Get annotation for upregulated genes 
Down.annotated <- f1_sub_annot[genes_Down_f1,]

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_Down_f1))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_f1_Down.C <- new("topGOdata", ontology = "CC", 
                               description = "f1 enrichment using kallisto quantification and dammit annotation",
                               allGenes = compared_genes, 
                               annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                               gene2GO = f1_GOCeC,
                               nodeSize = 5) 

#Run fisher test
resultFisher_f1_Down.C <- runTest(GoData_f1_Down.C, 
                                    algorithm = "weight01", 
                                    statistic = "fisher")
resultFisher_f1_Down.C

#819 GO terms scored, 3 with p < 0.01
#(347 out of 509 significant genes, and 14314 feasible genes, or gene universe)

allRes_f1_Down.CeC <- GenTable(GoData_f1_Down.C, 
                                         classicFisher = resultFisher_f1_Down.C,
                                         orderBy= "classicFisher", 
                                         numChar = 100, 
                                         topNodes = length(resultFisher_f1_Down.C@score))



#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the upregulated genes 
#were removed 
allRes_f1_Down.CeC_filtered <- allRes_f1_Down.CeC[which(allRes_f1_Down.CeC$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_f1_Down.CeC_filtered$classicFisher, method = "fdr")

allRes_f1_Down.CeC_filtered$padj <- adj.pvalues

write.xlsx(allRes_f1_Down.CeC_filtered, file = paste(results.dir_2,"CeC_Down_resTable.xlsx", sep = "/"), row.names = F) 

n_term_sig <- sum(allRes_f1_Down.CeC_filtered$padj <= 0.05)

#Plot the sig Go terms 
f1_CC.down <- ggplot(allRes_f1_Down.CeC_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular compartement", 
       subtitle = "f1 downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.05, "CC_Down_sig_FDR_0.05.png", sep = "/"), width = 1100, height = 600)
f1_CC.down
dev.off()


n_term_sig <- sum(allRes_f1_Down.CeC_filtered$padj <= 0.1)

#Plot the sig Go terms 
f1_CC.down <- ggplot(allRes_f1_Down.CeC_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular compartement", 
       subtitle = "f1 downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))




png(paste(results.dir_2.1, "CC_Down_sig_FDR_0.1.png", sep = "/"), width = 1100, height = 600)
f1_CC.down
dev.off()


allRes_f1.CeC_ordered <- allRes_f1_Down.CeC_filtered[which(as.numeric(allRes_f1_Down.CeC_filtered$classicFisher) < 0.01),]
allRes_f1.CeC_ordered <- allRes_f1.CeC_ordered[order(allRes_f1.CeC_ordered$Significant, decreasing = T),]


f1_CC_2.down <- ggplot(allRes_f1.CeC_ordered[c(1:4),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular function", 
       subtitle = "F1in vs F1out downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_noFDR, "CC_Down_top4.png", sep = "/"), width = 1000, height = 600)
f1_CC_2.down
dev.off()


###Diff expressed ---- 
#(over representation of terms that are in both under and over represented)

#concatenate datasets from up and down annotated 
diff_f1.CC <- rbind(UP.annotated, Down.annotated)
genes_diff_expressed <- row.names(diff_f1.CC)

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_diff_expressed))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_f1_diff_expressed.C <- new("topGOdata", ontology = "CC", 
                                  description = "f1in_vs_f1out enrichment using kallisto quantification and dammit annotation",
                                  allGenes = compared_genes, 
                                  annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                                  gene2GO = f1_GOCeC,
                                  nodeSize = 5) 

#Run fisher test
resultFisher_f1_diff_expressed.C <- runTest(GoData_f1_diff_expressed.C, 
                                            algorithm = "weight01", 
                                            statistic = "fisher")
resultFisher_f1_diff_expressed.C

#819 GO terms scored, 6 with p < 0.01
#(780 out of 1214 significant genes, and 14314 feasible genes, or gene universe)

allRes_f1in_diff_expressed.CeC <- GenTable(GoData_f1_diff_expressed.C, 
                                           classicFisher = resultFisher_f1_diff_expressed.C,
                                           orderBy= "classicFisher", 
                                           numChar = 100, 
                                           topNodes = length(resultFisher_f1_diff_expressed.C@score))


#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the diff expressed genes 
#were removed 
allRes_f1in_diff_expressed.CeC_filtered <- allRes_f1in_diff_expressed.CeC[which(allRes_f1in_diff_expressed.CeC$Significant > 0),]

#correct for multiple testing with FDR
adj.Cvalues <- p.adjust(allRes_f1in_diff_expressed.CeC_filtered$classicFisher, method = "fdr")

allRes_f1in_diff_expressed.CeC_filtered$padj <- adj.Cvalues

write.xlsx(allRes_f1in_diff_expressed.CeC_filtered, file = paste(results.dir_2,"CeC_diff_expressed_resTable.xlsx", sep = "/"), row.names = F) 


n_term_sig <- sum(allRes_f1in_diff_expressed.CeC_filtered$padj <= 0.05)

#Plot the sig Go terms 
f1_CC.diff_expressed <- ggplot(allRes_f1in_diff_expressed.CeC_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular component", 
       subtitle = "f1in vs f1out differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.05, "CC_diff_expressed_sig_FDR_0.05.png", sep = "/"), width = 1100, height = 600)
f1_CC.diff_expressed
dev.off()

n_term_sig <- sum(allRes_f1in_diff_expressed.CeC_filtered$padj <= 0.1)

#Plot the sig Go terms 
f1_CC.diff_expressed <- ggplot(allRes_f1in_diff_expressed.CeC_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular component", 
       subtitle = "f1in vs f1out differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.1, "CC_diff_expressed_sig_FDR_0.1.png", sep = "/"), width = 1100, height = 600)
f1_CC.diff_expressed
dev.off()


allRes_f1in.CeC_ordered <- allRes_f1in_diff_expressed.CeC_filtered[which(as.numeric(allRes_f1in_diff_expressed.CeC_filtered$classicFisher) < 0.01),]
allRes_f1in.CeC_ordered <- allRes_f1in.CeC_ordered[order(allRes_f1in.CeC_ordered$Significant, decreasing = T),]


f1_CC_2.diff_expressed <- ggplot(allRes_f1in.CeC_ordered[c(1:6),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular component", 
       subtitle = "F1in vs F1out differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_noFDR, "CC_diff_expressed_6.png", sep = "/"), width = 1000, height = 600)
f1_CC_2.diff_expressed
dev.off()



############################################
#### F4in vs f4out enrichement analysis ####
############################################

#make results directory
results.dir_2 <- paste(results.dir, "F4in_vs_f4out", sep = "/")
if (!dir.exists(results.dir_2)){
  dir.create(results.dir_2)
}


results.dir_2.1 <- paste(results.dir_2, "FDR_0.1", sep = "/")
if (!dir.exists(results.dir_2.1)){
  dir.create(results.dir_2.1)
}

results.dir_2.05 <- paste(results.dir_2, "FDR_0.05", sep = "/")
if (!dir.exists(results.dir_2.05)){
  dir.create(results.dir_2.05)
}


results.dir_noFDR <- paste(results.dir_2, "Before_muliple_testing", sep = "/")
if (!dir.exists(results.dir_noFDR)){
  dir.create(results.dir_noFDR)
}

f4_bg_genes <- read.csv(file = "DESeq2_tables/F4in_vs_F4out_Res.csv")$X
f4_sub_annot <- annotation_file[f4_bg_genes,]


###Make backgroud annotation ---- 
#For each compariton, the background annotation is different as we only consider the gene with
#gene counts > 0 for the compariton. For this filtering, the genes in the res table from DesEq will be used, 
#as genes with rowcounts = 0 were filtered out before 


##Subset annotation files for backgroud annotation 
##Go(P) <- biological process
f4_Go.P <- data.frame("gene" = row.names(f4_sub_annot), "GO.BP" = f4_sub_annot[,9])
f4_Go.P$'GO.BP' <- str_replace_all(f4_Go.P$'GO.BP', ";", ",")
write.table(f4_Go.P, file = "f4_GoP.txt", row.names = FALSE, col.names = F,quote = F, sep = "\t")

##Go(F) <- molecular function
f4_Go.F <- data.frame("gene" = row.names(f4_sub_annot), "GO.MF" = f4_sub_annot[,11])
f4_Go.F$'GO.MF' <- str_replace_all(f4_Go.F$'GO.MF', ";", ",")
write.table(f4_Go.F, file = "f4_GoF.txt", row.names = FALSE, col.names = F,quote = F, sep = "\t")

##Go(C) <- celular compartement
f4_Go.C <- data.frame("gene" = row.names(f4_sub_annot), "GO.CC" = f4_sub_annot[,13])
f4_Go.C$'GO.CC' <- str_replace_all(f4_Go.C$'GO.CC', ";", ",")
write.table(f4_Go.C, file = "f4_GoC.txt", row.names = FALSE, col.names = F,quote = F, sep = "\t")


##By Biological process ----
#Get background annotation 
f4_GOByP <- readMappings(file = "f4_GoP.txt")
bg_genes <- names(f4_GOByP)


###Upregulated genes ----
#get gene IDs for the enrichement 
genes_UP_f4 <- read.csv("DESeq2_tables/f4in_vs_f4out_sigUP.csv", header = TRUE)[,2]

#Get annotation for upregulated genes 
UP.annotated <- f4_sub_annot[genes_UP_f4,]

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_UP_f4))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_f4_up.P <- new("topGOdata", ontology = "BP", 
                      description = "f4in vs f4out enrichment using kallisto quantification and dammit annotation",
                      allGenes = compared_genes, 
                      annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                      gene2GO = f4_GOByP,
                      nodeSize = 5) 

#Run fisher test
resultFisher_f4_up.P <- runTest(GoData_f4_up.P, 
                                algorithm = "weight01", 
                                statistic = "fisher")
resultFisher_f4_up.P

#4393 GO terms scored, 21 with p < 0.01
#(670 out of 1133 significant genes, and 13406 feasible genes, or gene universe)

allRes_f4_up.ByP <- GenTable(GoData_f4_up.P, 
                             classicFisher = resultFisher_f4_up.P,
                             orderBy= "classicFisher", 
                             numChar = 100, 
                             topNodes = length(resultFisher_f4_up.P@score))

#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the upregulated genes 
#were removed 
allRes_f4_up.ByP_filtered <- allRes_f4_up.ByP[which(allRes_f4_up.ByP$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_f4_up.ByP_filtered$classicFisher, method = "fdr")

allRes_f4_up.ByP_filtered$padj <- adj.pvalues

write.xlsx(allRes_f4_up.ByP_filtered, file = paste(results.dir_2,"byP_UP_resTable.xlsx", sep = "/"), row.names = F) 

n_term_sig <- sum(allRes_f4_up.ByP_filtered$padj <= 0.05)

#Plot the sig Go terms 
f4_BP.up <- ggplot(allRes_f4_up.ByP_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "f4in vs f4out upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.05, "BP_UP_sig_FDR_0.05.png", sep = "/"), width = 1000, height = 600)
f4_BP.up
dev.off()

n_term_sig <- sum(allRes_f4_up.ByP_filtered$padj <= 0.1)

#Plot the sig Go terms 
f4_BP.up <- ggplot(allRes_f4_up.ByP_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "f4in vs f4out upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.1, "BP_UP_sig_FDR_0.1.png", sep = "/"), width = 1000, height = 600)
f4_BP.up
dev.off()


allRes_f4.ByP_ordered <- allRes_f4_up.ByP_filtered[which(as.numeric(allRes_f4_up.ByP_filtered$classicFisher) < 0.01),]
allRes_f4.ByP_ordered <- allRes_f4.ByP_ordered[order(allRes_f4.ByP_ordered$Significant, decreasing = T),]


f4_BP_2.up <- ggplot(allRes_f4.ByP_ordered[c(1:20),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "f4in vs f4out upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_noFDR, "BP_UP_top20.png", sep = "/"), width = 1000, height = 600)
f4_BP_2.up
dev.off()



###Downregulated genes ----

#get gene IDs for the enrichement 
genes_Down_f4 <- read.csv("DESeq2_tables/f4in_vs_f4out_sigDOWN.csv", header = TRUE)[,2]

#Get annotation for upregulated genes 
Down.annotated <- f4_sub_annot[genes_Down_f4,]

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_Down_f4))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_f4_Down.P <- new("topGOdata", ontology = "BP", 
                             description = "f4 enrichment using kallisto quantification and dammit annotation",
                             allGenes = compared_genes, 
                             annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                             gene2GO = f4_GOByP,
                             nodeSize = 5) 

#Run fisher test
resultFisher_f4_Down.P <- runTest(GoData_f4_Down.P, 
                                  algorithm = "weight01", 
                                  statistic = "fisher")
resultFisher_f4_Down.P

#4394 GO terms scored, 19 with p < 0.01
#(304 out of 509 significant genes, and 13122 feasible genes, or gene universe)

allRes_f4_Down.ByP <- GenTable(GoData_f4_Down.P, 
                               classicFisher = resultFisher_f4_Down.P,
                               orderBy= "classicFisher", 
                               numChar = 100, 
                               topNodes = length(resultFisher_f4_Down.P@score))


#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the upregulated genes 
#were removed 
allRes_f4_Down.ByP_filtered <- allRes_f4_Down.ByP[which(allRes_f4_Down.ByP$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_f4_Down.ByP_filtered$classicFisher, method = "fdr")

allRes_f4_Down.ByP_filtered$padj <- adj.pvalues

write.xlsx(allRes_f4_Down.ByP_filtered, file = paste(results.dir_2,"byP_Down_resTable.xlsx", sep = "/"), row.names = F) 

n_term_sig <- sum(allRes_f4_Down.ByP_filtered$padj <= 0.05)

#Plot the sig Go terms 
f4_BP.down <- ggplot(allRes_f4_Down.ByP_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "f4in vs f4out downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.05, "BP_Down_sig_FDR_0.05.png", sep = "/"), width = 1100, height = 600)
f4_BP.down
dev.off()


n_term_sig <- sum(allRes_f4_Down.ByP_filtered$padj <= 0.1)

#Plot the sig Go terms 
f4_BP.down <- ggplot(allRes_f4_Down.ByP_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "f4in vs f4out downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.1, "BP_Down_sig_FDR_0.1.png", sep = "/"), width = 1100, height = 600)
f4_BP.down
dev.off()


allRes_f4.ByP_ordered <- allRes_f4_Down.ByP_filtered[which(as.numeric(allRes_f4_Down.ByP_filtered$classicFisher) < 0.01),]
allRes_f4.ByP_ordered <- allRes_f4.ByP_ordered[order(allRes_f4.ByP_ordered$Significant, decreasing = T),]


f4_BP_2.down <- ggplot(allRes_f4.ByP_ordered[c(1:20),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "f4in vs f4out downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))




png(paste(results.dir_noFDR, "BP_Down_top20.png", sep = "/"), width = 1000, height = 600)
f4_BP_2.down
dev.off()


###Diff expressed ---- 
#(over representation of terms that are in both under and over represented)

#concatenate datasets from up and down annotated 
diff_f4.BP <- rbind(UP.annotated, Down.annotated)
genes_diff_expressed <- row.names(diff_f4.BP)

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_diff_expressed))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_f4_diff_expressed.P <- new("topGOdata", ontology = "BP", 
                                  description = "f4in_vs_f4out enrichment using kallisto quantification and dammit annotation",
                                  allGenes = compared_genes, 
                                  annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                                  gene2GO = f4_GOByP,
                                  nodeSize = 5) 

#Run fisher test
resultFisher_f4_diff_expressed.P <- runTest(GoData_f4_diff_expressed.P, 
                                            algorithm = "weight01", 
                                            statistic = "fisher")
resultFisher_f4_diff_expressed.P

#4498 GO terms scored,  with p < 0.01
#( out of  significant genes, and 13122 feasible genes, or gene universe)

allRes_f4in_diff_expressed.ByP <- GenTable(GoData_f4_diff_expressed.P, 
                                           classicFisher = resultFisher_f4_diff_expressed.P,
                                           orderBy= "classicFisher", 
                                           numChar = 100, 
                                           topNodes = length(resultFisher_f4_diff_expressed.P@score))


#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the diff expressed genes 
#were removed 
allRes_f4in_diff_expressed.ByP_filtered <- allRes_f4in_diff_expressed.ByP[which(allRes_f4in_diff_expressed.ByP$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_f4in_diff_expressed.ByP_filtered$classicFisher, method = "fdr")

allRes_f4in_diff_expressed.ByP_filtered$padj <- adj.pvalues

write.xlsx(allRes_f4in_diff_expressed.ByP_filtered, file = paste(results.dir_2,"byP_diff_expressed_resTable.xlsx", sep = "/"), row.names = F) 

n_term_sig <- sum(allRes_f4in_diff_expressed.ByP_filtered$padj <= 0.05)

#Plot the sig Go terms 
f4_BP.diff_expressed <- ggplot(allRes_f4in_diff_expressed.ByP_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "f4in vs f4out differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.05, "BP_diff_expressed_sig_FDR_0.05.png", sep = "/"), width = 1100, height = 600)
f4_BP.diff_expressed
dev.off()



n_term_sig <- sum(allRes_f4in_diff_expressed.ByP_filtered$padj <= 0.1)

#Plot the sig Go terms 
f4_BP.diff_expressed <- ggplot(allRes_f4in_diff_expressed.ByP_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "f4in vs f4out differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.1, "BP_diff_expressed_sig_FDR_0.1.png", sep = "/"), width = 1100, height = 600)
f4_BP.diff_expressed
dev.off()

allRes_f4in.ByP_ordered <- allRes_f4in_diff_expressed.ByP_filtered[which(as.numeric(allRes_f4in_diff_expressed.ByP_filtered$classicFisher) < 0.01),]
allRes_f4in.ByP_ordered <- allRes_f4in.ByP_ordered[order(allRes_f4in.ByP_ordered$Significant, decreasing = T),]


f4_BP_2.diff_expressed <- ggplot(allRes_f4in.ByP_ordered[c(1:20),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Biological Process", 
       subtitle = "f4in vs f4out differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_noFDR, "BP_diff_expressed_top20.png", sep = "/"), width = 1100, height = 600)
f4_BP_2.diff_expressed
dev.off()


##By Molecular function ----
#Get background annotation 
f4_GOMoF <- readMappings(file = "f4_GoF.txt")
bg_genes <- names(f4_GOMoF)


###Upregulated genes ----
#get gene IDs for the enrichement 
genes_UP_f4 <- read.csv("DESeq2_tables/f4in_vs_f4out_sigUP.csv", header = TRUE)[,2]

#Get annotation for upregulated genes 
UP.annotated <- f4_sub_annot[genes_UP_f4,]

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_UP_f4))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_f4_up.F <- new("topGOdata", ontology = "MF", 
                      description = "f4 enrichment using kallisto quantification and dammit annotation",
                      allGenes = compared_genes, 
                      annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                      gene2GO = f4_GOMoF,
                      nodeSize = 5) 

#Run fisher test
resultFisher_f4_up.F <- runTest(GoData_f4_up.F, 
                                algorithm = "weight01", 
                                statistic = "fisher")
resultFisher_f4_up.F

#1317 GO terms scored, 16 with p < 0.01
#(454 out of 705 significant genes, and 13470 feasible genes, or gene universe)

allRes_f4_up.MoF <- GenTable(GoData_f4_up.F, 
                             classicFisher = resultFisher_f4_up.F,
                             orderBy= "classicFisher", 
                             numChar = 100, 
                             topNodes = length(resultFisher_f4_up.F@score))

#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the upregulated genes 
#were removed 
allRes_f4_up.MoF_filtered <- allRes_f4_up.MoF[which(allRes_f4_up.MoF$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_f4_up.MoF_filtered$classicFisher, method = "fdr")

allRes_f4_up.MoF_filtered$padj <- adj.pvalues

write.xlsx(allRes_f4_up.MoF_filtered, file = paste(results.dir_2,"MoF_UP_resTable.xlsx", sep = "/"), row.names = F) 

n_term_sig <- sum(allRes_f4_up.MoF_filtered$padj <= 0.05)

#Plot the sig Go terms 
f4_MF.up <- ggplot(allRes_f4_up.MoF_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular Function", 
       subtitle = "f4 upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.05, "MF_UP_sig_FDR_0.05.png", sep = "/"), width = 1200, height = 600)
f4_MF.up
dev.off()


n_term_sig <- sum(allRes_f4_up.MoF_filtered$padj <= 0.1)

#Plot the sig Go terms 
f4_MF.up <- ggplot(allRes_f4_up.MoF_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular Function", 
       subtitle = "f4 upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.1, "MF_UP_sig_FDR_0.1.png", sep = "/"), width = 1200, height = 600)
f4_MF.up
dev.off()


allRes_f4.MoF_ordered <- allRes_f4_up.MoF_filtered[which(as.numeric(allRes_f4_up.MoF_filtered$classicFisher) < 0.01),]
allRes_f4.MoF_ordered <- allRes_f4.MoF_ordered[order(allRes_f4.MoF_ordered$Significant, decreasing = T),]

allRes_f4.MoF_ordered$Term[3] <- "oxidoreductase activity, acting on paired donors,\nwith incorporation or reduction of molecular oxygen"

f4_MF_2.up <- ggplot(allRes_f4.MoF_ordered[c(1:18),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular function", 
       subtitle = "f4in vs f4out upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_noFDR, "MF_UP_top18.png", sep = "/"), width = 1000, height = 600)
f4_MF_2.up
dev.off()



###Downregulated genes ----

#get gene IDs for the enrichement 
genes_Down_f4 <- read.csv("DESeq2_tables/f4in_vs_f4out_sigDOWN.csv", header = TRUE)[,2]

#Get annotation for upregulated genes 
Down.annotated <- f4_sub_annot[genes_Down_f4,]

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_Down_f4))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_f4_Down.F <- new("topGOdata", ontology = "MF", 
                        description = "f4 enrichment using kallisto quantification and dammit annotation",
                        allGenes = compared_genes, 
                        annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                        gene2GO = f4_GOMoF,
                        nodeSize = 5) 

#Run fisher test
resultFisher_f4_Down.F <- runTest(GoData_f4_Down.F, 
                                  algorithm = "weight01", 
                                  statistic = "fisher")
resultFisher_f4_Down.F

#1317 GO terms scored, 15 with p < 0.01
#(336  out of 509 significant genes, and 13406 feasible genes, or gene universe)

allRes_f4_Down.MoF <- GenTable(GoData_f4_Down.F, 
                               classicFisher = resultFisher_f4_Down.F,
                               orderBy= "classicFisher", 
                               numChar = 100, 
                               topNodes = length(resultFisher_f4_Down.F@score))


#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the upregulated genes 
#were removed 
allRes_f4_Down.MoF_filtered <- allRes_f4_Down.MoF[which(allRes_f4_Down.MoF$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_f4_Down.MoF_filtered$classicFisher, method = "fdr")

allRes_f4_Down.MoF_filtered$padj <- adj.pvalues

write.xlsx(allRes_f4_Down.MoF_filtered, file = paste(results.dir_2,"MoF_Down_resTable.xlsx", sep = "/"), row.names = F) 

n_term_sig <- sum(allRes_f4_Down.MoF_filtered$padj <= 0.05)

#Plot the sig Go terms 
f4_MF.down <- ggplot(allRes_f4_Down.MoF_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular function", 
       subtitle = "f4 downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.05, "MF_Down_sig_FDR_0.05.png", sep = "/"), width = 1100, height = 600)
f4_MF.down
dev.off()

n_term_sig <- sum(allRes_f4_Down.MoF_filtered$padj <= 0.1)

#Plot the sig Go terms 
f4_MF.down <- ggplot(allRes_f4_Down.MoF_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular function", 
       subtitle = "f4 downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.1, "MF_Down_sig_FDR_0.1.png", sep = "/"), width = 1100, height = 600)
f4_MF.down
dev.off()

allRes_f4.MoF_ordered <- allRes_f4_Down.MoF_filtered[which(as.numeric(allRes_f4_Down.MoF_filtered$classicFisher) < 0.01),]
allRes_f4.MoF_ordered <- allRes_f4.MoF_ordered[order(allRes_f4.MoF_ordered$Significant, decreasing = T),]

allRes_f4.MoF_ordered$Term[6] <- "oxidoreductase activity, acting on paired donors,\nwith incorporation or reduction of molecular oxygen"

f4_MF_2.down <- ggplot(allRes_f4.MoF_ordered[c(1:20),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular function", 
       subtitle = "f4in vs f4out downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_noFDR, "MF_Down_top20.png", sep = "/"), width = 1100, height = 600)
f4_MF_2.down
dev.off()

###Dif expressed ---- 
#(over representation of terms that are in both under and over represented)

#concatenate datasets from up and down annotated 
diff_f4.MF <- rbind(UP.annotated, Down.annotated)
genes_diff_expressed <- row.names(diff_f4.MF)

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_diff_expressed))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_f4_diff_expressed.M <- new("topGOdata", ontology = "MF", 
                                  description = "f4in_vs_f4out enrichment using kallisto quantification and dammit annotation",
                                  allGenes = compared_genes, 
                                  annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                                  gene2GO = f4_GOMoF,
                                  nodeSize = 5) 

#Run fisher test
resultFisher_f4_diff_expressed.M <- runTest(GoData_f4_diff_expressed.M, 
                                            algorithm = "weight01", 
                                            statistic = "fisher")
resultFisher_f4_diff_expressed.M

#1317 GO terms scored, 27 with p < 0.01
#(790 out of 1214 significant genes, and 13406 feasible genes, or gene universe)

allRes_f4in_diff_expressed.MoF <- GenTable(GoData_f4_diff_expressed.M, 
                                           classicFisher = resultFisher_f4_diff_expressed.M,
                                           orderBy= "classicFisher", 
                                           numChar = 100, 
                                           topNodes = length(resultFisher_f4_diff_expressed.M@score))


#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the diff expressed genes 
#were removed 
allRes_f4in_diff_expressed.MoF_filtered <- allRes_f4in_diff_expressed.MoF[which(allRes_f4in_diff_expressed.MoF$Significant > 0),]

#correct for multiple testing with FDR
adj.Mvalues <- p.adjust(allRes_f4in_diff_expressed.MoF_filtered$classicFisher, method = "fdr")

allRes_f4in_diff_expressed.MoF_filtered$padj <- adj.Mvalues

write.xlsx(allRes_f4in_diff_expressed.MoF_filtered, file = paste(results.dir_2,"MoF_diff_expressed_resTable.xlsx", sep = "/"), row.names = F) 


n_term_sig <- sum(allRes_f4in_diff_expressed.MoF_filtered$padj <= 0.05)

#Plot the sig Go terms 
f4_MF.diff_expressed <- ggplot(allRes_f4in_diff_expressed.MoF_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular function", 
       subtitle = "f4in vs f4out differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.05, "MF_diff_expressed_sig_FDR_0.05.png", sep = "/"), width = 1100, height = 600)
f4_MF.diff_expressed
dev.off()

n_term_sig <- sum(allRes_f4in_diff_expressed.MoF_filtered$padj <= 0.1)

#Plot the sig Go terms 
f4_MF.diff_expressed <- ggplot(allRes_f4in_diff_expressed.MoF_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular function", 
       subtitle = "f4in vs f4out differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.1, "MF_diff_expressed_sig_FDR_0.1.png", sep = "/"), width = 1100, height = 600)
f4_MF.diff_expressed
dev.off()


allRes_f4in.MoF_ordered <- allRes_f4in_diff_expressed.MoF_filtered[which(as.numeric(allRes_f4in_diff_expressed.MoF_filtered$classicFisher) < 0.01),]
allRes_f4in.MoF_ordered <- allRes_f4in.MoF_ordered[order(allRes_f4in.MoF_ordered$Significant, decreasing = T),]

allRes_f4in.MoF_ordered$Term[5] <- "oxidoreductase activity, acting on paired donors,\nwith incorporation or reduction of molecular oxygen"

f4_MF_2.diff_expressed <- ggplot(allRes_f4in.MoF_ordered[c(1:20),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Molecular function", 
       subtitle = "f4in vs f4out differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_noFDR, "MF_diff_expressed_top20.png", sep = "/"), width = 1100, height = 600)
f4_MF_2.diff_expressed
dev.off()




##By cellular_component ----
#Get background annotation 
f4_GOCeC <- readMappings(file = "f4_GoC.txt")
bg_genes <- names(f4_GOCeC)


###Upregulated genes ----
#get gene IDs for the enrichement 
genes_UP_f4 <- read.csv("DESeq2_tables/f4in_vs_f4out_sigUP.csv", header = TRUE)[,2]

#Get annotation for upregulated genes 
UP.annotated <- f4_sub_annot[genes_UP_f4,]

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_UP_f4))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_f4_up.C <- new("topGOdata", ontology = "CC", 
                      description = "f4 enrichment using kallisto quantification and dammit annotation",
                      allGenes = compared_genes, 
                      annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                      gene2GO = f4_GOCeC,
                      nodeSize = 5) 

#Run fisher test
resultFisher_f4_up.C <- runTest(GoData_f4_up.C, 
                                algorithm = "weight01", 
                                statistic = "fisher")
resultFisher_f4_up.C

#823 GO terms scored, 6 with p < 0.01
#(707 out of 1133 significant genes, and 13406 feasible genes, or gene universe)

allRes_f4_up.CeC <- GenTable(GoData_f4_up.C, 
                             classicFisher = resultFisher_f4_up.C,
                             orderBy= "classicFisher", 
                             numChar = 100, 
                             topNodes = length(resultFisher_f4_up.C@score))

#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the upregulated genes 
#were removed 
allRes_f4_up.CeC_filtered <- allRes_f4_up.CeC[which(allRes_f4_up.CeC$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_f4_up.CeC_filtered$classicFisher, method = "fdr")

allRes_f4_up.CeC_filtered$padj <- adj.pvalues

write.xlsx(allRes_f4_up.CeC_filtered, file = paste(results.dir_2,"CeC_UP_resTable.xlsx", sep = "/"), row.names = F) 


n_term_sig <- sum(allRes_f4_up.CeC_filtered$padj <= 0.05)

#Plot the sig Go terms 
f4_CC.up <- ggplot(allRes_f4_up.CeC_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular compartement", 
       subtitle = "f4 upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.05, "CC_UP_sig_FDR_0.05.png", sep = "/"), width = 1000, height = 600)
f4_CC.up
dev.off()

n_term_sig <- sum(allRes_f4_up.CeC_filtered$padj <= 0.1)

#Plot the sig Go terms 
f4_CC.up <- ggplot(allRes_f4_up.CeC_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular compartement", 
       subtitle = "f4 upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.1, "CC_UP_sig_FDR_0.1.png", sep = "/"), width = 1000, height = 600)
f4_CC.up
dev.off()

allRes_f4.CeC_ordered <- allRes_f4_up.CeC_filtered[which(as.numeric(allRes_f4_up.CeC_filtered$classicFisher) < 0.01),]
allRes_f4.CeC_ordered <- allRes_f4.CeC_ordered[order(allRes_f4.CeC_ordered$Significant, decreasing = T),]


f4_CC_2.up <- ggplot(allRes_f4.CeC_ordered[c(1:6),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "royalblue") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular compartement", 
       subtitle = "f4in vs f4out upregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_noFDR, "CC_UP_top6.png", sep = "/"), width = 1000, height = 600)
f4_CC_2.up
dev.off()



###Downregulated genes ----

#get gene IDs for the enrichement 
genes_Down_f4 <- read.csv("DESeq2_tables/f4in_vs_f4out_sigDOWN.csv", header = TRUE)[,2]

#Get annotation for upregulated genes 
Down.annotated <- f4_sub_annot[genes_Down_f4,]

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_Down_f4))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_f4_Down.C <- new("topGOdata", ontology = "CC", 
                        description = "f4 enrichment using kallisto quantification and dammit annotation",
                        allGenes = compared_genes, 
                        annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                        gene2GO = f4_GOCeC,
                        nodeSize = 5) 

#Run fisher test
resultFisher_f4_Down.C <- runTest(GoData_f4_Down.C, 
                                  algorithm = "weight01", 
                                  statistic = "fisher")
resultFisher_f4_Down.C

#819 GO terms scored, 3 with p < 0.01
#(1188 out of 1667 significant genes, and 14314 feasible genes, or gene universe)

allRes_f4_Down.CeC <- GenTable(GoData_f4_Down.C, 
                               classicFisher = resultFisher_f4_Down.C,
                               orderBy= "classicFisher", 
                               numChar = 100, 
                               topNodes = length(resultFisher_f4_Down.C@score))



#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the upregulated genes 
#were removed 
allRes_f4_Down.CeC_filtered <- allRes_f4_Down.CeC[which(allRes_f4_Down.CeC$Significant > 0),]

#correct for multiple testing with FDR
adj.pvalues <- p.adjust(allRes_f4_Down.CeC_filtered$classicFisher, method = "fdr")

allRes_f4_Down.CeC_filtered$padj <- adj.pvalues

write.xlsx(allRes_f4_Down.CeC_filtered, file = paste(results.dir_2,"CeC_Down_resTable.xlsx", sep = "/"), row.names = F) 


n_term_sig <- sum(allRes_f4_Down.CeC_filtered$padj <= 0.05)

#Plot the sig Go terms 
f4_CC.down <- ggplot(allRes_f4_Down.CeC_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular compartement", 
       subtitle = "f4 downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.05, "CC_Down_sig_FDR_0.05.png", sep = "/"), width = 1100, height = 600)
f4_CC.down
dev.off()


n_term_sig <- sum(allRes_f4_Down.CeC_filtered$padj <= 0.1)

#Plot the sig Go terms 
f4_CC.down <- ggplot(allRes_f4_Down.CeC_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular compartement", 
       subtitle = "f4 downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_2.1, "CC_Down_sig_FDR_0.1.png", sep = "/"), width = 1100, height = 600)
f4_CC.down
dev.off()

allRes_f4.CeC_ordered <- allRes_f4_Down.CeC_filtered[which(as.numeric(allRes_f4_Down.CeC_filtered$classicFisher) < 0.01),]
allRes_f4.CeC_ordered <- allRes_f4.CeC_ordered[order(allRes_f4.CeC_ordered$Significant, decreasing = T),]


f4_CC_2.down <- ggplot(allRes_f4.CeC_ordered[c(1:12),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "red") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular function", 
       subtitle = "f4in vs f4out downregulated", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() + 
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))


png(paste(results.dir_noFDR, "CC_Down_top12.png", sep = "/"), width = 1000, height = 600)
f4_CC_2.down
dev.off()


###Dif expressed ---- 
#(over representation of terms that are in both under and over represented)

#concatenate datasets from up and down annotated 
diff_f4.CC <- rbind(UP.annotated, Down.annotated)
genes_diff_expressed <- row.names(diff_f4.CC)

#get interesting genes
compared_genes <- factor(as.integer(bg_genes %in% genes_diff_expressed))
names(compared_genes) <- bg_genes

###NOTE: Node size = 5 will be used to prune GO hierarchy from the terms 
#which have less than 5 annotated genes, as these lead to more stable results 

GoData_f4_diff_expressed.C <- new("topGOdata", ontology = "CC", 
                                  description = "f4in_vs_f4out enrichment using kallisto quantification and dammit annotation",
                                  allGenes = compared_genes, 
                                  annot = annFUN.gene2GO, #this argument specifies which annotation type will be used 
                                  gene2GO = f4_GOCeC,
                                  nodeSize = 5) 

#Run fisher test
resultFisher_f4_diff_expressed.C <- runTest(GoData_f4_diff_expressed.C, 
                                            algorithm = "weight01", 
                                            statistic = "fisher")
resultFisher_f4_diff_expressed.C

#819 GO terms scored, 6 with p < 0.01
#(780 out of 1214 significant genes, and 14314 feasible genes, or gene universe)

allRes_f4in_diff_expressed.CeC <- GenTable(GoData_f4_diff_expressed.C, 
                                           classicFisher = resultFisher_f4_diff_expressed.C,
                                           orderBy= "classicFisher", 
                                           numChar = 100, 
                                           topNodes = length(resultFisher_f4_diff_expressed.C@score))


#correct for multiple testing 
#Before correcting the p.values, GO terms that do not appear in the diff expressed genes 
#were removed 
allRes_f4in_diff_expressed.CeC_filtered <- allRes_f4in_diff_expressed.CeC[which(allRes_f4in_diff_expressed.CeC$Significant > 0),]

#correct for multiple testing with FDR
adj.Cvalues <- p.adjust(allRes_f4in_diff_expressed.CeC_filtered$classicFisher, method = "fdr")

allRes_f4in_diff_expressed.CeC_filtered$padj <- adj.Cvalues

write.xlsx(allRes_f4in_diff_expressed.CeC_filtered, file = paste(results.dir_2,"CeC_diff_expressed_resTable.xlsx", sep = "/"), row.names = F) 


n_term_sig <- sum(allRes_f4in_diff_expressed.CeC_filtered$padj <= 0.05)

#Plot the sig Go terms 
f4_CC.diff_expressed <- ggplot(allRes_f4in_diff_expressed.CeC_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular component", 
       subtitle = "f4in vs f4out differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.05, "CC_diff_expressed_sig_0.05_FDR.png", sep = "/"), width = 1100, height = 600)
f4_CC.diff_expressed
dev.off()


n_term_sig <- sum(allRes_f4in_diff_expressed.CeC_filtered$padj <= 0.1)

#Plot the sig Go terms 
f4_CC.diff_expressed <- ggplot(allRes_f4in_diff_expressed.CeC_filtered[c(1:n_term_sig),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular component", 
       subtitle = "f4in vs f4out differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_2.1, "CC_diff_expressed_sig_0.1_FDR.png", sep = "/"), width = 1100, height = 600)
f4_CC.diff_expressed
dev.off()

allRes_f4in.CeC_ordered <- allRes_f4in_diff_expressed.CeC_filtered[which(as.numeric(allRes_f4in_diff_expressed.CeC_filtered$classicFisher) < 0.01),]
allRes_f4in.CeC_ordered <- allRes_f4in.CeC_ordered[order(allRes_f4in.CeC_ordered$Significant, decreasing = T),]


f4_CC_2.diff_expressed <- ggplot(allRes_f4in.CeC_ordered[c(1:10),], aes(x=reorder(Term, Significant), y = Significant)) +
  geom_bar(position="dodge",stat="identity", width= 0.6, fill = "gray") + coord_flip() +
  geom_text(aes(label=paste(round(Significant/Annotated*100), "%")), hjust = -0.4) + 
  labs(title = "Cellular component", 
       subtitle = "f4in vs f4out differential expressed", 
       x = "GO Term", y = "Significant genes number") + 
  theme_classic() +   
  theme(axis.text.x =element_text(size=14),  
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.title = element_text(size=14,face="bold"), 
        title = element_text(size=16,face="bold"))



png(paste(results.dir_noFDR, "CC_diff_expressed_10.png", sep = "/"), width = 1000, height = 600)
f4_CC_2.diff_expressed
dev.off()


##Make summary_table----
gene_set <- c("meso vs idio Upregulated", "meso vs idio downregulated", "meso vs idio diff expressed", 
              "meso vs idio Upregulated", "meso vs idio downregulated", "meso vs idio diff expressed", 
              "meso vs idio Upregulated", "meso vs idio downregulated", "meso vs idio diff expressed", 
              "F1in vs F1out Upregulated", "F1in vs F1out downregulated", "F1in vs F1out diff expressed", 
              "F1in vs F1out Upregulated", "F1in vs F1out downregulated", "F1in vs F1out diff expressed", 
              "F1in vs F1out Upregulated", "F1in vs F1out downregulated", "F1in vs F1out diff expressed", 
              "F4in vs F4out Upregulated", "F4in vs F4out downregulated", "F4in vs F4out diff expressed", 
              "F4in vs F4out Upregulated", "F4in vs F4out downregulated", "F4in vs F4out diff expressed", 
              "F4in vs F4out Upregulated", "F4in vs F4out downregulated", "F4in vs F4out diff expressed")

Ontology <- rep(c(rep("BP", 3), rep("MF", 3), rep("CC", 3)),3)

Significant_genes <- c(rep(c(length(genes_UP_meso_vs_idio), length(genes_Down_meso_vs_idio), length(c(genes_UP_meso_vs_idio, genes_Down_meso_vs_idio))), 3),
                       rep(c(length(genes_UP_f1), length(genes_Down_f1), length(c(genes_UP_f1, genes_Down_f1))), 3), 
                       rep(c(length(genes_UP_f4), length(genes_Down_f4), length(c(genes_UP_f4, genes_Down_f4))), 3)) 

#feasible genes
feasible_genes <- c(resultFisher_idio_up.P@geneData[2], resultFisher_idio_Down.P@geneData[2], resultFisher_idio_diff_expressed.P@geneData[2], 
                    resultFisher_idio_up.F@geneData[2], resultFisher_idio_Down.F@geneData[2], resultFisher_idio_diff_expressed.F@geneData[2],
                    resultFisher_idio_up.C@geneData[2], resultFisher_idio_Down.C@geneData[2], resultFisher_idio_diff_expressed.C@geneData[2], 
                    resultFisher_f1_up.P@geneData[2], resultFisher_f1_Down.P@geneData[2], resultFisher_f1_diff_expressed.P@geneData[2], 
                    resultFisher_f1_up.F@geneData[2], resultFisher_f1_Down.F@geneData[2], resultFisher_f1_diff_expressed.M@geneData[2],
                    resultFisher_f1_up.C@geneData[2], resultFisher_f1_Down.C@geneData[2], resultFisher_f1_diff_expressed.C@geneData[2], 
                    resultFisher_f4_up.P@geneData[2], resultFisher_f4_Down.P@geneData[2], resultFisher_f4_diff_expressed.P@geneData[2], 
                    resultFisher_f4_up.F@geneData[2], resultFisher_f4_Down.F@geneData[2], resultFisher_f4_diff_expressed.M@geneData[2],
                    resultFisher_f4_up.C@geneData[2], resultFisher_f4_Down.C@geneData[2], resultFisher_f4_diff_expressed.C@geneData[2])



#GO terms scored
GO_terms_score <- c(rep(length(allRes_meso_vs_idio_up.ByP$GO.ID), 3), rep(length(allRes_meso_vs_idio_up.MoF$GO.ID), 3), rep(length(allRes_meso_vs_idio_up.CeC$GO.ID), 3), 
                    rep(length(allRes_f1_up.ByP$GO.ID), 3), rep(length(allRes_f1_up.MoF$GO.ID), 3), rep(length(allRes_f1_up.CeC$GO.ID), 3),
                    rep(length(allRes_f4_up.ByP$GO.ID), 3), rep(length(allRes_f4_up.MoF$GO.ID), 3), rep(length(allRes_f4_up.CeC$GO.ID), 3))

#go terms significant at 0.01 (before multiple testing)

Go_terms_less_0.01 <- c(length(which(resultFisher_idio_up.P@score <= 0.01)), length(which(resultFisher_idio_Down.P@score <= 0.01)), length(which(resultFisher_idio_diff_expressed.P@score <= 0.01)), 
                        length(which(resultFisher_idio_up.F@score <= 0.01)), length(which(resultFisher_idio_Down.F@score <= 0.01)), length(which(resultFisher_idio_diff_expressed.F@score <= 0.01)), 
                        length(which(resultFisher_idio_up.C@score <= 0.01)), length(which(resultFisher_idio_Down.C@score <= 0.01)), length(which(resultFisher_idio_diff_expressed.C@score <= 0.01)),
                        length(which(resultFisher_f1_up.P@score <= 0.01)), length(which(resultFisher_f1_Down.P@score <= 0.01)), length(which(resultFisher_f1_diff_expressed.P@score <= 0.01)), 
                        length(which(resultFisher_f1_up.F@score <= 0.01)), length(which(resultFisher_f1_Down.F@score <= 0.01)), length(which(resultFisher_f1_diff_expressed.M@score <= 0.01)), 
                        length(which(resultFisher_f1_up.C@score <= 0.01)), length(which(resultFisher_f1_Down.C@score <= 0.01)), length(which(resultFisher_f1_diff_expressed.C@score <= 0.01)), 
                        length(which(resultFisher_f4_up.P@score <= 0.01)), length(which(resultFisher_f4_Down.P@score <= 0.01)), length(which(resultFisher_f4_diff_expressed.P@score <= 0.01)), 
                        length(which(resultFisher_f4_up.F@score <= 0.01)), length(which(resultFisher_f4_Down.F@score <= 0.01)), length(which(resultFisher_f4_diff_expressed.M@score <= 0.01)), 
                        length(which(resultFisher_f4_up.C@score <= 0.01)), length(which(resultFisher_f4_Down.C@score <= 0.01)), length(which(resultFisher_f4_diff_expressed.C@score <= 0.01)))


#go terms significant at 0.1 and 0.05 after FDR 

Go_terms_less_0.1_FDR <- c(length(which(allRes_meso_vs_idio_up.ByP_filtered$padj <= 0.1)),length(which(allRes_meso_vs_idio_Down.ByP_filtered$padj <= 0.1)), length(which(allRes_meso_vs_idio_diff_expressed.ByP_filtered$padj <= 0.1)), 
                           length(which(allRes_meso_vs_idio_up.MoF_filtered$padj <= 0.1)),length(which(allRes_meso_vs_idio_Down.MoF_filtered$padj <= 0.1)), length(which(allRes_meso_vs_idio_diff_expressed.MoF_filtered$padj <= 0.1)),
                           length(which(allRes_meso_vs_idio_up.CeC_filtered$padj <= 0.1)),length(which(allRes_meso_vs_idio_Down.CeC_filtered$padj <= 0.1)), length(which(allRes_meso_vs_idio_diff_expressed.CeC_filtered$padj <= 0.1)), 
                           length(which(allRes_f1_up.ByP_filtered$padj <= 0.1)),length(which(allRes_f1_Down.ByP_filtered$padj <= 0.1)), length(which(allRes_f1in_diff_expressed.ByP_filtered$padj <= 0.1)), 
                           length(which(allRes_f1_up.MoF_filtered$padj <= 0.1)),length(which(allRes_f1_Down.MoF_filtered$padj <= 0.1)), length(which(allRes_f1in_diff_expressed.MoF_filtered$padj <= 0.1)),
                           length(which(allRes_f1_up.CeC_filtered$padj <= 0.1)),length(which(allRes_f1_Down.CeC_filtered$padj <= 0.1)), length(which(allRes_f1in_diff_expressed.CeC_filtered$padj <= 0.1)), 
                           length(which(allRes_f4_up.ByP_filtered$padj <= 0.1)),length(which(allRes_f4_Down.ByP_filtered$padj <= 0.1)), length(which(allRes_f4in_diff_expressed.ByP_filtered$padj <= 0.1)), 
                           length(which(allRes_f4_up.MoF_filtered$padj <= 0.1)),length(which(allRes_f4_Down.MoF_filtered$padj <= 0.1)), length(which(allRes_f4in_diff_expressed.MoF_filtered$padj <= 0.1)),
                           length(which(allRes_f4_up.CeC_filtered$padj <= 0.1)),length(which(allRes_f4_Down.CeC_filtered$padj <= 0.1)), length(which(allRes_f4in_diff_expressed.CeC_filtered$padj <= 0.1)))


Go_terms_less_0.05_FDR <- c(length(which(allRes_meso_vs_idio_up.ByP_filtered$padj <= 0.05)),length(which(allRes_meso_vs_idio_Down.ByP_filtered$padj <= 0.05)), length(which(allRes_meso_vs_idio_diff_expressed.ByP_filtered$padj <= 0.05)), 
                           length(which(allRes_meso_vs_idio_up.MoF_filtered$padj <= 0.05)),length(which(allRes_meso_vs_idio_Down.MoF_filtered$padj <= 0.05)), length(which(allRes_meso_vs_idio_diff_expressed.MoF_filtered$padj <= 0.05)),
                           length(which(allRes_meso_vs_idio_up.CeC_filtered$padj <= 0.05)),length(which(allRes_meso_vs_idio_Down.CeC_filtered$padj <= 0.05)), length(which(allRes_meso_vs_idio_diff_expressed.CeC_filtered$padj <= 0.05)), 
                           length(which(allRes_f1_up.ByP_filtered$padj <= 0.1)),length(which(allRes_f1_Down.ByP_filtered$padj <= 0.1)), length(which(allRes_f1in_diff_expressed.ByP_filtered$padj <= 0.1)), 
                           length(which(allRes_f1_up.MoF_filtered$padj <= 0.1)),length(which(allRes_f1_Down.MoF_filtered$padj <= 0.1)), length(which(allRes_f1in_diff_expressed.MoF_filtered$padj <= 0.1)),
                           length(which(allRes_f1_up.CeC_filtered$padj <= 0.1)),length(which(allRes_f1_Down.CeC_filtered$padj <= 0.1)), length(which(allRes_f1in_diff_expressed.CeC_filtered$padj <= 0.1)), 
                           length(which(allRes_f4_up.ByP_filtered$padj <= 0.1)),length(which(allRes_f4_Down.ByP_filtered$padj <= 0.1)), length(which(allRes_f4in_diff_expressed.ByP_filtered$padj <= 0.1)), 
                           length(which(allRes_f4_up.MoF_filtered$padj <= 0.1)),length(which(allRes_f4_Down.MoF_filtered$padj <= 0.1)), length(which(allRes_f4in_diff_expressed.MoF_filtered$padj <= 0.1)),
                           length(which(allRes_f4_up.CeC_filtered$padj <= 0.1)),length(which(allRes_f4_Down.CeC_filtered$padj <= 0.1)), length(which(allRes_f4in_diff_expressed.CeC_filtered$padj <= 0.1)))

#Gene universe 
gene_universe <- c(resultFisher_idio_up.P@geneData[1], resultFisher_idio_Down.P@geneData[1], resultFisher_idio_diff_expressed.P@geneData[1], 
                    resultFisher_idio_up.F@geneData[1], resultFisher_idio_Down.F@geneData[1], resultFisher_idio_diff_expressed.F@geneData[1],
                    resultFisher_idio_up.C@geneData[1], resultFisher_idio_Down.C@geneData[1], resultFisher_idio_diff_expressed.C@geneData[1], 
                    resultFisher_f1_up.P@geneData[1], resultFisher_f1_Down.P@geneData[1], resultFisher_f1_diff_expressed.P@geneData[1], 
                    resultFisher_f1_up.F@geneData[1], resultFisher_f1_Down.F@geneData[1], resultFisher_f1_diff_expressed.M@geneData[1],
                    resultFisher_f1_up.C@geneData[1], resultFisher_f1_Down.C@geneData[1], resultFisher_f1_diff_expressed.C@geneData[1], 
                    resultFisher_f4_up.P@geneData[1], resultFisher_f4_Down.P@geneData[1], resultFisher_f4_diff_expressed.P@geneData[1], 
                    resultFisher_f4_up.F@geneData[1], resultFisher_f4_Down.F@geneData[1], resultFisher_f4_diff_expressed.M@geneData[1],
                    resultFisher_f4_up.C@geneData[1], resultFisher_f4_Down.C@geneData[1], resultFisher_f4_diff_expressed.C@geneData[1])

algorithms <- rep("weight01", 27)


##make the final table
summary_file <- data.frame("gene_set" = gene_set, "ontology" = Ontology,"Significant genes" = Significant_genes,
                           "significant genes with go terms" = feasible_genes, "Go terms scored" = GO_terms_score, 
                           "Go terms with p < 0.01" = Go_terms_less_0.01, "Go terms with p < 0.1 after FDR" = Go_terms_less_0.1_FDR, 
                           "Go terms with p < 0.05 after FDR" =  Go_terms_less_0.05_FDR, "Gene universe nº" = gene_universe, "algortihm" = algorithms)

write.xlsx(summary_file, file = "summary_file.xlsx", sheetName = "Folha 1", row.names = F)




###Make geneByTerms tables----

##FOR meso vs idio ----
##BP.up 
genes_UP_idio_vs_meso <- read.csv("DESeq2_tables/meso_vs_idio_sigUP.csv", header = TRUE)[,2]
UP.annotated <- meso_vs_idio_sub_annot[genes_UP_idio_vs_meso,]

#get list of genes with the GO_term 
allRes_meso_vs_idio_up.ByP_filtered$genes <- sapply(allRes_meso_vs_idio_up.ByP_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_meso_idio_up.P, x) 
  genes[[1]][genes[[1]] %in% row.names(UP.annotated)]
})


GoTerm_genes_up_file_name.B <- "IDs_byBP_UP.xlsx"
n_terms <- 50

for (i in 1:n_terms){
  genes_id <- allRes_meso_vs_idio_up.ByP_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_meso_vs_idio_up.ByP_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/meso_vs_idio/", GoTerm_genes_up_file_name.B), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/meso_vs_idio/", GoTerm_genes_up_file_name.B), sheetName = str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}


##BP.down 
genes_Down_idio_vs_meso <- read.csv("DESeq2_tables/meso_vs_idio_sigDOWN.csv", header = TRUE)[,2]
Down.annotated <- meso_vs_idio_sub_annot[genes_Down_idio_vs_meso,]


allRes_meso_vs_idio_Down.ByP_filtered$genes <- sapply(allRes_meso_vs_idio_Down.ByP_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_meso_idio_Down.P, x) 
  genes[[1]][genes[[1]] %in% row.names(Down.annotated)]
})

GoTerm_genes_down_file_name.B <- "IDs_byBP_Down.xlsx"
n_terms <- 50


for (i in 1:n_terms){
  genes_id <- allRes_meso_vs_idio_Down.ByP_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_meso_vs_idio_Down.ByP_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/meso_vs_idio/", GoTerm_genes_down_file_name.B), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/meso_vs_idio/", GoTerm_genes_down_file_name.B), sheetName = str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}


##BP.diff_expressed  
diff_idio.BP <- rbind(UP.annotated, Down.annotated)
genes_diff_expressed <- row.names(diff_idio.BP)


allRes_meso_vs_idio_diff_expressed.ByP_filtered$genes <- sapply(allRes_meso_vs_idio_diff_expressed.ByP_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_meso_idio_diff_expressed.P, x) 
  genes[[1]][genes[[1]] %in% row.names(diff_idio.BP)]
})

GoTerm_genes_diff_file_name.B <- "IDs_byBP_Diff_expressed.xlsx"
n_terms <- 50


#Fix some terms due to xlsx incompatibilities
allRes_meso_vs_idio_diff_expressed.ByP_filtered$Term[31] <- " photosynthetic electron transport chain"


for (i in 1:n_terms){
  genes_id <- allRes_meso_vs_idio_diff_expressed.ByP_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_meso_vs_idio_diff_expressed.ByP_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/meso_vs_idio/", GoTerm_genes_diff_file_name.B), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/meso_vs_idio/", GoTerm_genes_diff_file_name.B), sheetName = str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}



##MF.up 
genes_UP_idio_vs_meso <- read.csv("DESeq2_tables/meso_vs_idio_sigUP.csv", header = TRUE)[,2]
UP.annotated <- meso_vs_idio_sub_annot[genes_UP_idio_vs_meso,]

#get list of genes with the GO_term 
allRes_meso_vs_idio_up.MoF_filtered$genes <- sapply(allRes_meso_vs_idio_up.MoF_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_meso_idio_up.F, x) 
  genes[[1]][genes[[1]] %in% row.names(UP.annotated)]
})


GoTerm_genes_up_file_name.F <- "IDs_byMF_UP.xlsx"
n_terms <- 50

for (i in 1:n_terms){
  genes_id <- allRes_meso_vs_idio_up.MoF_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_meso_vs_idio_up.MoF_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/meso_vs_idio/", GoTerm_genes_up_file_name.F), sheetName = str_replace(Go_term, ":", " ") )
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/meso_vs_idio/", GoTerm_genes_up_file_name.F), sheetName = str_replace(Go_term, ":", " ") , append = T)
  }
  gc()
}


##MF.down 
genes_Down_idio_vs_meso <- read.csv("DESeq2_tables/meso_vs_idio_sigDOWN.csv", header = TRUE)[,2]
Down.annotated <- meso_vs_idio_sub_annot[genes_Down_idio_vs_meso,]


allRes_meso_vs_idio_Down.MoF_filtered$genes <- sapply(allRes_meso_vs_idio_Down.MoF_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_meso_idio_Down.F, x) 
  genes[[1]][genes[[1]] %in% row.names(Down.annotated)]
})

GoTerm_genes_down_file_name.B <- "IDs_byMF_Down.xlsx"
n_terms <- 50

for (i in 1:n_terms){
  genes_id <- allRes_meso_vs_idio_Down.MoF_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_meso_vs_idio_Down.MoF_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/meso_vs_idio/", GoTerm_genes_down_file_name.B), sheetName = str_replace(Go_term, ":", " ") )
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/meso_vs_idio/", GoTerm_genes_down_file_name.B), sheetName = str_replace(Go_term, ":", " ") , append = T)
  }
  gc()
}


##MF.diff_expressed  
diff_idio.MF <- rbind(UP.annotated, Down.annotated)
genes_diff_expressed <- row.names(diff_idio.MF)


allRes_meso_vs_idio_diff_expressed.MoF_filtered$genes <- sapply(allRes_meso_vs_idio_diff_expressed.MoF_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_meso_idio_diff_expressed.F, x) 
  genes[[1]][genes[[1]] %in% row.names(diff_idio.MF)]
})

GoTerm_genes_diff_file_name.B <- "IDs_byMF_Diff_expressed.xlsx"
n_terms <- 50


for (i in 1:n_terms){
  genes_id <- allRes_meso_vs_idio_diff_expressed.MoF_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_meso_vs_idio_diff_expressed.MoF_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/meso_vs_idio/", GoTerm_genes_diff_file_name.B), sheetName = str_replace(Go_term, ":", " ") )
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/meso_vs_idio/", GoTerm_genes_diff_file_name.B), sheetName =str_replace(Go_term, ":", " ") , append = T)
  }
  gc()
}



##CC.up 
genes_UP_idio_vs_meso <- read.csv("DESeq2_tables/meso_vs_idio_sigUP.csv", header = TRUE)[,2]
UP.annotated <- meso_vs_idio_sub_annot[genes_UP_idio_vs_meso,]

#get list of genes with the GO_term 
allRes_meso_vs_idio_up.CeC_filtered$genes <- sapply(allRes_meso_vs_idio_up.CeC_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_meso_idio_up.C, x) 
  genes[[1]][genes[[1]] %in% row.names(UP.annotated)]
})


GoTerm_genes_up_file_name.C <- "IDs_byCC_UP.xlsx"
n_terms <- 50

#Fix some GOTerms names that are invalid for xlsx

for (i in 1:n_terms){
  genes_id <- allRes_meso_vs_idio_up.CeC_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_meso_vs_idio_up.CeC_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/meso_vs_idio/", GoTerm_genes_up_file_name.C), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/meso_vs_idio/", GoTerm_genes_up_file_name.C), str_replace(Go_term, ":", " ") , append = T)
  }
  gc()
}


##CC.down 
genes_Down_idio_vs_meso <- read.csv("DESeq2_tables/meso_vs_idio_sigDOWN.csv", header = TRUE)[,2]
Down.annotated <- meso_vs_idio_sub_annot[genes_Down_idio_vs_meso,]


allRes_meso_vs_idio_Down.CeC_filtered$genes <- sapply(allRes_meso_vs_idio_Down.CeC_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_meso_idio_Down.C, x) 
  genes[[1]][genes[[1]] %in% row.names(Down.annotated)]
})

GoTerm_genes_down_file_name.B <- "IDs_byCC_Down.xlsx"
n_terms <- 50


#Fix some terms due to xlsx incompatibilities


for (i in 1:n_terms){
  genes_id <- allRes_meso_vs_idio_Down.CeC_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_meso_vs_idio_Down.CeC_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/meso_vs_idio/", GoTerm_genes_down_file_name.B), sheetName = str_replace(Go_term, ":", " ") )
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/meso_vs_idio/", GoTerm_genes_down_file_name.B), sheetName = str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}


##CC.diff_expressed  
diff_idio.CC <- rbind(UP.annotated, Down.annotated)
genes_diff_expressed <- row.names(diff_idio.CC)


allRes_meso_vs_idio_diff_expressed.CeC_filtered$genes <- sapply(allRes_meso_vs_idio_diff_expressed.CeC_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_meso_idio_diff_expressed.C, x) 
  genes[[1]][genes[[1]] %in% row.names(diff_idio.CC)]
})

GoTerm_genes_diff_file_name.B <- "IDs_byCC_Diff_expressed.xlsx"
n_terms <- 50



for (i in 1:n_terms){
  genes_id <- allRes_meso_vs_idio_diff_expressed.CeC_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_meso_vs_idio_diff_expressed.CeC_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/meso_vs_idio/", GoTerm_genes_diff_file_name.B), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/meso_vs_idio/", GoTerm_genes_diff_file_name.B), sheetName =str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}





##For F1in vs F1out ----
##BP.up 
genes_UP_f1 <- read.csv("DESeq2_tables/F1in_vs_F1out_sigUP.csv", header = TRUE)[,2]
UP.annotated <- f1_sub_annot[genes_UP_f1,]

#get list of genes with the GO_term 
allRes_f1_up.ByP_filtered$genes <- sapply(allRes_f1_up.ByP_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_f1_up.P, x) 
  genes[[1]][genes[[1]] %in% row.names(UP.annotated)]
})


GoTerm_genes_up_file_name.B <- "IDs_byBP_UP.xlsx"
n_terms <- 50

for (i in 1:n_terms){
  genes_id <- allRes_f1_up.ByP_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_f1_up.ByP_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/F1in_vs_f1out/", GoTerm_genes_up_file_name.B), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/F1in_vs_f1out/", GoTerm_genes_up_file_name.B), sheetName = str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}


##BP.down 
genes_Down_idio_vs_meso <- read.csv("DESeq2_tables/F1in_vs_F1out_sigDOWN.csv", header = TRUE)[,2]
Down.annotated <- f1_sub_annot[genes_Down_idio_vs_meso,]


allRes_f1_Down.ByP_filtered$genes <- sapply(allRes_f1_Down.ByP_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_f1_Down.P, x) 
  genes[[1]][genes[[1]] %in% row.names(Down.annotated)]
})

GoTerm_genes_down_file_name.B <- "IDs_byBP_Down.xlsx"
n_terms <- 50


for (i in 1:n_terms){
  genes_id <- allRes_f1_Down.ByP_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_f1_Down.ByP_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/F1in_vs_f1out/", GoTerm_genes_down_file_name.B), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/F1in_vs_f1out/", GoTerm_genes_down_file_name.B), sheetName = str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}


##BP.diff_expressed 
diff_f1.BP <- rbind(UP.annotated, Down.annotated)
genes_diff_expressed <- row.names(diff_f1.BP)


allRes_f1in_diff_expressed.ByP_filtered$genes <- sapply(allRes_f1in_diff_expressed.ByP_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_f1_diff_expressed.P, x) 
  genes[[1]][genes[[1]] %in% row.names(diff_f1.BP)]
})

GoTerm_genes_diff_file_name.B <- "IDs_byBP_Diff_expressed.xlsx"
n_terms <- 50



for (i in 1:n_terms){
  genes_id <- allRes_f1in_diff_expressed.ByP_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_f1in_diff_expressed.ByP_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/F1in_vs_f1out/", GoTerm_genes_diff_file_name.B), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/F1in_vs_f1out/", GoTerm_genes_diff_file_name.B), sheetName = str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}



##MF.up 
genes_UP_idio_vs_meso <- read.csv("DESeq2_tables/F1in_vs_F1out_sigUP.csv", header = TRUE)[,2]
UP.annotated <- f1_sub_annot[genes_UP_idio_vs_meso,]

#get list of genes with the GO_term 
allRes_f1_up.MoF_filtered$genes <- sapply(allRes_f1_up.MoF_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_f1_up.F, x) 
  genes[[1]][genes[[1]] %in% row.names(UP.annotated)]
})


GoTerm_genes_up_file_name.F <- "IDs_byMF_UP.xlsx"
n_terms <- 50

for (i in 1:n_terms){
  genes_id <- allRes_f1_up.MoF_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_f1_up.MoF_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/F1in_vs_f1out/", GoTerm_genes_up_file_name.F), sheetName = str_replace(Go_term, ":", " ") )
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/F1in_vs_f1out/", GoTerm_genes_up_file_name.F), sheetName = str_replace(Go_term, ":", " ") , append = T)
  }
  gc()
}


##MF.down 
genes_Down_idio_vs_meso <- read.csv("DESeq2_tables/F1in_vs_F1out_sigDOWN.csv", header = TRUE)[,2]
Down.annotated <- f1_sub_annot[genes_Down_idio_vs_meso,]


allRes_f1_Down.MoF_filtered$genes <- sapply(allRes_f1_Down.MoF_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_f1_Down.F, x) 
  genes[[1]][genes[[1]] %in% row.names(Down.annotated)]
})

GoTerm_genes_down_file_name.B <- "IDs_byMF_Down.xlsx"
n_terms <- 50

for (i in 1:n_terms){
  genes_id <- allRes_f1_Down.MoF_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_f1_Down.MoF_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/F1in_vs_f1out/", GoTerm_genes_down_file_name.B), sheetName = str_replace(Go_term, ":", " ") )
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/F1in_vs_f1out/", GoTerm_genes_down_file_name.B), sheetName = str_replace(Go_term, ":", " ") , append = T)
  }
  gc()
}


##MF.diff_expressed  
diff_f1.MF <- rbind(UP.annotated, Down.annotated)
genes_diff_expressed <- row.names(diff_f1.MF)


allRes_f1in_diff_expressed.MoF_filtered$genes <- sapply(allRes_f1in_diff_expressed.MoF_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_f1_diff_expressed.M, x) 
  genes[[1]][genes[[1]] %in% row.names(diff_f1.MF)]
})

GoTerm_genes_diff_file_name.B <- "IDs_byMF_Diff_expressed.xlsx"
n_terms <- 50


for (i in 1:n_terms){
  genes_id <- allRes_f1in_diff_expressed.MoF_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_f1in_diff_expressed.MoF_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/F1in_vs_f1out/", GoTerm_genes_diff_file_name.B), sheetName = str_replace(Go_term, ":", " ") )
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/F1in_vs_f1out/", GoTerm_genes_diff_file_name.B), sheetName =str_replace(Go_term, ":", " ") , append = T)
  }
  gc()
}



##CC.up 
genes_UP_idio_vs_meso <- read.csv("DESeq2_tables/F1in_vs_F1out_sigUP.csv", header = TRUE)[,2]
UP.annotated <- f1_sub_annot[genes_UP_idio_vs_meso,]

#get list of genes with the GO_term 
allRes_f1_up.CeC_filtered$genes <- sapply(allRes_f1_up.CeC_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_f1_up.C, x) 
  genes[[1]][genes[[1]] %in% row.names(UP.annotated)]
})


GoTerm_genes_up_file_name.C <- "IDs_byCC_UP.xlsx"
n_terms <- 50

#Fix some GOTerms names that are invalid for xlsx

for (i in 1:n_terms){
  genes_id <- allRes_f1_up.CeC_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_f1_up.CeC_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/F1in_vs_f1out/", GoTerm_genes_up_file_name.C), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/F1in_vs_f1out/", GoTerm_genes_up_file_name.C), str_replace(Go_term, ":", " ") , append = T)
  }
  gc()
}


##CC.down 
genes_Down_idio_vs_meso <- read.csv("DESeq2_tables/F1in_vs_F1out_sigDOWN.csv", header = TRUE)[,2]
Down.annotated <- f1_sub_annot[genes_Down_idio_vs_meso,]


allRes_f1_Down.CeC_filtered$genes <- sapply(allRes_f1_Down.CeC_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_f1_Down.C, x) 
  genes[[1]][genes[[1]] %in% row.names(Down.annotated)]
})

GoTerm_genes_down_file_name.B <- "IDs_byCC_Down.xlsx"
n_terms <- 50


#Fix some terms due to xlsx incompatibilities


for (i in 1:n_terms){
  genes_id <- allRes_f1_Down.CeC_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_f1_Down.CeC_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/F1in_vs_f1out/", GoTerm_genes_down_file_name.B), sheetName = str_replace(Go_term, ":", " ") )
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/F1in_vs_f1out/", GoTerm_genes_down_file_name.B), sheetName = str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}


##CC.diff_expressed  
diff_f1.CC <- rbind(UP.annotated, Down.annotated)
genes_diff_expressed <- row.names(diff_f1.CC)


allRes_f1in_diff_expressed.CeC_filtered$genes <- sapply(allRes_f1in_diff_expressed.CeC_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_f1_diff_expressed.C, x) 
  genes[[1]][genes[[1]] %in% row.names(diff_f1.CC)]
})

GoTerm_genes_diff_file_name.B <- "IDs_byCC_Diff_expressed.xlsx"
n_terms <- 50



for (i in 1:n_terms){
  genes_id <- allRes_f1in_diff_expressed.CeC_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_f1in_diff_expressed.CeC_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/F1in_vs_f1out/", GoTerm_genes_diff_file_name.B), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/F1in_vs_f1out/", GoTerm_genes_diff_file_name.B), sheetName =str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}






##For f4in vs f4out ----
##BP.up 
genes_UP_f4 <- read.csv("DESeq2_tables/f4in_vs_f4out_sigUP.csv", header = TRUE)[,2]
UP.annotated <- f4_sub_annot[genes_UP_f4,]

#get list of genes with the GO_term 
allRes_f4_up.ByP_filtered$genes <- sapply(allRes_f4_up.ByP_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_f4_up.P, x) 
  genes[[1]][genes[[1]] %in% row.names(UP.annotated)]
})


GoTerm_genes_up_file_name.B <- "IDs_byBP_UP.xlsx"
n_terms <- 50

for (i in 1:n_terms){
  genes_id <- allRes_f4_up.ByP_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_f4_up.ByP_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/f4in_vs_f4out/", GoTerm_genes_up_file_name.B), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/f4in_vs_f4out/", GoTerm_genes_up_file_name.B), sheetName = str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}


##BP.down 
genes_Down_idio_vs_meso <- read.csv("DESeq2_tables/f4in_vs_f4out_sigDOWN.csv", header = TRUE)[,2]
Down.annotated <- f4_sub_annot[genes_Down_idio_vs_meso,]


allRes_f4_Down.ByP_filtered$genes <- sapply(allRes_f4_Down.ByP_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_f4_Down.P, x) 
  genes[[1]][genes[[1]] %in% row.names(Down.annotated)]
})

GoTerm_genes_down_file_name.B <- "IDs_byBP_Down.xlsx"
n_terms <- 50


for (i in 1:n_terms){
  genes_id <- allRes_f4_Down.ByP_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_f4_Down.ByP_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/f4in_vs_f4out/", GoTerm_genes_down_file_name.B), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/f4in_vs_f4out/", GoTerm_genes_down_file_name.B), sheetName = str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}


##BP.diff_expressed 
diff_f4.BP <- rbind(UP.annotated, Down.annotated)
genes_diff_expressed <- row.names(diff_f4.BP)


allRes_f4in_diff_expressed.ByP_filtered$genes <- sapply(allRes_f4in_diff_expressed.ByP_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_f4_diff_expressed.P, x) 
  genes[[1]][genes[[1]] %in% row.names(diff_f4.BP)]
})

GoTerm_genes_diff_file_name.B <- "IDs_byBP_Diff_expressed.xlsx"
n_terms <- 50



for (i in 1:n_terms){
  genes_id <- allRes_f4in_diff_expressed.ByP_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_f4in_diff_expressed.ByP_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/f4in_vs_f4out/", GoTerm_genes_diff_file_name.B), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/f4in_vs_f4out/", GoTerm_genes_diff_file_name.B), sheetName = str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}



##MF.up 
genes_UP_idio_vs_meso <- read.csv("DESeq2_tables/f4in_vs_f4out_sigUP.csv", header = TRUE)[,2]
UP.annotated <- f4_sub_annot[genes_UP_idio_vs_meso,]

#get list of genes with the GO_term 
allRes_f4_up.MoF_filtered$genes <- sapply(allRes_f4_up.MoF_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_f4_up.F, x) 
  genes[[1]][genes[[1]] %in% row.names(UP.annotated)]
})


GoTerm_genes_up_file_name.F <- "IDs_byMF_UP.xlsx"
n_terms <- 50

for (i in 1:n_terms){
  genes_id <- allRes_f4_up.MoF_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_f4_up.MoF_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/f4in_vs_f4out/", GoTerm_genes_up_file_name.F), sheetName = str_replace(Go_term, ":", " ") )
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/f4in_vs_f4out/", GoTerm_genes_up_file_name.F), sheetName = str_replace(Go_term, ":", " ") , append = T)
  }
  gc()
}


##MF.down 
genes_Down_idio_vs_meso <- read.csv("DESeq2_tables/f4in_vs_f4out_sigDOWN.csv", header = TRUE)[,2]
Down.annotated <- f4_sub_annot[genes_Down_idio_vs_meso,]


allRes_f4_Down.MoF_filtered$genes <- sapply(allRes_f4_Down.MoF_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_f4_Down.F, x) 
  genes[[1]][genes[[1]] %in% row.names(Down.annotated)]
})

GoTerm_genes_down_file_name.B <- "IDs_byMF_Down.xlsx"
n_terms <- 50

for (i in 1:n_terms){
  genes_id <- allRes_f4_Down.MoF_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_f4_Down.MoF_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/f4in_vs_f4out/", GoTerm_genes_down_file_name.B), sheetName = str_replace(Go_term, ":", " ") )
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/f4in_vs_f4out/", GoTerm_genes_down_file_name.B), sheetName = str_replace(Go_term, ":", " ") , append = T)
  }
  gc()
}


##MF.diff_expressed  
diff_f4.MF <- rbind(UP.annotated, Down.annotated)
genes_diff_expressed <- row.names(diff_f4.MF)


allRes_f4in_diff_expressed.MoF_filtered$genes <- sapply(allRes_f4in_diff_expressed.MoF_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_f4_diff_expressed.M, x) 
  genes[[1]][genes[[1]] %in% row.names(diff_f4.MF)]
})

GoTerm_genes_diff_file_name.B <- "IDs_byMF_Diff_expressed.xlsx"
n_terms <- 50


for (i in 1:n_terms){
  genes_id <- allRes_f4in_diff_expressed.MoF_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_f4in_diff_expressed.MoF_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/f4in_vs_f4out/", GoTerm_genes_diff_file_name.B), sheetName = str_replace(Go_term, ":", " ") )
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/f4in_vs_f4out/", GoTerm_genes_diff_file_name.B), sheetName =str_replace(Go_term, ":", " ") , append = T)
  }
  gc()
}



##CC.up 
genes_UP_idio_vs_meso <- read.csv("DESeq2_tables/f4in_vs_f4out_sigUP.csv", header = TRUE)[,2]
UP.annotated <- f4_sub_annot[genes_UP_idio_vs_meso,]

#get list of genes with the GO_term 
allRes_f4_up.CeC_filtered$genes <- sapply(allRes_f4_up.CeC_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_f4_up.C, x) 
  genes[[1]][genes[[1]] %in% row.names(UP.annotated)]
})


GoTerm_genes_up_file_name.C <- "IDs_byCC_UP.xlsx"
n_terms <- 50

#Fix some GOTerms names that are invalid for xlsx

for (i in 1:n_terms){
  genes_id <- allRes_f4_up.CeC_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_f4_up.CeC_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/f4in_vs_f4out/", GoTerm_genes_up_file_name.C), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/f4in_vs_f4out/", GoTerm_genes_up_file_name.C), str_replace(Go_term, ":", " ") , append = T)
  }
  gc()
}


##CC.down 
genes_Down_idio_vs_meso <- read.csv("DESeq2_tables/f4in_vs_f4out_sigDOWN.csv", header = TRUE)[,2]
Down.annotated <- f4_sub_annot[genes_Down_idio_vs_meso,]


allRes_f4_Down.CeC_filtered$genes <- sapply(allRes_f4_Down.CeC_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_f4_Down.C, x) 
  genes[[1]][genes[[1]] %in% row.names(Down.annotated)]
})

GoTerm_genes_down_file_name.B <- "IDs_byCC_Down.xlsx"
n_terms <- 50


#Fix some terms due to xlsx incompatibilities


for (i in 1:n_terms){
  genes_id <- allRes_f4_Down.CeC_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_f4_Down.CeC_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/f4in_vs_f4out/", GoTerm_genes_down_file_name.B), sheetName = str_replace(Go_term, ":", " ") )
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/f4in_vs_f4out/", GoTerm_genes_down_file_name.B), sheetName = str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}


##CC.diff_expressed  
diff_f4.CC <- rbind(UP.annotated, Down.annotated)
genes_diff_expressed <- row.names(diff_f4.CC)


allRes_f4in_diff_expressed.CeC_filtered$genes <- sapply(allRes_f4in_diff_expressed.CeC_filtered$GO.ID, function(x) { # applies this function to all GOterms in the results list
  genes <- genesInTerm(GoData_f4_diff_expressed.C, x) 
  genes[[1]][genes[[1]] %in% row.names(diff_f4.CC)]
})

GoTerm_genes_diff_file_name.B <- "IDs_byCC_Diff_expressed.xlsx"
n_terms <- 50



for (i in 1:n_terms){
  genes_id <- allRes_f4in_diff_expressed.CeC_filtered$genes[i] #get the genes from each term 
  Go_term <- names(genes_id) #Get the Go term code 
  Go_name <- allRes_f4in_diff_expressed.CeC_filtered$Term[i] #get the Go term name 
  genes_id <- as.character(unlist(genes_id)) #transform into a list to subset
  GO_term_table <- annotation_file[genes_id, ]
  GO_term_table$GOTERM <- Go_term
  GO_term_table$GOName <- Go_name
  GO_term_table <- GO_term_table[,c(17,18,1:16)] #reorder the collum 
  if (i == 1){
    #in the first iteration create the excel file 
    write.xlsx(GO_term_table, file = paste("results/f4in_vs_f4out/", GoTerm_genes_diff_file_name.B), sheetName = str_replace(Go_term, ":", " "))
  }
  else{
    #in interation >= 2 append a new table to the excel file 
    write.xlsx(GO_term_table, file = paste("results/f4in_vs_f4out/", GoTerm_genes_diff_file_name.B), sheetName =str_replace(Go_term, ":", " "), append = T)
  }
  gc()
}