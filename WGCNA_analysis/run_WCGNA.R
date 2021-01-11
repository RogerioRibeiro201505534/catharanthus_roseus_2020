####Use this script to build a network with WCGNA -> Using all samples and all genes 

#Misc options
setwd("C:/Users/Faculdade/Desktop/Dissertação/results 2020/mapping_to_genome/11_WCGNA")
rm(list = ls())

#Load libraries 
library(WGCNA)
library(DESeq2)
library("tximport")
library(xlsx)
library(ComplexHeatmap)
#library(sva)

#Make output dir 
out.dir <- "results"
if (!dir.exists(out.dir)){
  dir.create(out.dir)
}


#Load data (totalptt and leaves were discarded due to not being used in comparisons)
sampleNames <- c('F1P1_Int', 'F1P2_Int', 'F1P3_Int','F1P1_Ext', 'F1P2_Ext', 'F1P3_Ext', 'F4P1_Int', 'F4P2_Int', 'F4P3_Int','F4P1_Ext', 'F4P2_Ext', 'F4P3_Ext', 'idio_1', 'idio_2', 'idio_3', 'meso_1', 'meso_2', 'meso_3')
sampleCondition <- c("f1in_f1in","f1in_f1in","f1in_f1in","f1out_f1out","f1out_f1out","f1out_f1out", "f4in_f4in","f4in_f4in","f4in_f4in","f4out_f4out","f4out_f4out","f4out_f4out", "idio", "idio", "idio", "meso", "meso", "meso")
sampleFiles <- "abundance.tsv"

tx2gene <- read.csv(file = "transcript2gene.tsv", sep = "\t", header = F)


#### Load expression data and normalize ####
# create dataframe that becomes a deseq table.
directory <- "C:/Users/Faculdade/Desktop/Dissertação/results 2020/mapping_to_genome/10_WCGNA/counts/"
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
sampleTable$condition <- as.factor(sampleTable$condition)
files <-  paste("counts", sampleTable$sampleName, sampleTable$fileName, sep = "/")
names(files) <- paste0(sampleTable$sampleName)

kallisto_Counts <- tximport(files, type = "kallisto", tx2gene = tx2gene, dropInfReps = TRUE)

ddkallisto <- DESeqDataSetFromTximport(kallisto_Counts, colData = sampleTable, design = ~ condition)


#expression_data_frame <- assay(ddkallisto)

#normalize for batch effects:
#batch <- c(rep(2,12), rep(1,6))

#adjustedCounts <- ComBat_seq(expression_data_frame, batch = batch)



#create a DESeq object (used to normalize date and to apply the rlog transform data)
#coldata <- data.frame("condition" = sampleTable$condition, "type" = rep("PE", 18))
#row.names(coldata) <- sampleTable$sampleName
#ddsmatrix <- DESeqDataSetFromMatrix(countData = adjustedCounts, colData = coldata, 
                                    #design = ~ condition)

dds <- DESeq(ddkallisto) #37320


#Filtering ---- 
#Filter genes that have a low expression value. Keep genes with > 20 counts in at least 
#50% of the samples(9). This is donne for a) keep a manageble dataset and b) remove low expressed 
#Which also removes in some way low varience genes, which are most likely to be noise

expression_unfiltered <- data.frame(assay(dds))
filtering_vector <- c()

for (i in 1:length(expression_unfiltered$F1P1_Int)){
  if (sum(expression_unfiltered[i,] > 20) > 9){
    filtering_vector[i] <- T
  }
  else{
    filtering_vector[i] <- F
  }
}

dds <- dds[filtering_vector, ] 
length(row.names(dds)) #15203




#apply the rlog transformation 
rld <- rlog(dds, blind=TRUE) # betaPriorVar:  if missing is estimated from the data

#perform QC tests on the samples---- 
#get the dataframe for analysis 
expression_matrix <- t(as.data.frame(assay(rld)))

#check genes and samples with too many missing values 
gsg = goodSamplesGenes(expression_matrix, verbose = 3)

#cluster samples too see if there are any obvious outlier 




sampleTree = hclust(dist(expression_matrix), method = "average")

sampleTree$labels <- c("F1P1In", "F1P2In", "F1P3In",
                       "F1P1Out", "F1P2Out", "F1P3Out",
                       "F4P1In", "F4P2In", "F4P3In",
                       "F4P1Out", "F4P2Out", "F4P3Out",
                       "idio 1", "idio 2", "idio 3", 
                       "meso 1", "meso 2", "meso 3")

par(cex = 0.6)
par(mar = c(0,4,2,0))
png(paste(out.dir, "sample_clustering_plot.png", sep = "/"), res = 300, height = 1500, width = 2500)
plot(sampleTree, main = "", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

#Build network (step-wise)----

#pick the more meaningfull Beta variable value -> test values from 1:30
powers = seq(10,40, by = 2)
sft <- pickSoftThreshold(expression_matrix, powerVector = powers, verbose = 5, networkType = "signed")

#plot mean conncectivity and R^2 to the free
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

png(paste(out.dir, "scale_independence.png", sep = "/"), res = 300, height = 1200, width = 1500)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
dev.off()

# Mean connectivity as a function of the soft-thresholding power
png(paste(out.dir, "Mean_connectivity.png", sep = "/"), height = 1200, width = 1500)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


#calculate adjancency and Topological Overlap Matrix (TOM)
softPower = 30
adjacency = adjacency(expression_matrix, power = softPower, type = "signed")

TOM = TOMsimilarity(adjacency,TOMType = "signed")
dissTOM = 1-TOM

row.names(dissTOM) <- row.names(adjacency)
colnames(dissTOM) <- colnames(adjacency)


#Clustering using TOM 
geneTree = hclust(as.dist(dissTOM), method = "average")
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)



#cluster the gene into modules of correlated genes----

#cluster the dendogram into modules 
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)

#plot the dendogram with the colour association 
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
png(paste(out.dir, "Gene_dendogram_with_colours.png", sep = "/"), res = 300, height = 2000, width = 3000)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


##Calculate each module eigengene to cluster closely related modules ---- 
# Calculate eigengenes
MEList = moduleEigengenes(expression_matrix, colors = dynamicColors, excludeGrey = T)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

#plot Module eigengenes hieraquical clustering 
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")


#Cluster similar modules 
# Plot the cut line into the dendrogram
MEDissThres = 0.18
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(expression_matrix, dynamicColors, 
                          cutHeight = MEDissThres, getNewUnassdME = FALSE, 
                          verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs


##save old dendogram of module eigengenes
png(paste(out.dir, "eigengeneclusterBeforeMerge.png", sep = "/"), res = 300, width = 3000, height = 2000)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()

#save dendogram of module eigengenes

# Calculate dissimilarity of merged module eigengenes
ME_merged.diss = 1-cor(mergedMEs)
# Cluster module eigengenes
METree_merged = hclust(as.dist(ME_merged.diss), method = "average")

#plot Module eigengenes hieraquical clustering 
png(paste(out.dir, "eigengeneclusterAfterMerge.png", sep = "/"), res = 300, width = 2000, height = 2000)
plot(METree_merged, main = "Clustering of module eigengenes after merge",
     xlab = "", sub = "")
dev.off()


png(paste(out.dir, "Gene_dendogram_with_colours.png", sep = "/"), width = 700)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

table(mergedColors)

png(paste(out.dir, "eigengene network.png", sep = "/"), res = 300, height = 1200, width = 1500) 
plotEigengeneNetworks(mergedMEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90, signed = TRUE, greyLabel = 0)

dev.off()

#save.image(file='WGCNA.RData')

#Correlate these modules with external traits (alkaloids levels)----

#Load metabolomic data
traitData_alkaloids <- read.csv(file = "alkaloids_rogerio.tsv", sep ="\t")

#test for Vindoline, catharanthina, AVLB (also for vinblastina)
traitData_alkaloids <- traitData_alkaloids[,c(1,4,6,8,10)]
traitData_alkaloids$Vinblastine[which(traitData_alkaloids$Vinblastine == 0)] <- NA

#Match row.names
traitData_alkaloids$Sample <- c("F1P1_Ext", "F1P1_Int", "F4P1_Ext", "F4P1_Int", 
                                "F1P2_Ext", "F1P2_Int", "F4P2_Ext", "F4P2_Int",
                                "F1P3_Ext", "F1P3_Int", "F4P3_Ext", "F4P3_Int",
                                "idio_1", "idio_2", "idio_3", 
                                "meso_1", "meso_2", "meso_3")

my_sample <- rownames(expression_matrix)
trait_Rows <- match(my_sample, traitData_alkaloids$Sample)
datTraits <- traitData_alkaloids[trait_Rows, -1]
row.names(datTraits) <- traitData_alkaloids[trait_Rows, 1]

##correlate the alkaloids levels with expression values and eigengenes
#mergedMEs are the ME of the merged modules == MEs in the tutorial 3

#First step calculates correlation, second the p.value
moduleTraitCor <- cor(mergedMEs, datTraits, use = "pairwise.complete.obs", method = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(expression_matrix))
moduleTraitPvalue_list <- moduleTraitPvalue[1:60]


#Adjust for pvalue
moduleTraitPvadj_list <- p.adjust(moduleTraitPvalue_list, method="fdr")

moduleTraitPadj <- moduleTraitPvalue
moduleTraitPadj[1:60] <- moduleTraitPvadj_list[1:60]

row.names(moduleTraitCor) <- paste("Module", seq(1,15))

#Plot the module trait correlation in a heatmap

png(paste(out.dir, "module_traitCorrelationPearson.png", sep = "/"), res = 300, width = 2800, height = 2000)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPadj, 1), ")", sep = "")



dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = row.names(moduleTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()

#and again without p value and just with *, with corelation rounded up
#p value < 0.0.5

png(paste(out.dir, "module_traitCorrelationPearson_1.png", sep = "/"), res = 300, width = 2800, height = 2000)
textMatrix = paste(signif(moduleTraitCor, 1),
                   ifelse(moduleTraitPadj <= 0.05,"*", " "), sep = "")



dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = row.names(moduleTraitCor),
               ySymbols = names(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()



png(paste(out.dir, "module_traitCorrelationPearson_2.png", sep = "/"), res = 300, width = 2800, height = 2000)
textMatrix = paste(ifelse(moduleTraitPadj <= 0.05,"*", " "), sep = "")



dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = row.names(moduleTraitCor),
               ySymbols = names(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()



#Calculate Gene significance and module Membership ----- 
#gene significance is defined as a the cor between a gene and a external trait
#module Membership is defined as the cor between a gene and the module eigengenes


modNames <-substring(names(mergedMEs), 3)

geneModuleMembership <- as.data.frame(cor(expression_matrix, mergedMEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(expression_matrix)))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#correct pvalue for multiple testing 
for (i in 1:length(MMPvalue)){
  if (i == 1){
    #initialize the vector
    MMPvalue_vector <- MMPvalue[,i]
  }
  else{
    MMPvalue_vector <- c(MMPvalue_vector, MMPvalue[,i])
  }
}

MMPadj_vector <- p.adjust(MMPvalue_vector, method="fdr")

MMadj <- MMPvalue

for (i in 0:(length(MMadj)-1)){
  MMadj[,i+1] <- MMPadj_vector[(i*length(MMadj[,1]) + 1):(length(MMadj[,1]) * (i+1))]
}



vindoline <- as.data.frame(datTraits$Vindoline)
names(vindoline) <- "vindoline"
catharanthine <- as.data.frame(datTraits$Catharanthine)
names(catharanthine) <- "Catharanthine"
vinblastine <- as.data.frame(datTraits$Vinblastine)
names(vinblastine) <- "vinblastine"
AVLB <- as.data.frame(datTraits$AVBL)
names(AVLB) <- "AVLB"

geneTraitSignificance <- data.frame("Vindoline" = cor(expression_matrix, vindoline, use = "p"), 
                                    "catharanthine" = cor(expression_matrix, catharanthine, use = "p"),
                                    "Vinblastine" = cor(expression_matrix, vinblastine, use = "p"), 
                                    "AVLB" = cor(expression_matrix, AVLB, use = "p"))
        

GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(expression_matrix)))

names(geneTraitSignificance) = paste("GS.", names(geneTraitSignificance), sep="")
names(GSPvalue) = paste("p.GS.", names(GSPvalue), sep="")


#correct p.value for multiple testing
for (i in 1:length(GSPvalue)){
  if (i == 1){
    #initialize the vector
    GSPvalue_vector <- GSPvalue[,i]
  }
  else{
    GSPvalue_vector <- c(GSPvalue_vector, GSPvalue[,i])
  }
}

GSPadj_vector <- p.adjust(GSPvalue_vector, method="fdr")
GSadj <- GSPvalue


for (i in 0:(length(GSadj)-1)){
  GSadj[,i+1] <- GSPadj_vector[(i*length(GSadj[,1]) + 1):(length(GSadj[,1]) * (i+1))]
}


#Get an annotation table with the module colour 

#load annotation from dammit with manua curation of the pathway genes 
dammit_annotation.final <- read.delim(file = "annotation/dammit_annotation_final.tsv")

#Make a table with each gene colour, ordered by colour used in downstream analsyis ---- 
#add the MM for the module it clustered and GS for all traits (+ adjusted p.values)
#could add MM for other modules later if need

modules_merged <- names(table(mergedColors))
flag = T
for (i in modules_merged){
  temp_genes <- colnames(expression_matrix)[mergedColors == i]
  temp_gene_table <- dammit_annotation.final[as.character(temp_genes), ]
  temp_gene_table$module <- i
  if (i != "grey"){
    temp_gene_table$MM <- geneModuleMembership[row.names(temp_gene_table),paste("MM", i, sep = "")]
    temp_gene_table$MMPadj <- MMadj[row.names(temp_gene_table),paste("p.MM", i, sep = "")]
    temp_gene_table$GSvindoline <-  geneTraitSignificance[row.names(temp_gene_table),1]
    temp_gene_table$GSCatharanthine <-  geneTraitSignificance[row.names(temp_gene_table),2]
    temp_gene_table$GSvinblastine <-  geneTraitSignificance[row.names(temp_gene_table),3]
    temp_gene_table$GSAVBL <-  geneTraitSignificance[row.names(temp_gene_table),4]
    temp_gene_table$GSVindolinePadj<-  GSadj[row.names(temp_gene_table),1]
    temp_gene_table$GSCatharantinePadj <-  GSadj[row.names(temp_gene_table),2]
    temp_gene_table$GSvinblastinePadj <-  GSadj[row.names(temp_gene_table),3]
    temp_gene_table$GSAVBLPadj <-  GSadj[row.names(temp_gene_table),4]
  }
  else{
    temp_gene_table$MM <- NA
    temp_gene_table$MMPadj <- NA
    temp_gene_table$GSvindoline <- NA
    temp_gene_table$GSCatharanthine <- NA
    temp_gene_table$GSvinblastine <- NA
    temp_gene_table$GSAVBL <-  NA
    temp_gene_table$GSVindolinePadj<-  NA
    temp_gene_table$GSCatharantinePadj <-  NA
    temp_gene_table$GSvinblastinePadj <-  NA
    temp_gene_table$GSAVBLPadj <-  NA
  }
  temp_gene_table <- temp_gene_table[,c(2,3,17:27,8:16)]
  if (flag == T){
    final_df <- temp_gene_table
    flag = F
  }
  else{
    final_df <- rbind(final_df, temp_gene_table)
  }
}


write.table(final_df, file = "results/genes_per_modules.tsv", sep = "\t")



#plot heatmap 

#Plot heatmap with the genes and the similarity measure
plot_tom <- dissTOM
diag(plot_tom) <- NA


TOMplot(plot_tom, geneTree, mergedColors, main = "Network heatmap plot, meso vs idio\nDEG genes")

#elevate the plot_tom to a number to create better figure

plot_tom_2 <- plot_tom ^ 20
png(paste(out.dir, "Network heatmap plot.png", sep = "/"))
TOMplot(plot_tom_2, geneTree, mergedColors, main = "Network heatmap plot, meso vs idio\nDEG genes")
dev.off()





##Make heatmap of x  module (sample code)
#royal_blue_genes <- row.names(final_df)[which(final_df$module == "royalblue")]

#using expression levels (after rlog)
#expression_matrix2 <- assay(rld)
#expression_matrix_royalBlue <- subset(expression_matrix2, row.names(expression_matrix2) %in% royal_blue_genes)


#rename the rows with the manual annotated genes
#row.names(expression_matrix_royalBlue)[which(row.names(expression_matrix_royalBlue) == "CRO_127167")] <- "D4H"
#row.names(expression_matrix_royalBlue)[which(row.names(expression_matrix_royalBlue) == "CRO_120021")] <- "DAT"
#row.names(expression_matrix_royalBlue)[which(row.names(expression_matrix_royalBlue) == "CRO_113666")] <- "THAS1"
#row.names(expression_matrix_royalBlue)[which(row.names(expression_matrix_royalBlue) == "CRO_113150")] <- "THAS2s"


#png("test.png", height = 1800, width = 500)
#Heatmap(expression_matrix_royalBlue, clustering_distance_columns = "spearman", clustering_distance_rows = "spearman", rect_gp = gpar(col = "white", lwd = 2))
#dev.off()


##This tables was created after this script was ran 
#temp code
#library(xlsx)
#expression_matrix3 <- expression_matrix2

#royal_blue_genes_table <- read.xlsx(file = "results/genes_per_module.xlsx", sheetName = "royalblue")
#row.names(royal_blue_genes_table) <- royal_blue_genes_table$gene_id 
#royal_blue_genes_table <- royal_blue_genes_table[,-1]


#png("test_2_names.png", height = 1800, width = 4000)
#Heatmap(expression_matrix_royalBlue, clustering_distance_columns = "spearman", clustering_distance_rows = "spearman", rect_gp = gpar(col = "white", lwd = 2), row_labels = royal_blue_genes_table$Swissprot.name)
#dev.off()


##prepare images for final plot
label_id <- paste0("Module ", seq(1:15))
METree_merged$labels <- label_id

png("results/module_eigengene_clustering_final.png", res = 300, width = 2000, height = 1100) 
plot(METree_merged, main = "Hieraquical clustering of merged module eigengens",
     xlab = "", sub = "")
dev.off()

colnames(mergedMEs) <- METree_merged$labels

png("results/eigengene_heatmap_clustering_final.png", res = 300, width = 2200, height = 2000) 
plotEigengeneNetworks(mergedMEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90, signed = TRUE, greyLabel = 0)
dev.off()


