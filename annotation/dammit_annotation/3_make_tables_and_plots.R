#this script will be used to add fold change and p.value to the annotation table and also will 
#be used to make plots for each indidivual gene

#set wrkdir
setwd("C:/Users/Faculdade/Desktop/Dissertação/results 2020/mapping_to_genome/annotation/03_dammit")
#Misc 
rm(list = ls())

#Load libraries 
library(DESeq2)
library("readxl")
library ("xlsx")
library ("ggplot2")
library(gridExtra)
library(EnhancedVolcano)
library(scales)

#make results dir

plot.dir <- "./plots"
if (!dir.exists(plot.dir)){
  dir.create(plot.dir)
}


############################
## Make expressions plots ##
############################

kallisto_normalized_table <- read.table(file = "DESeq2_tables/kallisto_normalizedCountsTable.txt")
kallisto_tpm <- read.table(file = "kallisto_TPM/kallisto_TPM.tsv")

#get subset with the upregulated genes in the 3 conditions 
kallisto_genes <- read.xlsx(file = "results/intersect_3/dammit_annotation_geral_kallisto.xlsx", sheetIndex = 1)
kallisto_genes <- as.character(kallisto_genes$NA.)

#get table without total protoplasts and leaf samples (as these were not used in any comparitons)
kallisto_Up_genes_normalized_table_partial <- kallisto_normalized_table[kallisto_genes,][,-c(16,17,21,22,23)]
kallisto_Up_genes_TPM_table_partial <- kallisto_tpm[kallisto_genes,][-c(19,20,21,22,23)]

#get name vector to use in the title of the plots 

geneName_vector <- c("HSP70", "ABAH4", "Phytosulfokines", "Small HSP", 
                     "FtsZ homolog|BIC1", "CATHA.8139", "CRO_119698", "CATHA.10013", 
                     "PAP17", "NTMC2T1.1", "Calmodulin-like protein", "CRO_125795",
                     "Putative ribonuclease H", "PAP2", "AtPabN2 | ABC transporter", 
                     "Serine/threonine-protein kinase", "AtPSK5", "NFP4.6", "GA2-oxidase", 
                     "HSP18.2", "SCPL8", "HIPP21", "Multicopper oxidase", "RD29A", 
                     "RFP1", "Cytochrome b-c1 complex subunit", "RFS", "rehalose-phosphate phosphatase A",
                     "CrCYP76T24", "ght5", "HSP18.2", "mtnN", "HSP17.4")


##Generate plots with expression, sample name and treatement 

samples <- c("P1F1in", "P2F1in", "P3F1in", "P1F1out", "P2F1out", "P3F1out", "P1F4in", "P2F4in", "P3F4in", "P1F4out", "P2F4out", "P3F4out", "idio_1", "idio_2", "idio_3", "meso_1", "meso_2", "meso_3")
condition <- c("F1in", "F1in", "F1in", "F1out", "F1out", "F1out", "F4in", "F4in", "F4in", "F4out", "F4out", "F4out", "idio", "idio","idio", "meso", "meso", "meso")

#Generate plots with raw counts 
colour_blind_pallete <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#ffcc33", "#0072B2", "#D55E00", "#CC79A7")
names(colour_blind_pallete) <-  c("F1in","F1out", "F4in","F4out", "idio", "leaf",  "meso", "totalptt")



genes_plot = list()
for (i in 1:length(kallisto_genes)){
  #generate the data set by updating the row we pulling the data from
  expression <- as.numeric(kallisto_Up_genes_normalized_table_partial[i,])
  gene_expression_dataFrame <- data.frame("expression" = expression, "samples" = samples, "condition" = condition)
  
  plot <- ggplot(gene_expression_dataFrame, aes(x=condition, y = expression, fill = condition)) + 
    geom_dotplot(binaxis = 'y', stackdir='center', stackratio = 1.2, 
                 dotsize = 1, show.legend = FALSE) + 
    scale_fill_manual(breaks = names(colour_blind_pallete), values=colour_blind_pallete) + 
    scale_x_discrete(limits=c("meso", "idio", "F1in", "F1out", "F4in", "F4out")) + 
    ggtitle(geneName_vector[i]) + 
    theme_classic() + 
    ylab("kallisto normalized expression")
  
  
  genes_plot[[i]] <- ggplotGrob(plot)
   
}


#Generate plot for TPM counts and add to the gene_plot list


for (i in 1:length(kallisto_genes)){
  #generate the data set by updating the row we pulling the data from
  expression <- as.numeric(kallisto_Up_genes_TPM_table_partial[i,])
  gene_expression_dataFrame <- data.frame("expression" = expression, "samples" = samples, "condition" = condition)
  
  plot <- ggplot(gene_expression_dataFrame, aes(x=condition, y = expression, fill = condition)) + 
    geom_dotplot(binaxis = 'y', stackdir='center', stackratio = 1.2, 
                 dotsize = 1, show.legend = FALSE) + 
    scale_fill_manual(breaks = names(colour_blind_pallete), values=colour_blind_pallete) +  
    scale_x_discrete(limits=c("meso", "idio", "F1in", "F1out", "F4in", "F4out")) + 
    ggtitle(geneName_vector[i]) + 
    theme_classic() + 
    ylab("kallisto TPM")
  
  j <- i + 33
  genes_plot[[j]] <- ggplotGrob(plot)
  
}


#MAKE PLOTS, 4 per pdf (i, i+35, i+1, i+36)

for (i in seq(1, 33, 2)){
  pdf_name <- paste("kallisto_gene_expression_part", (i+1)/2, ".pdf", sep  = "") 
  pdf(paste(plot.dir, pdf_name, sep = "/"), width = 10)
  grid.arrange(grobs = genes_plot[c(i, i+33, i+1, i+34)],
               ncol = 2)
  dev.off()
}

pdf(paste(plot.dir, "kallisto_gene_expression_part17.pdf", sep = "/"), width = 10)
grid.arrange(grobs = genes_plot[c(33, 33*2)],
             ncol = 2)
dev.off()



#plot CrPRX1 levels in all samples

colour_blind_pallete_temp <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#ffcc33", "#0072B2", "#D55E00", "#CC79A7")
names(colour_blind_pallete_temp) <-  c("F1in","F1out", "F4in","F4out", "idio", "leaf",  "meso", "totalptt")


samples_temp <- c("P1F1in", "P2F1in", "P3F1in", "P1F1out", "P2F1out","P3F1out", 
                  "P1F4in","P2F4in","P3F4in", "P1F4out","P2F4out","P3F4out", 
                  "idio_1", "idio_2", "idio_3", "leaf_1", "leaf_2",
                  "meso_1", "meso_2","meso_3", "totalptt_1", "totalptt_2", "totalptt_3") 
                  
condition_temp <- c("F1in", "F1in", "F1in", "F1out", "F1out","F1out", 
                    "F4in","F4in","F4in", "F4out","F4out","F4out", 
                    "idio", "idio", "idio", "leaf", "leaf",
                    "meso", "meso","meso", "totalptt", "totalptt", "totalptt")

expression <- as.numeric(kallisto_normalized_table[which(row.names(kallisto_normalized_table) == "CRO_141131"),])
gene_expression_dataFrame <- data.frame("expression" = expression, "samples" = samples_temp, "condition" = condition_temp)


png(paste(plot.dir, "prx_1_expression.png", sep = "/"), res = 300, width = 1500, height = 1000)
ggplot(gene_expression_dataFrame, aes(x=condition_temp, y = expression, fill = condition_temp)) + 
  geom_dotplot(binaxis = 'y', stackdir='center', stackratio = 1.2, 
               dotsize = 1, show.legend = FALSE) + 
  scale_fill_manual(breaks = names(colour_blind_pallete_temp), values=colour_blind_pallete_temp) +  
  scale_x_discrete(limits=c("meso", "idio", "leaf", "totalptt","F1in", "F1out", "F4in", "F4out")) + 
  ggtitle("PRX1 expression level") + 
  theme_classic() + 
  ylab("kallisto normalized Counts")

dev.off()

expression <- as.numeric(kallisto_tpm[which(row.names(kallisto_normalized_table) == "CRO_141131"),])
gene_expression_dataFrame <- data.frame("expression" = expression, "samples" = samples_temp, "condition" = condition_temp)

png(paste(plot.dir, "prx_1_expression_TMP.png", sep = "/"), res = 300, width = 1500, height = 1000)
ggplot(gene_expression_dataFrame, aes(x=condition_temp, y = expression, fill = condition_temp)) + 
  geom_dotplot(binaxis = 'y', stackdir='center', stackratio = 1.2, 
               dotsize = 1, show.legend = FALSE) + 
  scale_fill_manual(breaks = names(colour_blind_pallete_temp), values=colour_blind_pallete_temp) +  
  scale_x_discrete(limits=c("meso", "idio", "leaf", "totalptt","F1in", "F1out", "F4in", "F4out")) + 
  ggtitle("PRX1 expression level") + 
  theme_classic() + 
  ylab("kallisto TPM")


####################################################
#### Make volcano plot with the annotated genes ####
####################################################

####For meso vs idio----

idio_kallisto <- read.table("DESeq2_tables/meso_vs_idio_kallisto_Res.csv", row.names = 1, sep = ",", header = T)
idio_kallisto <- idio_kallisto[complete.cases(idio_kallisto), ]

intersect_3_kallisto <- read.xlsx("results/intersect_3/dammit_annotation_geral_kallisto.xlsx", sheetIndex = 1)
labels_id <- as.character(intersect_3_kallisto$NA.)

#change the pvalue from kallisto to be max -log10(10^-20)
lower_limit <- 10^-20
idio_kallisto$padj <- ifelse(idio_kallisto$padj < lower_limit, lower_limit, idio_kallisto$padj)


color_vector <- rep(NA, length(idio_kallisto$baseMean))
for (i in 1:length(idio_kallisto$padj)){
  if (abs(idio_kallisto$log2FoldChange[i]) > 1 & idio_kallisto$padj[i] < 0.05){
    color_vector[i] <- "tomato3"
  }
  else if (abs(idio_kallisto$log2FoldChange[i]) < 1 & idio_kallisto$padj[i] < 0.05){
    color_vector[i] <- "royalblue"
  }
  else if (abs(idio_kallisto$log2FoldChange[i]) > 1 & idio_kallisto$padj[i] > 0.05){
    color_vector[i] <- "green1"
  }
  else if (abs(idio_kallisto$log2FoldChange[i]) < 1 & idio_kallisto$padj[i] > 0.05){
    color_vector[i] <- "gray30"
  }
}


names(color_vector)[color_vector == 'tomato3'] <- "Fold change > 2 & adjusted p.value < 0.05"
names(color_vector)[color_vector == 'royalblue'] <- "adjusted p.value < 0.05"
names(color_vector)[color_vector == 'green1'] <- "Fold Change > 2"
names(color_vector)[color_vector == 'gray30'] <- "Non sig"


#plot a volcano plot 

#Change the row names for annotation purposes 
which(rownames(idio_kallisto) == "CRO_138972") # 1897
which(rownames(idio_kallisto) == "CRO_131122") # 1750
which(rownames(idio_kallisto) == "CRO_120534") # 2619
new_vector <- rownames(idio_kallisto)
new_vector[1750] <- "NFP 4.6"
new_vector[1897] <- "CYP450 76T24"
#new_vector[2619] <- "PAP17" 
rownames(idio_kallisto) <- new_vector
labels_id[18] <- "NFP 4.6"
labels_id[29] <- "CYP450 76T24"
#labels_id[9] <- "PAP17"


plot <- EnhancedVolcano(idio_kallisto,
                          lab = rownames(idio_kallisto),
                          x = 'log2FoldChange',
                          y = 'padj',
                          xlim = c(-10, 10), 
                          ylim = c(0, 20),
                          xlab = bquote(~Log[2]~ 'fold change'),
                          selectLab = c("NFP 4.6", "CYP450 76T24", "PAP17"),
                          title = 'Meso vs idio',
                          subtitle =  '',
                          pCutoff = 0.05,
                          FCcutoff = 1.0,
                          pointSize = ifelse(row.names(idio_kallisto) %in% labels_id, 1, 1),
                          labSize = 3.1,
                          colCustom = color_vector,
                          gridlines.major = TRUE,
                          gridlines.minor = FALSE,
                          border = 'partial',
                          borderWidth = 1.5,
                          borderColour = 'black',
                          drawConnectors = TRUE,
                          widthConnectors = 0.5, 
                          colConnectors = 'black', 
                          axisLabSize = 22,
                        titleLabSize = 28)

plot2 <- plot + theme_classic() + theme(legend.position="none") + 
  geom_point(data = subset(idio_kallisto, row.names(idio_kallisto) %in% labels_id), 
             aes(x = log2FoldChange, y = -log10(padj), colour = "relevant genes"), 
             color = '#17301C', size = 1.3)

png(paste(plot.dir, "/Volcano_plot_meso_vs_idio_33.png", sep=""), res = 300, width = 2500, height = 2000)
plot2
dev.off()



####For F1in vs F1out---- 

f1_kallisto <- read.table("DESeq2_tables/F1in_vs_F1out_kallisto_Res.csv", row.names = 1, sep = ",", header = T)
f1_kallisto <- f1_kallisto[complete.cases(f1_kallisto), ]

intersect_3_kallisto <- read.xlsx("results/intersect_3/dammit_annotation_geral_kallisto.xlsx", sheetIndex = 1)
labels_id <- as.character(intersect_3_kallisto$NA.)

intersect_F1_F4_kallisto <- read.xlsx("results/F1_F4_intersect/dammit_annotation_F1_F4.xlsx", sheetIndex = 1)
labels_2_id <- as.character(intersect_F1_F4_kallisto$NA.)

#change the pvalue from kallisto to be max -log10(10^-50)
lower_limit <- 10^-20
f1_kallisto$padj <- ifelse(f1_kallisto$padj < lower_limit, lower_limit, f1_kallisto$padj)


color_vector <- rep(NA, length(f1_kallisto$baseMean))
for (i in 1:length(f1_kallisto$padj)){
  if (abs(f1_kallisto$log2FoldChange[i]) > 1 & f1_kallisto$padj[i] < 0.05){
    color_vector[i] <- "tomato3"
  }
  else if (abs(f1_kallisto$log2FoldChange[i]) < 1 & f1_kallisto$padj[i] < 0.05){
    color_vector[i] <- "royalblue"
  }
  else if (abs(f1_kallisto$log2FoldChange[i]) > 1 & f1_kallisto$padj[i] > 0.05){
    color_vector[i] <- "green1"
  }
  else if (abs(f1_kallisto$log2FoldChange[i]) < 1 & f1_kallisto$padj[i] > 0.05){
    color_vector[i] <- "gray30"
  }
}


names(color_vector)[color_vector == 'tomato3'] <- "Fold change > 2 & adjusted p.value < 0.05"
names(color_vector)[color_vector == 'royalblue'] <- "adjusted p.value < 0.05"
names(color_vector)[color_vector == 'green1'] <- "Fold Change > 2"
names(color_vector)[color_vector == 'gray30'] <- "Non sig"


#plot a volcano plot 

#Change the row names for annotation purposes 
which(rownames(f1_kallisto) == "CRO_120534") # 3269
which(rownames(f1_kallisto) == "CRO_131122") # 262
which(rownames(f1_kallisto) == "CRO_138972") # 3598

new_vector <- rownames(f1_kallisto)
new_vector[262] <- "NFP4.6"
new_vector[3598] <- "CYP450 76T24"
new_vector[3269] <- "PAP17"
rownames(f1_kallisto) <- new_vector
labels_id[9] <- "PAP17"
labels_id[18] <- "NFP 4.6"
labels_id[29] <- "CYP450 76T24"




plot <- EnhancedVolcano(f1_kallisto,
                        lab = rownames(f1_kallisto),
                        x = 'log2FoldChange',
                        y = 'padj',
                        xlim = c(-10, 10), 
                        ylim = c(0, 20),
                        xlab = bquote(~Log[2]~ 'fold change'),
                        selectLab = c("NFP4.6", "CYP450 76T24", "PAP17"),
                        title = 'f1in vs f1out',
                        subtitle =  'kallisto',
                        pCutoff = 0.05,
                        FCcutoff = 1.0,
                        pointSize = ifelse(row.names(f1_kallisto) %in% labels_id, 1.3, 1),
                        labSize = 3,
                        colCustom = color_vector,
                        gridlines.major = TRUE,
                        gridlines.minor = FALSE,
                        border = 'partial',
                        borderWidth = 1.5,
                        borderColour = 'black',
                        drawConnectors = TRUE,
                        widthConnectors = 0.5, 
                        colConnectors = 'black')

plot2 <- plot + theme_classic() + theme(legend.position="none") + 
  geom_point(data = subset(f1_kallisto, row.names(f1_kallisto) %in% labels_id), 
             aes(x = log2FoldChange, y = -log10(padj), colour = "relevant genes"), 
             color = '#17301C', size = 1.3)

png(paste(plot.dir, "/Volcano_plot_f1in_vs_f1out_33.png", sep=""), width=12,height=8, units='in', res=500)
plot2
dev.off()




#Change the row names for annotation purposes 
#insert relevant annotation here
which(rownames(f1_kallisto) == "CRO_112949") # 172
new_vector <- rownames(f1_kallisto)
new_vector[172] <- "CrMATE3"
rownames(f1_kallisto) <- new_vector
labels_2_id[189] <- "CrMATE3"



plot <- EnhancedVolcano(f1_kallisto,
                        lab = rownames(f1_kallisto),
                        x = 'log2FoldChange',
                        y = 'padj',
                        xlim = c(-10, 10), 
                        ylim = c(0, 20),
                        xlab = bquote(~Log[2]~ 'fold change'),
                        selectLab = "CrMATE3",
                        title = 'f1in vs f1out',
                        subtitle =  'kallisto',
                        pCutoff = 0.05,
                        FCcutoff = 1.0,
                        pointSize = ifelse(row.names(f1_kallisto) %in% labels_2_id, 1, 1),
                        labSize = 3,
                        colCustom = color_vector,
                        gridlines.major = TRUE,
                        gridlines.minor = FALSE,
                        border = 'partial',
                        borderWidth = 1.5,
                        borderColour = 'black',
                        drawConnectors = TRUE,
                        widthConnectors = 0.5, 
                        colConnectors = 'black')

plot2 <- plot + theme_classic() + theme(legend.position="none") + 
  geom_point(data = subset(f1_kallisto, row.names(f1_kallisto) %in% labels_2_id), 
             aes(x = log2FoldChange, y = -log10(padj), colour = "relevant genes"), 
             color = '#17301C', size = 1)

png(paste(plot.dir, "/Volcano_plot_f1in_vs_f1out_452.png", sep=""), width=12,height=8, units='in', res=500)
plot2
dev.off()



####For f4in vs f4out---- 

f4_kallisto <- read.table("DESeq2_tables/f4in_vs_f4out_kallisto_Res.csv", row.names = 1, sep = ",", header = T)
f4_kallisto <- f4_kallisto[complete.cases(f4_kallisto), ]

intersect_3_kallisto <- read.xlsx("results/intersect_3/dammit_annotation_geral_kallisto.xlsx", sheetIndex = 1)
labels_id <- as.character(intersect_3_kallisto$NA.)

intersect_F1_F4_kallisto <- read.xlsx("results/F1_F4_intersect/dammit_annotation_F1_F4.xlsx", sheetIndex = 1)
labels_2_id <- as.character(intersect_F1_F4_kallisto$NA.)

#change the pvalue from kallisto to be max -log10(10^-50)
lower_limit <- 10^-20
f4_kallisto$padj <- ifelse(f4_kallisto$padj < lower_limit, lower_limit, f4_kallisto$padj)


color_vector <- rep(NA, length(f4_kallisto$baseMean))
for (i in 1:length(f4_kallisto$padj)){
  if (abs(f4_kallisto$log2FoldChange[i]) > 1 & f4_kallisto$padj[i] < 0.05){
    color_vector[i] <- "tomato3"
  }
  else if (abs(f4_kallisto$log2FoldChange[i]) < 1 & f4_kallisto$padj[i] < 0.05){
    color_vector[i] <- "royalblue"
  }
  else if (abs(f4_kallisto$log2FoldChange[i]) > 1 & f4_kallisto$padj[i] > 0.05){
    color_vector[i] <- "green1"
  }
  else if (abs(f4_kallisto$log2FoldChange[i]) < 1 & f4_kallisto$padj[i] > 0.05){
    color_vector[i] <- "gray30"
  }
}


names(color_vector)[color_vector == 'tomato3'] <- "Fold change > 2 & adjusted p.value < 0.05"
names(color_vector)[color_vector == 'royalblue'] <- "adjusted p.value < 0.05"
names(color_vector)[color_vector == 'green1'] <- "Fold Change > 2"
names(color_vector)[color_vector == 'gray30'] <- "Non sig"


#plot a volcano plot 

#Change the row names for annotation purposes 
which(rownames(f4_kallisto) == "CRO_120534") # 6350
which(rownames(f4_kallisto) == "CRO_131122") # 1377
which(rownames(f4_kallisto) == "CRO_138972") # 5464

new_vector <- rownames(f4_kallisto)
new_vector[1377] <- "NFP4.6"
new_vector[5464] <- "CYP450 76T24"
new_vector[6350] <- "PAP17"
rownames(f4_kallisto) <- new_vector
labels_id[9] <- "PAP17"
labels_id[18] <- "NFP 4.6"
labels_id[29] <- "CYP450 76T24"



plot <- EnhancedVolcano(f4_kallisto,
                        lab = rownames(f4_kallisto),
                        x = 'log2FoldChange',
                        y = 'padj',
                        xlim = c(-10, 10), 
                        ylim = c(0, 20),
                        xlab = bquote(~Log[2]~ 'fold change'),
                        selectLab = c("PAP17", "NFP4.6", "CYP450 76T24"),
                        title = 'f4in vs f4out',
                        subtitle =  'kallisto',
                        pCutoff = 0.05,
                        FCcutoff = 1.0,
                        pointSize = ifelse(row.names(f4_kallisto) %in% labels_id, 1, 1),
                        labSize = 4,
                        colCustom = color_vector,
                        gridlines.major = TRUE,
                        gridlines.minor = FALSE,
                        border = 'partial',
                        borderWidth = 1.5,
                        borderColour = 'black',
                        drawConnectors = TRUE,
                        widthConnectors = 0.5, 
                        colConnectors = 'black')

plot2 <- plot + theme_classic() + theme(legend.position="none") + 
  geom_point(data = subset(f4_kallisto, row.names(f4_kallisto) %in% labels_id), 
             aes(x = log2FoldChange, y = -log10(padj), colour = "relevant genes"), 
             color = '#17301C', size = 1.3)

png(paste(plot.dir, "/Volcano_plot_f4in_vs_f4out_33.png", sep=""), width=12,height=8, units='in', res=500)
plot2
dev.off()



#Change the row names for annotation purposes 
#insert relevant annotation here
which(rownames(f4_kallisto) == "CRO_112949") # 64
new_vector <- rownames(f4_kallisto)
new_vector[64] <- "CrMate3"
rownames(f4_kallisto) <- new_vector
labels_2_id[190] <- "CrMate3"



plot <- EnhancedVolcano(f4_kallisto,
                        lab = rownames(f4_kallisto),
                        x = 'log2FoldChange',
                        y = 'padj',
                        xlim = c(-10, 10), 
                        ylim = c(0, 20),
                        xlab = bquote(~Log[2]~ 'fold change'),
                        selectLab = "CrMate3",
                        title = 'f4in vs f4out',
                        subtitle =  'kallisto',
                        pCutoff = 0.05,
                        FCcutoff = 1.0,
                        pointSize = ifelse(row.names(f4_kallisto) %in% labels_2_id, 1, 1),
                        labSize = 3,
                        colCustom = color_vector,
                        gridlines.major = TRUE,
                        gridlines.minor = FALSE,
                        border = 'partial',
                        borderWidth = 1.5,
                        borderColour = 'black',
                        drawConnectors = TRUE,
                        widthConnectors = 0.5, 
                        colConnectors = 'black')

plot2 <- plot + theme_classic() + theme(legend.position="none") + 
  geom_point(data = subset(f4_kallisto, row.names(f4_kallisto) %in% labels_2_id), 
             aes(x = log2FoldChange, y = -log10(padj), colour = "relevant genes"), 
             color = '#17301C', size = 1)

png(paste(plot.dir, "/Volcano_plot_f4in_vs_f4out_452.png", sep=""), width=12,height=8, units='in', res=500)
plot2
dev.off()




###############################################
#### Make MA plot with the annotated genes ####
###############################################

#for meso vs idio ----
idio_kallisto <- read.csv("DESeq2_tables/meso_vs_idio_kallisto_Res.csv", row.names = 1)

#Build a dataframe with the significance of each gene (upregulated, downregulated, none, not tested)
idio_kallisto_df.MA <- idio_kallisto
idio_kallisto_df.MA$sig <- c()
idio_kallisto_df.MA$sig <- as.factor(ifelse(is.na(idio_kallisto_df.MA$padj), "Not tested",
                                  ifelse(idio_kallisto_df.MA$padj > 0.05, "Non Sig", 
                                  ifelse(idio_kallisto$log2FoldChange >= 1, "up", 
                                  ifelse(idio_kallisto_df.MA$log2FoldChange <= -1, "down", "Non Sig")))))


y_upper_limit <- 5
y_lower_limit <- -5

idio_kallisto_df.MA$shape <- ifelse(idio_kallisto_df.MA$log2FoldChange >= y_upper_limit, "Tup", ifelse(idio_kallisto_df.MA$log2FoldChange <= y_lower_limit, "Tdown", "normal"))
idio_kallisto_df.MA$adjLog2change <- ifelse(idio_kallisto_df.MA$log2FoldChange >= y_upper_limit, y_upper_limit, ifelse(idio_kallisto_df.MA$log2FoldChange <= y_lower_limit, y_lower_limit, idio_kallisto_df.MA$log2FoldChange))


MAplot <- ggplot(idio_kallisto_df.MA, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.5) + 
  theme_classic() + 
  ggtitle("Idio vs meso") + xlab("Mean Expression") + ylab("Log2FoldChange") + 
  theme(plot.title = element_text(size=22, face = "bold"), 
        axis.text=element_text(size=12))


png(paste(plot.dir, "base_MA_plot_idio_vs_meso.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MAplot
dev.off()

#plot pathway specific MA plots 
manual_annotation <- read.xlsx(file = "StepsToFinalAnnotation/manual_annotation_modified.xlsx", sheetIndex = 1)


#iridoid Pathway
idio_kallisto_df.iridoid <- idio_kallisto_df.MA
iridoid_pathway_genes <- manual_annotation[1:6, c(2,4)]
iridoid_pathway_genes$Enzyme <- c("G10H", "10HGO", "IS", "IO", "7DLGT", "7DLH")


#make manual annotation collum
idio_kallisto_df.iridoid$label <- row.names(idio_kallisto_df.iridoid)
for (i in 1:length(iridoid_pathway_genes$Gene_id)){
  car_id <- as.character(iridoid_pathway_genes$Gene_id[i])
  idio_kallisto_df.iridoid[car_id,10] <- as.character(iridoid_pathway_genes$Enzyme[i])
}

up_label_genes <- iridoid_pathway_genes$Gene_id[c(3,5)]
down_label_genes <- iridoid_pathway_genes$Gene_id[-c(3,5)]

MA_iridoid <- ggplot(idio_kallisto_df.iridoid, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("Idio vs meso", subtitle = "Iridoid Pathway") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(ylim = c(-5,-3), box.padding= 0.5, size = 6, aes(label = ifelse(row.names(idio_kallisto_df.iridoid) %in% down_label_genes, label, ""))) +
  geom_label_repel(ylim = c(2,3), box.padding= 0.5, size = 6, aes(label = ifelse(row.names(idio_kallisto_df.iridoid) %in% up_label_genes, label, ""))) + 
  geom_point(data = subset(idio_kallisto_df.iridoid, row.names(idio_kallisto_df.iridoid) %in% iridoid_pathway_genes$Gene_id), 
           aes(x = baseMean, y = adjLog2change), color = "black", 
           size = 1.8, shape = 21, stroke = 1.5) + 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18), axis.text=element_text(size=12))


png(paste(plot.dir, "iridoid_MA_plot_idio_vs_meso.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MA_iridoid
dev.off()


# #central Pathway
# idio_kallisto_df.central <- idio_kallisto_df.MA
# central_pathway_genes <- manual_annotation[c(7:11,19:23,27:32,73,77:79), c(2,4)]
# central_pathway_genes$Enzyme <- c("LAMT","SLS2", "TDC", "STR", "SGD", "ASO", "DPAS", "GS1","CS","TS", "Asa", "GO_1","GO_2", "Redox 1", "Redox 2", "SAT", "GS2", "SLS1", "SLS3", "SLS3")
# 
# 
# #make manual annotation collum
# idio_kallisto_df.central$label <- row.names(idio_kallisto_df.central)
# for (i in 1:length(central_pathway_genes$Gene_id)){
#   car_id <- as.character(central_pathway_genes$Gene_id[i])
#   idio_kallisto_df.central[car_id,10] <- as.character(central_pathway_genes$Enzyme[i])
# }
# 
# 
# MA_central <- ggplot(idio_kallisto_df.central, aes(x = baseMean, y = adjLog2change, group = shape))+ 
#   geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
#   scale_shape_manual(values=c(19, 6, 2)) + 
#   scale_color_manual(values = c("green1", "gray30","grey", "tomato3")) + 
#   scale_x_continuous(trans = 'log10', 
#                      breaks = trans_breaks("log10", function(x) 10^x),
#                      labels = trans_format("log10", math_format(10^.x))) + 
#   geom_hline(yintercept = 0, color = "black", size = 0.3) + 
#   theme_classic() + 
#   ggtitle("Idio vs meso", subtitle =  "Central Pathway") + xlab("Mean Expression") + ylab("Log2FoldChange") +
#   geom_label_repel(box.padding= 0.51, size = 6, aes(label = ifelse(row.names(idio_kallisto_df.central) %in% central_pathway_genes$Gene_id, label, ""))) + 
#   geom_point(data = subset(idio_kallisto_df.central, row.names(idio_kallisto_df.central) %in% central_pathway_genes$Gene_id), 
#              aes(x = baseMean, y = adjLog2change), color = "black", 
#              size = 1.8, shape = 21, stroke = 1.5)+ 
#   theme(plot.title = element_text(size=22, face = "bold"), 
#         plot.subtitle = element_text(size=18),   axis.text=element_text(size=12))
# 
# 
# 
# png(paste(plot.dir, "central_MA_plot_idio_vs_meso.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
# MA_central
# dev.off()

##Make first plot of the central pathway 
idio_kallisto_df.central_1 <- idio_kallisto_df.MA
central_1_pathway_genes <- manual_annotation[c(7:11,27,77:79), c(2,4)]
central_1_pathway_genes$Enzyme <- c("LAMT","SLS2", "TDC", "STR", "SGD", "Asa", "SLS1", "SLS3", "SLS4")


#make manual annotation collum
idio_kallisto_df.central_1$label <- row.names(idio_kallisto_df.central_1)
for (i in 1:length(central_1_pathway_genes$Gene_id)){
  car_id <- as.character(central_1_pathway_genes$Gene_id[i])
  idio_kallisto_df.central_1[car_id,10] <- as.character(central_1_pathway_genes$Enzyme[i])
}


MA_central_1 <- ggplot(idio_kallisto_df.central_1, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("idioin vs idioout", subtitle =  "central Pathway") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(box.padding= 2, size = 6, aes(label = ifelse(row.names(idio_kallisto_df.central_1) %in% central_1_pathway_genes$Gene_id, label, ""))) + 
  geom_point(data = subset(idio_kallisto_df.central_1, row.names(idio_kallisto_df.central_1) %in% central_1_pathway_genes$Gene_id), 
             aes(x = baseMean, y = adjLog2change), color = "black", 
             size = 1.8, shape = 21, stroke = 1.5)+ 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18), axis.text=element_text(size=12))



png(paste(plot.dir, "central_1_MA_plot_idioin_vs_idioout.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MA_central_1
dev.off()

#Plot the second part of the central pathway
idio_kallisto_df.central_2 <- idio_kallisto_df.MA
central_2_pathway_genes <- manual_annotation[c(19:23,28:32,73), c(2,4)]
central_2_pathway_genes$Enzyme <- c("ASO", "DPAS", "GS1","CS","TS", "GO_1","GO_2", "Redox 1", "Redox 2", "SAT", "GS2")


idio_kallisto_df.central_2$label <- row.names(idio_kallisto_df.central_2)
for (i in 1:length(central_2_pathway_genes$Gene_id)){
  car_id <- as.character(central_2_pathway_genes$Gene_id[i])
  idio_kallisto_df.central_2[car_id,10] <- as.character(central_2_pathway_genes$Enzyme[i])
}


MA_central_2 <- ggplot(idio_kallisto_df.central_2, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("idioin vs idioout", subtitle =  "central Pathway") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(box.padding= 3.1, size = 6, aes(label = ifelse(row.names(idio_kallisto_df.central_2) %in% central_2_pathway_genes$Gene_id, label, ""))) + 
  geom_point(data = subset(idio_kallisto_df.central_2, row.names(idio_kallisto_df.central_2) %in% central_2_pathway_genes$Gene_id), 
             aes(x = baseMean, y = adjLog2change), color = "black", 
             size = 1.8, shape = 21, stroke = 1.5)+ 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18), axis.text=element_text(size=12))



png(paste(plot.dir, "central_2_MA_plot_idioin_vs_idioout.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MA_central_2
dev.off()

#tamberstonine to vindoline pathway
idio_kallisto_df.vindoline <- idio_kallisto_df.MA

vindoline_pathway_genes <- manual_annotation[c(12:18,33, 76), c(2,4)]
vindoline_pathway_genes$Enzyme <- c("T16H2","16OMT", "T3O", "T3R", "16NMT", "D4H", "DAT", "PRX1", "T16H1")

#make manual annotation collum
idio_kallisto_df.vindoline$label <- row.names(idio_kallisto_df.vindoline)
for (i in 1:length(vindoline_pathway_genes$Gene_id)){
  car_id <- as.character(vindoline_pathway_genes$Gene_id[i])
  idio_kallisto_df.vindoline[car_id,10] <- as.character(vindoline_pathway_genes$Enzyme[i])
}

#Remove two genes as they are over 5 foldchange 
vindoline_pathway_genes_cut <- vindoline_pathway_genes[-c(1,7),]


MA_vindoline <- ggplot(idio_kallisto_df.vindoline, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("Idio vs meso", subtitle =  "Vindoline Pathway") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(box.padding= 2, size = 6, aes(label = ifelse(row.names(idio_kallisto_df.vindoline) %in% vindoline_pathway_genes$Gene_id, label, ""))) + 
  geom_point(data = subset(idio_kallisto_df.vindoline, row.names(idio_kallisto_df.vindoline) %in% vindoline_pathway_genes_cut$Gene_id), 
             aes(x = baseMean, y = adjLog2change), color = "black", 
             size = 1.8, shape = 21, stroke = 1.5)+ 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18),   axis.text=element_text(size=12))



png(paste(plot.dir, "vindoline_MA_plot_idio_vs_meso.png", sep = "/"), width = 2000, height = 1800, res = 300)
MA_vindoline
dev.off()


#transporters 
idio_kallisto_df.transporters <- idio_kallisto_df.MA

transporters_pathway_genes <- manual_annotation[c(61:72), c(2,4)]
transporters_pathway_genes$Enzyme <- c("CrTPT2", "CrNPF2.9", "CrNPF2.4", "CrNPF2.5", "CrNPF2.6", "MATE1_KX372304", "MATE2_KX372305", "CrMate2PlantBio", "CrMate3PlantBio", "CrMate4PlantBio", "CrMate5PlantBio", "CrMate6PlantBio")

#make manual annotation collum
idio_kallisto_df.transporters$label <- row.names(idio_kallisto_df.transporters)
for (i in 1:length(transporters_pathway_genes$Gene_id)){
  car_id <- as.character(transporters_pathway_genes$Gene_id[i])
  idio_kallisto_df.transporters[car_id,10] <- as.character(transporters_pathway_genes$Enzyme[i])
}

up_up_label_genes <- transporters_pathway_genes$Gene_id[c(1,10)]
up_label_genes <- transporters_pathway_genes$Gene_id[c(2,4,6,11)]
mate_2 <- transporters_pathway_genes$Gene_id[c(7)]
down_label_genes <- transporters_pathway_genes$Gene_id[-c(1,2,4,6,7,10,11)]

MA_transporters <- ggplot(idio_kallisto_df.transporters, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 2) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("Idio vs meso", subtitle =  "Transporters") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(ylim = c(4,5), box.padding= 3.3, size = 6, aes(label = ifelse(row.names(idio_kallisto_df.transporters) %in% up_up_label_genes, label, ""))) + 
  geom_label_repel(ylim =c(1,3), box.padding= 1, size = 6, aes(label = ifelse(row.names(idio_kallisto_df.transporters) %in% up_label_genes, label, ""))) +
  geom_label_repel(xlim = c(2, 3), ylim =c(-2.4, -2.5), box.padding = 1, size=6, aes(label = ifelse(row.names(idio_kallisto_df.transporters) %in% mate_2, label, ""))) + 
  geom_label_repel(box.padding= 2.1, size = 6, aes(label = ifelse(row.names(idio_kallisto_df.transporters) %in% down_label_genes, label, ""))) +
  geom_point(data = subset(idio_kallisto_df.transporters, row.names(idio_kallisto_df.transporters) %in% transporters_pathway_genes$Gene_id), 
             aes(x = baseMean, y = adjLog2change), color = "black", 
             size = 2.3, shape = 21, stroke = 0.9)+ 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18),   axis.text=element_text(size=12))



png(paste(plot.dir, "transporters_MA_plot_idio_vs_meso.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MA_transporters
dev.off()

#TFs
idio_kallisto_df.Tfactors <- idio_kallisto_df.MA

Tfactors_pathway_genes <- manual_annotation[c(43:60), c(2,4)]
Tfactors_pathway_genes$Enzyme <- c( "Bis1", "Bis2", "ORCA2", "ORCA3", "ORCA4", "ORCA5", "ORCA6", "CR1", "WRKY1", "MYC2", "ZCT1", "ZCT2", "ZCT3", "GBF1", "GBF2", "BPF1", "CrGATA1", "ERF5")

#make manual annotation collum
idio_kallisto_df.Tfactors$label <- row.names(idio_kallisto_df.Tfactors)
for (i in 1:length(Tfactors_pathway_genes$Gene_id)){
  car_id <- as.character(Tfactors_pathway_genes$Gene_id[i])
  idio_kallisto_df.Tfactors[car_id,10] <- as.character(Tfactors_pathway_genes$Enzyme[i])
}


MA_Tfactors <- ggplot(idio_kallisto_df.Tfactors, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("Idio vs meso", subtitle =  "Trancription factors") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(box.padding= 1.2, size = 5, aes(label = ifelse(row.names(idio_kallisto_df.Tfactors) %in% Tfactors_pathway_genes$Gene_id, label, ""))) + 
  geom_point(data = subset(idio_kallisto_df.Tfactors, row.names(idio_kallisto_df.Tfactors) %in% Tfactors_pathway_genes$Gene_id), 
             aes(x = baseMean, y = adjLog2change), color = "black", 
             size = 1.8, shape = 21, stroke = 1)+ 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18),
        axis.text=element_text(size=12))



png(paste(plot.dir, "Tfactors_MA_plot_idio_vs_meso.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MA_Tfactors
dev.off()

#secondary pathway
idio_kallisto_df.secPathways <- idio_kallisto_df.MA

secPathways_pathway_genes <- manual_annotation[c(24:26,34:42, 75, 80,81), c(2,4)]
secPathways_pathway_genes$Enzyme <- c("HL3", "HL4", "As", "T19H", "TEX1", "TAT", "MAT", "HYS", "VAS_1", "VAS_2", "THAS", "V19H", "TEX2", "THAS2", "THAS3")

#make manual annotation collum
idio_kallisto_df.secPathways$label <- row.names(idio_kallisto_df.secPathways)
for (i in 1:length(secPathways_pathway_genes$Gene_id)){
  car_id <- as.character(secPathways_pathway_genes$Gene_id[i])
  idio_kallisto_df.secPathways[car_id,10] <- as.character(secPathways_pathway_genes$Enzyme[i])
}

#remove genes with the logfoldchange over 5 
secPathways_pathway_genes_cut <- secPathways_pathway_genes[-c(11),]


MA_secPathways <- ggplot(idio_kallisto_df.secPathways, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("Idio vs meso", subtitle =  "Secondary pathways enzymes") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(box.padding= 2, size = 6, aes(label = ifelse(row.names(idio_kallisto_df.secPathways) %in% secPathways_pathway_genes$Gene_id, label, ""))) + 
  geom_point(data = subset(idio_kallisto_df.secPathways, row.names(idio_kallisto_df.secPathways) %in% secPathways_pathway_genes_cut$Gene_id), 
             aes(x = baseMean, y = adjLog2change), color = "black", 
             size = 1.8, shape = 21, stroke = 1)+ 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18),   axis.text=element_text(size=12))


png(paste(plot.dir, "secPathways_MA_plot_idio_vs_meso.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MA_secPathways
dev.off()



#for f1in vs F1out ----
f1_kallisto <- read.csv("DESeq2_tables/F1in_vs_F1out_kallisto_Res.csv", row.names = 1)

#Build a dataframe with the significance of each gene (upregulated, downregulated, none)
f1_kallisto_df.MA <- f1_kallisto
f1_kallisto_df.MA$sig <- c()
f1_kallisto_df.MA$sig <- as.factor(ifelse(is.na(f1_kallisto_df.MA$padj), "Not tested",
                                  ifelse(f1_kallisto_df.MA$padj > 0.05, "Non Sig", 
                                  ifelse(f1_kallisto$log2FoldChange >= 1, "up", 
                                         ifelse(f1_kallisto_df.MA$log2FoldChange <= -1, "down", "Non Sig")))))


y_upper_limit <- 5
y_lower_limit <- -5

f1_kallisto_df.MA$shape <- ifelse(f1_kallisto_df.MA$log2FoldChange >= y_upper_limit, "Tup", ifelse(f1_kallisto_df.MA$log2FoldChange <= y_lower_limit, "Tdown", "normal"))
f1_kallisto_df.MA$adjLog2change <- ifelse(f1_kallisto_df.MA$log2FoldChange >= y_upper_limit, y_upper_limit, ifelse(f1_kallisto_df.MA$log2FoldChange <= y_lower_limit, y_lower_limit, f1_kallisto_df.MA$log2FoldChange))


MAplot <- ggplot(f1_kallisto_df.MA, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("F1in vs F1out") + xlab("Mean Expression") + ylab("Log2FoldChange") + 
  theme(plot.title = element_text(size=22, face = "bold"), axis.text=element_text(size=12))


png(paste(plot.dir, "base_MA_plot_f1in_vs_f1out.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MAplot
dev.off()


#plot pathway specific MA plots 
manual_annotation <- read.xlsx(file = "StepsToFinalAnnotation/manual_annotation_modified.xlsx", sheetIndex = 1)


#iridoid Pathway
f1_kallisto_df.iridoid <- f1_kallisto_df.MA
iridoid_pathway_genes <- manual_annotation[1:6, c(2,4)]
iridoid_pathway_genes$Enzyme<- c("G10H", "10HGO", "IS", "IO", "7DLGT", "7DLH")


#make manual annotation collum
f1_kallisto_df.iridoid$label <- row.names(f1_kallisto_df.iridoid)
for (i in 1:length(iridoid_pathway_genes$Gene_id)){
  car_id <- as.character(iridoid_pathway_genes$Gene_id[i])
  f1_kallisto_df.iridoid[car_id,10] <- as.character(iridoid_pathway_genes$Enzyme[i])
}


MA_iridoid <- ggplot(f1_kallisto_df.iridoid, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("F1in vs F1out", subtitle = "Iridoid Pathway") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(box.padding= .5, size = 6, aes(label = ifelse(row.names(f1_kallisto_df.iridoid) %in% iridoid_pathway_genes$Gene_id, label, ""))) + 
  geom_point(data = subset(f1_kallisto_df.iridoid, row.names(f1_kallisto_df.iridoid) %in% iridoid_pathway_genes$Gene_id), 
             aes(x = baseMean, y = adjLog2change), color = "black", 
             size = 1.8, shape = 21, stroke = 1.5) + 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18), axis.text=element_text(size=12))


png(paste(plot.dir, "iridoid_MA_plot_f1in_vs_f1out.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MA_iridoid
dev.off()


# #central Pathway
# #central pathway had to be divided into two parts
# f1_kallisto_df.central <- f1_kallisto_df.MA
# central_pathway_genes <- manual_annotation[c(7:11,19:23,27:32,73,77:79), c(2,4)]
# central_pathway_genes$Enzyme <- c("LAMT","SLS2", "TDC", "STR", "SGD", "ASO", "DPAS", "GS1","CS","TS", "Asa", "GO_1","GO_2", "Redox 1", "Redox 2", "SAT", "GS2", "SLS1", "SLS3", "SLS3")
# 
# 
# #make manual annotation collum
# f1_kallisto_df.central$label <- row.names(f1_kallisto_df.central)
# for (i in 1:length(central_pathway_genes$Gene_id)){
#   car_id <- as.character(central_pathway_genes$Gene_id[i])
#   f1_kallisto_df.central[car_id,10] <- as.character(central_pathway_genes$Enzyme[i])
# }
# 
# 
# MA_central <- ggplot(f1_kallisto_df.central, aes(x = baseMean, y = adjLog2change, group = shape))+ 
#   geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
#   scale_shape_manual(values=c(19, 6, 2)) + 
#   scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
#   scale_x_continuous(trans = 'log10', 
#                      breaks = trans_breaks("log10", function(x) 10^x),
#                      labels = trans_format("log10", math_format(10^.x))) + 
#   geom_hline(yintercept = 0, color = "black", size = 0.3) + 
#   theme_classic() + 
#   ggtitle("F1in vs F1out", subtitle =  "Central Pathway") + xlab("Mean Expression") + ylab("Log2FoldChange") +
#   geom_label_repel(box.padding= 2, size = 6, aes(label = ifelse(row.names(f1_kallisto_df.central) %in% central_pathway_genes$Gene_id, label, ""))) + 
#   geom_point(data = subset(f1_kallisto_df.central, row.names(f1_kallisto_df.central) %in% central_pathway_genes$Gene_id), 
#              aes(x = baseMean, y = adjLog2change), color = "black", 
#              size = 1.8, shape = 21, stroke = 1.5)+ 
#   theme(plot.title = element_text(size=22, face = "bold"), 
#         plot.subtitle = element_text(size=18), axis.text=element_text(size=12))



#png(paste(plot.dir, "central_MA_plot_f1in_vs_f1out.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
#MA_central
#dev.off()


#Plot the first part of the central pathway
f1_kallisto_df.central_1 <- f1_kallisto_df.MA
central_1_pathway_genes <- manual_annotation[c(7:11,27,77:79), c(2,4)]
central_1_pathway_genes$Enzyme <- c("LAMT","SLS2", "TDC", "STR", "SGD", "Asa", "SLS1", "SLS3", "SLS4")


#make manual annotation collum
f1_kallisto_df.central_1$label <- row.names(f1_kallisto_df.central_1)
for (i in 1:length(central_1_pathway_genes$Gene_id)){
  car_id <- as.character(central_1_pathway_genes$Gene_id[i])
  f1_kallisto_df.central_1[car_id,10] <- as.character(central_1_pathway_genes$Enzyme[i])
}


MA_central_1 <- ggplot(f1_kallisto_df.central_1, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("F1in vs F1out", subtitle =  "central Pathway") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(box.padding= 2, size = 6, aes(label = ifelse(row.names(f1_kallisto_df.central_1) %in% central_1_pathway_genes$Gene_id, label, ""))) + 
  geom_point(data = subset(f1_kallisto_df.central_1, row.names(f1_kallisto_df.central_1) %in% central_1_pathway_genes$Gene_id), 
             aes(x = baseMean, y = adjLog2change), color = "black", 
             size = 1.8, shape = 21, stroke = 1.5)+ 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18), axis.text=element_text(size=12))



png(paste(plot.dir, "central_1_MA_plot_f1in_vs_f1out.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MA_central_1
dev.off()

#Plot the second part of the central pathway
f1_kallisto_df.central_2 <- f1_kallisto_df.MA
central_2_pathway_genes <- manual_annotation[c(19:23,28:32,73), c(2,4)]
central_2_pathway_genes$Enzyme <- c("ASO", "DPAS", "GS1","CS","TS", "GO_1","GO_2", "Redox 1", "Redox 2", "SAT", "GS2")


f1_kallisto_df.central_2$label <- row.names(f1_kallisto_df.central_2)
for (i in 1:length(central_2_pathway_genes$Gene_id)){
  car_id <- as.character(central_2_pathway_genes$Gene_id[i])
  f1_kallisto_df.central_2[car_id,10] <- as.character(central_2_pathway_genes$Enzyme[i])
}


MA_central_2 <- ggplot(f1_kallisto_df.central_2, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("F1in vs F1out", subtitle =  "central Pathway") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(box.padding= 3.1, size = 6, aes(label = ifelse(row.names(f1_kallisto_df.central_2) %in% central_2_pathway_genes$Gene_id, label, ""))) + 
  geom_point(data = subset(f1_kallisto_df.central_2, row.names(f1_kallisto_df.central_2) %in% central_2_pathway_genes$Gene_id), 
             aes(x = baseMean, y = adjLog2change), color = "black", 
             size = 1.8, shape = 21, stroke = 1.5)+ 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18), axis.text=element_text(size=12))



png(paste(plot.dir, "central_2_MA_plot_f1in_vs_f1out.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MA_central_2
dev.off()



#tamberstonine to vindoline pathway
f1_kallisto_df.vindoline <- f1_kallisto_df.MA

vindoline_pathway_genes <- manual_annotation[c(12:18,33, 76), c(2,4)]
vindoline_pathway_genes$Enzyme <- c("T16H2","16OMT", "T3O", "T3R", "16NMT", "D4H", "DAT", "PRX1", "T16H1")

#make manual annotation collum
f1_kallisto_df.vindoline$label <- row.names(f1_kallisto_df.vindoline)
for (i in 1:length(vindoline_pathway_genes$Gene_id)){
  car_id <- as.character(vindoline_pathway_genes$Gene_id[i])
  f1_kallisto_df.vindoline[car_id,10] <- as.character(vindoline_pathway_genes$Enzyme[i])
}



MA_vindoline <- ggplot(f1_kallisto_df.vindoline, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("F1in vs F1out", subtitle =  "Vindoline Pathway") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(box.padding= 2, size = 6, aes(label = ifelse(row.names(f1_kallisto_df.vindoline) %in% vindoline_pathway_genes$Gene_id, label, ""))) + 
  geom_point(data = subset(f1_kallisto_df.vindoline, row.names(f1_kallisto_df.vindoline) %in% vindoline_pathway_genes$Gene_id), 
             aes(x = baseMean, y = adjLog2change), color = "black", 
             size = 1.8, shape = 21, stroke = 1.5)+ 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18), axis.text=element_text(size=12))



png(paste(plot.dir, "vindoline_MA_plot_f1in_vs_f1out.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MA_vindoline
dev.off()


#transporters 
f1_kallisto_df.transporters <- f1_kallisto_df.MA

transporters_pathway_genes <- manual_annotation[c(61:72), c(2,4)]
transporters_pathway_genes$Enzyme <- c("CrTPT2", "CrNPF2.9", "CrNPF2.4", "CrNPF2.5", "CrNPF2.6", "MATE1_KX372304", "MATE2_KX372305", "CrMate2PlantBio", "CrMate3PlantBio", "CrMate4PlantBio", "CrMate5PlantBio", "CrMate6PlantBio")

#make manual annotation collum
f1_kallisto_df.transporters$label <- row.names(f1_kallisto_df.transporters)
for (i in 1:length(transporters_pathway_genes$Gene_id)){
  car_id <- as.character(transporters_pathway_genes$Gene_id[i])
  f1_kallisto_df.transporters[car_id,10] <- as.character(transporters_pathway_genes$Enzyme[i])
}

Mate4 <- transporters_pathway_genes$Gene_id[10]
remaining <- transporters_pathway_genes$Gene_id[-10]

MA_transporters <- ggplot(f1_kallisto_df.transporters, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 2) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("F1in vs F1out", subtitle =  "Transporters") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(box.padding= 2.1, size = 6, aes(label = ifelse(row.names(f1_kallisto_df.transporters) %in% remaining, label, ""))) + 
  geom_label_repel(xlim = c(4.5, 5), ylim = c(-2.5,-2.6), box.padding= 2.1, size = 6, aes(label = ifelse(row.names(f1_kallisto_df.transporters) %in% Mate4, label, ""))) + 
  geom_point(data = subset(f1_kallisto_df.transporters, row.names(f1_kallisto_df.transporters) %in% transporters_pathway_genes$Gene_id), 
             aes(x = baseMean, y = adjLog2change), color = "black", 
             size = 2, shape = 21, stroke = 0.9)+ 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18), axis.text=element_text(size=12))



png(paste(plot.dir, "transporters_MA_plot_f1in_vs_f1out.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MA_transporters
dev.off()

#TFs
f1_kallisto_df.Tfactors <- f1_kallisto_df.MA

Tfactors_pathway_genes <- manual_annotation[c(43:60), c(2,4)]
Tfactors_pathway_genes$Enzyme <- c( "Bis1", "Bis2", "ORCA2", "ORCA3", "ORCA4", "ORCA5", "ORCA6", "CR1", "WRKY1", "MYC2", "ZCT1", "ZCT2", "ZCT3", "GBF1", "GBF2", "BPF1", "CrGATA1", "ERF5")

#make manual annotation collum
f1_kallisto_df.Tfactors$label <- row.names(f1_kallisto_df.Tfactors)
for (i in 1:length(Tfactors_pathway_genes$Gene_id)){
  car_id <- as.character(Tfactors_pathway_genes$Gene_id[i])
  f1_kallisto_df.Tfactors[car_id,10] <- as.character(Tfactors_pathway_genes$Enzyme[i])
}


MA_Tfactors <- ggplot(f1_kallisto_df.Tfactors, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("F1in vs F1out", subtitle =  "Trancription factors") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(box.padding= 1.7, size = 6, aes(label = ifelse(row.names(f1_kallisto_df.Tfactors) %in% Tfactors_pathway_genes$Gene_id, label, ""))) + 
  geom_point(data = subset(f1_kallisto_df.Tfactors, row.names(f1_kallisto_df.Tfactors) %in% Tfactors_pathway_genes$Gene_id), 
             aes(x = baseMean, y = adjLog2change), color = "black", 
             size = 1.8, shape = 21, stroke = 1)+ 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18), axis.text=element_text(size=12))



png(paste(plot.dir, "Tfactors_MA_plot_f1in_vs_f1out.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MA_Tfactors
dev.off()

#secondary pathway
f1_kallisto_df.secPathways <- f1_kallisto_df.MA

secPathways_pathway_genes <- manual_annotation[c(24:26,34:42, 75, 80,81), c(2,4)]
secPathways_pathway_genes$Enzyme <- c("HL3", "HL4", "As", "T19H", "TEX1", "TAT", "MAT", "HYS", "VAS_1", "VAS_2", "THAS", "V19H", "TEX2", "THAS2", "THAS3")

#make manual annotation collum
f1_kallisto_df.secPathways$label <- row.names(f1_kallisto_df.secPathways)
for (i in 1:length(secPathways_pathway_genes$Gene_id)){
  car_id <- as.character(secPathways_pathway_genes$Gene_id[i])
  f1_kallisto_df.secPathways[car_id,10] <- as.character(secPathways_pathway_genes$Enzyme[i])
}

down_label_genes <- secPathways_pathway_genes$Gene_id[c(1,2,10:13)]
up_label_genes <- secPathways_pathway_genes$Gene_id[-c(1,2,10:13)]


MA_secPathways <- ggplot(f1_kallisto_df.secPathways, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("F1in vs F1out", subtitle =  "Secondary Pathways enzymes") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(box.padding= 3, size = 6, aes(label = ifelse(row.names(f1_kallisto_df.secPathways) %in% up_label_genes, label, ""))) +
  geom_label_repel(xlim = c(1,4.5), box.padding= 2.6, size = 6, aes(label = ifelse(row.names(f1_kallisto_df.secPathways) %in% down_label_genes, label, ""))) + 
  geom_point(data = subset(f1_kallisto_df.secPathways, row.names(f1_kallisto_df.secPathways) %in% secPathways_pathway_genes$Gene_id), 
             aes(x = baseMean, y = adjLog2change), color = "black", 
             size = 1.8, shape = 21, stroke = 1)+ 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18), axis.text=element_text(size=12))



png(paste(plot.dir, "secPathways_MA_plot_f1in_vs_f1out.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MA_secPathways
dev.off()

#for f4in vs f4out ----
f4_kallisto <- read.csv("DESeq2_tables/f4in_vs_f4out_kallisto_Res.csv", row.names = 1)

#Build a dataframe with the significance of each gene (upregulated, downregulated, none)

f4_kallisto_df.MA <- f4_kallisto
f4_kallisto_df.MA$sig <- c()
f4_kallisto_df.MA$sig <- as.factor(ifelse(is.na(f4_kallisto_df.MA$padj), "Not tested",
                                      ifelse(f4_kallisto_df.MA$padj > 0.05, "Non Sig", 
                                      ifelse(f4_kallisto$log2FoldChange >= 1, "up", 
                                      ifelse(f4_kallisto_df.MA$log2FoldChange <= -1, "down", "Non Sig")))))

y_upper_limit <- 5
y_lower_limit <- -5

f4_kallisto_df.MA$shape <- ifelse(f4_kallisto_df.MA$log2FoldChange >= y_upper_limit, "Tup", ifelse(f4_kallisto_df.MA$log2FoldChange <= y_lower_limit, "Tdown", "normal"))
f4_kallisto_df.MA$adjLog2change <- ifelse(f4_kallisto_df.MA$log2FoldChange >= y_upper_limit, y_upper_limit, ifelse(f4_kallisto_df.MA$log2FoldChange <= y_lower_limit, y_lower_limit, f4_kallisto_df.MA$log2FoldChange))


MAplot <- ggplot(f4_kallisto_df.MA, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("f4in vs f4out") + xlab("Mean Expression") + ylab("Log2FoldChange") + 
  theme(plot.title = element_text(size=22, face = "bold"))


png(paste(plot.dir, "base_MA_plot_f4in_vs_f4out.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MAplot
dev.off()


#plot pathway specific MA plots 
manual_annotation <- read.xlsx(file = "StepsToFinalAnnotation/manual_annotation_modified.xlsx", sheetIndex = 1)


#iridoid Pathway
f4_kallisto_df.iridoid <- f4_kallisto_df.MA
iridoid_pathway_genes <- manual_annotation[1:6, c(2,4)]
iridoid_pathway_genes$Enzyme <- c("G10H", "10HGO", "IS", "IO", "7DLGT", "7DLH")


#make manual annotation collum
f4_kallisto_df.iridoid$label <- row.names(f4_kallisto_df.iridoid)
for (i in 1:length(iridoid_pathway_genes$Gene_id)){
  car_id <- as.character(iridoid_pathway_genes$Gene_id[i])
  f4_kallisto_df.iridoid[car_id,10] <- as.character(iridoid_pathway_genes$Enzyme[i])
}


MA_iridoid <- ggplot(f4_kallisto_df.iridoid, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("f4in vs f4out", subtitle = "Iridoid Pathway") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(box.padding= .5, size = 6, aes(label = ifelse(row.names(f4_kallisto_df.iridoid) %in% iridoid_pathway_genes$Gene_id, label, ""))) + 
  geom_point(data = subset(f4_kallisto_df.iridoid, row.names(f4_kallisto_df.iridoid) %in% iridoid_pathway_genes$Gene_id), 
             aes(x = baseMean, y = adjLog2change), color = "black", 
             size = 1.8, shape = 21, stroke = 1.5) + 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18), 
        axis.text=element_text(size=12))


png(paste(plot.dir, "iridoid_MA_plot_f4in_vs_f4out.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MA_iridoid
dev.off()


#central Pathway
#Plot the first part of the central pathway
f4_kallisto_df.central_1 <- f4_kallisto_df.MA
central_1_pathway_genes <- manual_annotation[c(7:11,27,77:79), c(2,4)]
central_1_pathway_genes$Enzyme <- c("LAMT","SLS2", "TDC", "STR", "SGD", "Asa", "SLS1", "SLS3", "SLS4")


#make manual annotation collum
f4_kallisto_df.central_1$label <- row.names(f4_kallisto_df.central_1)
for (i in 1:length(central_1_pathway_genes$Gene_id)){
  car_id <- as.character(central_1_pathway_genes$Gene_id[i])
  f4_kallisto_df.central_1[car_id,10] <- as.character(central_1_pathway_genes$Enzyme[i])
}

up_label_genes <- central_1_pathway_genes$Gene_id[c(2,5,8,9)]
down_label_genes <- central_1_pathway_genes$Gene_id[-c(2,5,8,9)]

MA_central_1 <- ggplot(f4_kallisto_df.central_1, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("f4in vs f4out", subtitle =  "central Pathway") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(box.padding= 2, size = 6, aes(label = ifelse(row.names(f4_kallisto_df.central_1) %in% central_1_pathway_genes$Gene_id, label, ""))) + 
  geom_point(data = subset(f4_kallisto_df.central_1, row.names(f4_kallisto_df.central_1) %in% central_1_pathway_genes$Gene_id), 
             aes(x = baseMean, y = adjLog2change), color = "black", 
             size = 1.8, shape = 21, stroke = 1.5)+ 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18), axis.text=element_text(size=12))



png(paste(plot.dir, "central_1_MA_plot_f4in_vs_f4out.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MA_central_1
dev.off()

#Plot the second part of the central pathway
f4_kallisto_df.central_2 <- f4_kallisto_df.MA
central_2_pathway_genes <- manual_annotation[c(19:23,28:32,73), c(2,4)]
central_2_pathway_genes$Enzyme <- c("ASO", "DPAS", "GS1","CS","TS", "GO_1","GO_2", "Redox 1", "Redox 2", "SAT", "GS2")


f4_kallisto_df.central_2$label <- row.names(f4_kallisto_df.central_2)
for (i in 1:length(central_2_pathway_genes$Gene_id)){
  car_id <- as.character(central_2_pathway_genes$Gene_id[i])
  f4_kallisto_df.central_2[car_id,10] <- as.character(central_2_pathway_genes$Enzyme[i])
}

central_2_pathway_genes_cut <- central_2_pathway_genes[-7,]

up_label_genes <- central_2_pathway_genes$Gene_id[c(2,3,6,9,11)]
down_label_genes <- central_2_pathway_genes$Gene_id[-c(2,3,6,9,11)]


MA_central_2 <- ggplot(f4_kallisto_df.central_2, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("f4in vs f4out", subtitle =  "central Pathway") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(ylim = c(0,3), box.padding= 1.3, size = 6, aes(label = ifelse(row.names(f4_kallisto_df.central_2) %in% up_label_genes, label, ""))) + 
  geom_label_repel(box.padding= 2.6, size = 6, aes(label = ifelse(row.names(f4_kallisto_df.central_2) %in% down_label_genes, label, ""))) + 
  geom_point(data = subset(f4_kallisto_df.central_2, row.names(f4_kallisto_df.central_2) %in% central_2_pathway_genes_cut$Gene_id), 
             aes(x = baseMean, y = adjLog2change), color = "black", 
             size = 1.8, shape = 21, stroke = 1.5)+ 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18), axis.text=element_text(size=12))



png(paste(plot.dir, "central_2_MA_plot_f4in_vs_f4out.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MA_central_2
dev.off()


#tamberstonine to vindoline pathway
f4_kallisto_df.vindoline <- f4_kallisto_df.MA

vindoline_pathway_genes <- manual_annotation[c(12:18,33, 76), c(2,4)]
vindoline_pathway_genes$Enzyme <- c("T16H2","16OMT", "T3O", "T3R", "16NMT", "D4H", "DAT", "PRX1", "T16H1")

#make manual annotation collum
f4_kallisto_df.vindoline$label <- row.names(f4_kallisto_df.vindoline)
for (i in 1:length(vindoline_pathway_genes$Gene_id)){
  car_id <- as.character(vindoline_pathway_genes$Gene_id[i])
  f4_kallisto_df.vindoline[car_id,10] <- as.character(vindoline_pathway_genes$Enzyme[i])
}



MA_vindoline <- ggplot(f4_kallisto_df.vindoline, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 2) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("f4in vs f4out", subtitle =  "Vindoline Pathway") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(xlim = c(1,NA), ylim = c(NA, 3), box.padding= 2.3, size = 6, aes(label = ifelse(row.names(f4_kallisto_df.vindoline) %in% vindoline_pathway_genes$Gene_id, label, ""))) + 
  geom_point(data = subset(f4_kallisto_df.vindoline, row.names(f4_kallisto_df.vindoline) %in% vindoline_pathway_genes$Gene_id), 
             aes(x = baseMean, y = adjLog2change), color = "black", 
             size = 2, shape = 21, stroke = 1.5)+ 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18), axis.text=element_text(size=12))



png(paste(plot.dir, "vindoline_MA_plot_f4in_vs_f4out.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MA_vindoline
dev.off()


#transporters 
f4_kallisto_df.transporters <- f4_kallisto_df.MA

transporters_pathway_genes <- manual_annotation[c(61:72), c(2,4)]
transporters_pathway_genes$Enzyme <- c("CrTPT2", "CrNPF2.9", "CrNPF2.4", "CrNPF2.5", "CrNPF2.6", "MATE1_KX372304", "MATE2_KX372305", "CrMate2PlantBio", "CrMate3PlantBio", "CrMate4PlantBio", "CrMate5PlantBio", "CrMate6PlantBio")

#make manual annotation collum
f4_kallisto_df.transporters$label <- row.names(f4_kallisto_df.transporters)
for (i in 1:length(transporters_pathway_genes$Gene_id)){
  car_id <- as.character(transporters_pathway_genes$Gene_id[i])
  f4_kallisto_df.transporters[car_id,10] <- as.character(transporters_pathway_genes$Enzyme[i])
}


MA_transporters <- ggplot(f4_kallisto_df.transporters, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 2) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("f4in vs f4out", subtitle =  "Transporters") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(box.padding= 0.2, size = 6, aes(label = ifelse(row.names(f4_kallisto_df.transporters) %in% transporters_pathway_genes$Gene_id, label, ""))) + 
  geom_point(data = subset(f4_kallisto_df.transporters, row.names(f4_kallisto_df.transporters) %in% transporters_pathway_genes$Gene_id), 
             aes(x = baseMean, y = adjLog2change), color = "black", 
             size = 2, shape = 21, stroke = 0.9)+ 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18),  axis.text=element_text(size=12))



png(paste(plot.dir, "transporters_MA_plot_f4in_vs_f4out.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MA_transporters
dev.off()

#TFs
f4_kallisto_df.Tfactors <- f4_kallisto_df.MA

Tfactors_pathway_genes <- manual_annotation[c(43:60), c(2,4)]
Tfactors_pathway_genes$Enzyme <- c( "Bis1", "Bis2", "ORCA2", "ORCA3", "ORCA4", "ORCA5", "ORCA6", "CR1", "WRKY1", "MYC2", "ZCT1", "ZCT2", "ZCT3", "GBF1", "GBF2", "BPF1", "CrGATA1", "ERF5")

#make manual annotation collum
f4_kallisto_df.Tfactors$label <- row.names(f4_kallisto_df.Tfactors)
for (i in 1:length(Tfactors_pathway_genes$Gene_id)){
  car_id <- as.character(Tfactors_pathway_genes$Gene_id[i])
  f4_kallisto_df.Tfactors[car_id,10] <- as.character(Tfactors_pathway_genes$Enzyme[i])
}

Tfactors_pathway_genes_cut <- Tfactors_pathway_genes$Gene_id[-4]


MA_Tfactors <- ggplot(f4_kallisto_df.Tfactors, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("f4in vs f4out", subtitle =  "Trancription factors") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(box.padding= 0.1, size = 6, aes(label = ifelse(row.names(f4_kallisto_df.Tfactors) %in% Tfactors_pathway_genes_cut, label, ""))) + 
  geom_label_repel(ylim = c(-4, -4.2), box.padding= 0.1, size = 6, aes(label = ifelse(row.names(f4_kallisto_df.Tfactors) %in% "CRO_110360", label, ""))) + 
  geom_point(data = subset(f4_kallisto_df.Tfactors, row.names(f4_kallisto_df.Tfactors) %in% Tfactors_pathway_genes$Gene_id), 
             aes(x = baseMean, y = adjLog2change), color = "black", 
             size = 1.8, shape = 21, stroke = 1)+ 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18), 
        axis.text=element_text(size=12))



png(paste(plot.dir, "Tfactors_MA_plot_f4in_vs_f4out.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MA_Tfactors
dev.off()

#secondary pathway
f4_kallisto_df.secPathways <- f4_kallisto_df.MA

secPathways_pathway_genes <- manual_annotation[c(24:26,34:42, 75, 80,81), c(2,4)]
secPathways_pathway_genes$Enzyme <- c("HL3", "HL4", "As", "T19H", "TEX1", "TAT", "MAT", "HYS", "VAS_1", "VAS_2", "THAS", "V19H", "TEX2", "THAS2", "THAS3")

#make manual annotation collum
f4_kallisto_df.secPathways$label <- row.names(f4_kallisto_df.secPathways)
for (i in 1:length(secPathways_pathway_genes$Gene_id)){
  car_id <- as.character(secPathways_pathway_genes$Gene_id[i])
  f4_kallisto_df.secPathways[car_id,10] <- as.character(secPathways_pathway_genes$Enzyme[i])
}

secPathways_pathway_genes_cut <- secPathways_pathway_genes$Gene_id[-1]
  
  
MA_secPathways <- ggplot(f4_kallisto_df.secPathways, aes(x = baseMean, y = adjLog2change, group = shape))+ 
  geom_point(aes(shape=shape, color = sig), show.legend = F, size = 1.8) + 
  scale_shape_manual(values=c(19, 6, 2)) + 
  scale_color_manual(values = c("green1", "gray30", "grey", "tomato3")) + 
  scale_x_continuous(trans = 'log10', 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_hline(yintercept = 0, color = "black", size = 0.3) + 
  theme_classic() + 
  ggtitle("f4in vs f4out", subtitle =  "Secondary Pathway enzymes") + xlab("Mean Expression") + ylab("Log2FoldChange") +
  geom_label_repel(box.padding= 1.7, size = 5, aes(label = ifelse(row.names(f4_kallisto_df.secPathways) %in% secPathways_pathway_genes$Gene_id, label, ""))) + 
  geom_point(data = subset(f4_kallisto_df.secPathways, row.names(f4_kallisto_df.secPathways) %in% secPathways_pathway_genes_cut), 
             aes(x = baseMean, y = adjLog2change), color = "black", 
             size = 1.8, shape = 21, stroke = 1)+ 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18), axis.text=element_text(size=12))



png(paste(plot.dir, "secPathways_MA_plot_f4in_vs_f4out.png", sep = "/"), width = 400, height = 250, units = "mm", res = 500)
MA_secPathways
dev.off()

