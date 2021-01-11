#Use this script to build MM and GS plots

#misc 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())


library(ggplot2)
library(xlsx)
library(ggrepel)


#for  module2, show the alkaloid realted genes (and label DPAS and CrMATE5)----
#the PUP transporters

module2Module <- read.xlsx(file = "results/genesPermodule.xlsx", sheetName = "Module 2")#change to salmon

alkaloid_genesMagentaModule <- read.xlsx(file = "results/topGO_enrichement/ IDs_salmon.xlsx", sheetIndex = 3)

module2Module$GeneSet <- ifelse(module2Module$NA. %in% alkaloid_genesMagentaModule$NA., "Alkaloid-related enzymes", "normal")

module2Module$GeneSet[38] <- "Alkaloid-related enzymes"
module2Module$GeneSet[1] <- "Transporters"
module2Module$GeneSet[4] <- "Transporters"
module2Module$GeneSet[47] <- "Transporters"



#manual_annotation label
module2Module$label <- NA
module2Module$label[123] <- "DPAS (Dehydroprecondylocarpine acetate synthase) *"
module2Module$label[9] <- "DPAS-like (CRO_133495)"
module2Module$label[35] <- "DPAS-like (CRO_116006)"
module2Module$label[63] <- "Methyltransferase (CRO_104106)"
module2Module$label[4] <- "CrMATE5"
module2Module$label[1] <- "PUP (CRO_111285)"
module2Module$label[47] <- "PUP (CRO_129702)"



#for vindoline 
GSvindolineMMplot <- ggplot(module2Module, aes(x = MM, y= GSvindoline)) + 
  geom_point(aes(color = GeneSet, shape = GeneSet, size = GeneSet)) + 
  scale_colour_manual(values = c("red", "blue","black", "green")) +
  scale_shape_manual(values =  c(15,16,1,17)) + 
  scale_size_manual(values = c(3,3,1,3)) + 
  theme_classic() + 
  ggtitle("Gene Significance for vindoline VS module Membership", subtitle = "module2") + 
  geom_label_repel(box.padding = 0.1,aes(label = ifelse(!is.na(label), label, ""))) + 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18),   axis.text=element_text(size=12)) + 
  ylab("Gene significance Vindoline") + xlab("Module Membership")
 

png("results/GSMMplots/GSvindolineMMplot_module2.png", res = 300, width = 2800, height = 2000)
GSvindolineMMplot
dev.off()  



#for catharanthine
GScatharanthineMMplot <- ggplot(module2Module, aes(x = MM, y= GSCatharanthine)) + 
  geom_point(aes(color = GeneSet, shape = GeneSet, size = GeneSet)) + 
  scale_colour_manual(values = c("red", "blue","black","green")) +
  scale_shape_manual(values =  c(15,16,1,17)) + 
  scale_size_manual(values = c(3,3,1,3)) + 
  theme_classic() + 
  ggtitle("Gene Significance for catharanthine VS module Membership", subtitle = "module2") + 
  geom_label_repel(aes(label = ifelse(!is.na(label), label, ""))) + 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18),   axis.text=element_text(size=12)) + 
  ylab("Gene significance catharanthine") + xlab("Module Membership")


png("results/GSMMplots/GScatharanthineMMplot_module2.png", res = 300, width = 2800, height = 2000)
GScatharanthineMMplot
dev.off()  


#for AVBL
GSAVBLMMplot <- ggplot(module2Module, aes(x = MM, y= GSAVBL)) + 
  geom_point(aes(color = GeneSet, shape = GeneSet, size = GeneSet)) + 
  scale_colour_manual(values = c("red","black", "green")) +
  scale_shape_manual(values =  c(15,1,17)) + 
  scale_size_manual(values = c(3,1,3)) + 
  theme_classic() + 
  ggtitle("Gene Significance for AVLB VS module Membership", subtitle = "module2") + 
  geom_label_repel(aes(label = ifelse(!is.na(label), label, ""))) + 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18),   axis.text=element_text(size=12)) + 
  ylab("Gene significance AVLB") + xlab("Module Membership") + 
  xlim(0.7,1) + ylim(0.7,1)


png("results/GSMMplots/GSAVBLMMplot_module2.png", res = 300, width = 2800, height = 2000)
GSAVBLMMplot
dev.off()  


#for vinblastine
GSvinblastineMMplot <- ggplot(module2Module, aes(x = MM, y= GSvinblastine)) + 
  geom_point(aes(color = GeneSet, shape = GeneSet, size = GeneSet)) + 
  scale_colour_manual(values = c("red", "blue","black", "green")) +
  scale_shape_manual(values =  c(15,16,1,17)) + 
  scale_size_manual(values = c(3,3,1,3)) + 
  theme_classic() + 
  ggtitle("Gene Significance for vinblastine VS module Membership", subtitle = "module2") + 
  geom_label_repel(box.padding= 0.6, aes(label = ifelse(!is.na(label), label, ""))) + 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18),   axis.text=element_text(size=12)) + 
  ylab("Gene significance vinblastine") + xlab("Module Membership")


png("results/GSMMplots/GSvinblastineMMplot_module2.png", res = 300, width = 2800, height = 2000)
GSvinblastineMMplot
dev.off()  

#for darkgreenmodule, show alkaloids metabolic process (the other alkaloid terms are correlated?)----

darkGreenModule <- read.xlsx(file = "results/genesPermodule.xlsx", sheetName = "darkgreen")

alkaloid_genesdarkGreenModule <- read.xlsx(file = "results/topGO_enrichement/ IDs_darkgreen.xlsx", sheetIndex = 1)

darkGreenModule$GeneSet <- ifelse(darkGreenModule$NA. %in% alkaloid_genesdarkGreen$NA., "alkaloid", "normal")

#manual_annotation label
darkGreenModule$label <- NA
darkGreenModule$label[48] <- "16OMT"
darkGreenModule$label[49] <- "T16H2"
darkGreenModule$label[65] <- "THAS"
darkGreenModule$label[94] <- "DAT"
darkGreenModule$label[95] <- "TAT-like (CRO_120038)"
darkGreenModule$label[104] <- "T19H-like (CRO_123713)"
darkGreenModule$label[122] <- "D4H"
darkGreenModule$label[123] <- "D4H-like (CRO_127167)"
darkGreenModule$label[182] <- "Vinorine synthase-like (CRO_137612)"



#for vindoline 
GSvindolineMMplot <- ggplot(darkGreenModule, aes(x = MM, y= GSvindoline)) + 
  geom_point(aes(color = GeneSet, shape = GeneSet, size = GeneSet)) + 
  scale_colour_manual(values = c("red","black")) +
  scale_shape_manual(values =  c(15,1)) + 
  scale_size_manual(values = c(3,1)) + 
  theme_classic() + 
  ggtitle("Gene Significance for vindoline VS module Membership", subtitle = "DarkGreen") + 
  geom_label_repel(box.padding = 0.5,aes(label = ifelse(darkGreenModule$`NA.` %in% alkaloid_genesMagentaModule$`NA.`, label, ""))) + 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18),   axis.text=element_text(size=12)) + 
  ylab("Gene significance Vindoline") + xlab("Module Membership")


png("results/GSMMplots/GSvindolineMMplot_DarkGreen.png", res = 300, width = 2800, height = 2000)
GSvindolineMMplot
dev.off()  



#for catharanthine
GScatharanthineMMplot <- ggplot(darkGreenModule, aes(x = MM, y= GSCatharanthine)) + 
  geom_point(aes(color = GeneSet, shape = GeneSet, size = GeneSet)) + 
  scale_colour_manual(values = c("red","black")) +
  scale_shape_manual(values =  c(15,1)) + 
  scale_size_manual(values = c(3,1)) + 
  theme_classic() + 
  ggtitle("Gene Significance for catharanthine VS module Membership", subtitle = "DarkGreen") + 
  geom_label_repel(aes(label = ifelse(darkGreenModule$`NA.` %in% alkaloid_genesMagentaModule$`NA.`, label, ""))) + 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18),   axis.text=element_text(size=12)) + 
  ylab("Gene significance catharanthine") + xlab("Module Membership")


png("results/GSMMplots/GScatharanthineMMplot_DarkGreen.png", res = 300, width = 2800, height = 2000)
GScatharanthineMMplot
dev.off()  


#for AVBL
GSAVBLMMplot <- ggplot(darkGreenModule, aes(x = MM, y= GSAVBL)) + 
  geom_point(aes(color = GeneSet, shape = GeneSet, size = GeneSet)) + 
  scale_colour_manual(values = c("red","black")) +
  scale_shape_manual(values =  c(15,1)) + 
  scale_size_manual(values = c(3,1)) + 
  theme_classic() + 
  ggtitle("Gene Significance for AVBL VS module Membership", subtitle = "DarkGreen") + 
  geom_label_repel(box.padding= 0.5, aes(label = ifelse(darkGreenModule$`NA.` %in% alkaloid_genesMagentaModule$`NA.`, label, ""))) + 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18),   axis.text=element_text(size=12)) + 
  ylab("Gene significance AVBL") + xlab("Module Membership")


png("results/GSMMplots/GSAVBLMMplot_DarkGreen.png", res = 300, width = 2800, height = 2000)
GSAVBLMMplot
dev.off()  



#for Tan module, show terpenoid catabolic process + NTR4.6----

tan <- read.xlsx(file = "results/genesPermodule.xlsx", sheetName = "tan")

catabolicProcess_genesTanModule <- read.xlsx(file = "results/topGO_enrichement/ IDs_tan.xlsx", sheetIndex = 4)

tan$GeneSet <- ifelse(tan$NA. %in% catabolicProcess_genesTanModule$NA., "Terpenoid \ncatabolic \nProcess", "normal")
tan$GeneSet[165] <- "Transporter"

#manual_annotation label
tan$label <- NA
tan$label[165] <- "NPF4.6-homolog"
tan$label[68] <- "ABAH-homolog"
tan$label[169] <- "G2O-homolog"
tan$label[209] <- "CCD1-homolog"


#for vinblastine
GSvinblastineMMplot <- ggplot(tan, aes(x = MM, y= GSvinblastine)) + 
  geom_point(aes(color = GeneSet, shape = GeneSet, size = GeneSet)) + 
  scale_colour_manual(values = c("black","red", "green")) +
  scale_shape_manual(values =  c(1,15,17)) + 
  scale_size_manual(values = c(1,3,3)) + 
  theme_classic() + 
  ggtitle("Gene Significance for vinblastine VS module Membership", subtitle = "Tan") + 
  geom_label_repel(box.padding = 1, aes(label = ifelse(is.na(label), "", label))) + 
  theme(plot.title = element_text(size=22, face = "bold"), 
        plot.subtitle = element_text(size=18),   axis.text=element_text(size=12)) + 
  ylab("Gene significance vinblastine") + xlab("Module Membership")


png("results/GSMMplots/GSvinblastineMMplot_Tan.png", res = 300, width = 2800, height = 2000)
GSvinblastineMMplot
dev.off()  

