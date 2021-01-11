##Make geral plots 
setwd("C:/Users/Faculdade/Desktop/Dissertação/results 2020/nanopore")
rm(list = ls())

#Load libraries
library(xlsx)
library(ggplot2)
library(DESeq2)
library(dplyr)


#create folder 
out.folder <- "plots"
if (!dir.exists(out.folder)){
  dir.create(out.folder)
}


#Load number of reads in the raw reads 
read_numbers <- read.delim(file = "01_process_reads/fasqc_raw/multiqc_data/multiqc_general_stats.txt")[,c(1,6)]
colnames(read_numbers)[2] <- "Raw"

#Load number of reads after pychopper
read_number_porechop <- read.delim(file = "01_process_reads/pychopper_oriented_all_counts.txt", sep = "\t", header = F)
read_numbers$Porechop <- read_number_porechop$V1

#Load the number of reads after trimming 
read_numbers$trimmed <- read.delim(file = "01_process_reads/trimmed_annotation_fastqc/multiqc_data/multiqc_general_stats.txt")[,6]

#Load the number of reads aligned
read_numbers$aligned <- read.csv(file = "05_quantification/salmon_read_counts.out", sep = "\t")[,2]



#Load Rin 
sample_rin <- c("4.4", "0", "6.3", "5.7", "6.3", "7.2", "6.3", "6.6", "6.8", "6.7", "7.1", "6")
read_numbers$Rin <- sample_rin

read_numbers <- read_numbers[,c(1,5,2:4)]

read_numbers$Sample <- factor(read_numbers$Sample, 
                       levels = c("P1F2Int", "P2F2Int", "P3F2Int",
                                  "P1F3Int", "P2F3Int", "P3F3Int", 
                                  "P1F2Ext", "P2F2Ext", "P3F2Ext",
                                  "P1F3Ext", "P2F3Ext", "P3F3Ext"))

#add a condition for adding shapes in the ggplot
read_numbers$condition <- rep(c("F2", "F2", "F3", "F3"), 3)

png(paste(out.folder, "Rin_vs_read_number.png", sep="/"), res = 300, width = 2000, height = 1500)

ggplot(read_numbers, aes(x = Rin, y = Raw, color = Sample)) + 
  geom_point(size = 4, aes(shape = condition)) + theme_classic() + 
  ggtitle("Rin vs Read number") + 
  ylab("Read number") +  
  scale_color_brewer(palette="Set3") + 
  scale_shape_manual(values=c(15,17))
dev.off()

#Make a plot of each read number trough each analysis step---- 
read_numbers_plot.df <- read_numbers

aligned_reads <- read_numbers[,c(1,5)]
colnames(aligned_reads)[2] <- "Reads"
aligned_reads$condition <- "Aligned"

trimmed_reads_offset <- read_numbers_plot.df[,c(1,4)]
trimmed_reads_offset[,2] <- read_numbers_plot.df$trimmed - read_numbers_plot.df$aligned
colnames(trimmed_reads_offset)[2] <- "Reads"
trimmed_reads_offset$condition <- "Trimmed"

read_raw_offset <- read_numbers_plot.df[,c(1,2)]
read_raw_offset[,2] <- read_numbers_plot.df$Raw - read_numbers_plot.df$trimmed
colnames(read_raw_offset)[2] <- "Reads"
read_raw_offset$condition <- "Raw"


read_numbers_plot.df <- rbind(rbind(aligned_reads,trimmed_reads_offset),read_raw_offset)

read_numbers_plot.df$condition <- factor(read_numbers_plot.df$condition,
                      levels = c("Raw", "Trimmed", "Aligned"))

read_numbers_plot.df$Sample <- factor(read_numbers_plot.df$Sample, 
                                      levels = c("P1F2Int", "P2F2Int", "P3F2Int",
                                                 "P1F3Int", "P2F3Int", "P3F3Int", 
                                                 "P1F2Ext", "P2F2Ext", "P3F2Ext",
                                                 "P1F3Ext", "P2F3Ext", "P3F3Ext"))

png(paste(out.folder, "Nanopore_read_number.png", sep = "/"),res = 300, height = 1500, width = 2500)

ggplot(read_numbers_plot.df, aes(fill = condition, y = Reads, x = Sample)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = c("#56B4E9","#F0E442", "#F04442")) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 290, hjust = 1, vjust = 0.5)) + 
  ggtitle("Nanopore read number") + ylab("Reads")

dev.off()


png(paste(out.folder, "Nanopore_read_number_perct.png", sep = "/"),res = 300, height = 1500, width = 2500)

ggplot(read_numbers_plot.df, aes(fill = condition, y = Reads, x = Sample)) + 
  geom_bar(stat="identity", position = "fill") + 
  scale_fill_manual(values = c("#56B4E9","#F0E442", "#F04442")) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 290, hjust = 1, vjust = 0.5)) + 
  ggtitle("Nanopore Read") + ylab("Reads percentage") + xlab("")
dev.off()


#save data tables in xlsx
write.xlsx(read_numbers, file = "Nanopore_reads.xlsx", sheetName = "Geral", row.names = F)



# salmon.illumina.ref$sample_2 <- read_numbers[,1]
# salmon.illumina.ref$Unmapped <- read_numbers[,5]
# 
# salmon.illumina.ref.Unmapped <- data.frame("sample_2" = salmon.illumina.ref$sample_2, "reads" = salmon.illumina.ref$Unmapped - salmon.illumina.ref$Aligned)
# salmon.illumina.ref.Unmapped$condition <- "Unmapped"
# 
# salmon.illumina.ref.ambigous <- salmon.illumina.ref[,c(5,4)] 
# salmon.illumina.ref.ambigous$condition <- "ambigous"
# colnames(salmon.illumina.ref.ambigous)[2] <- "reads"
# 
# salmon.illumina.ref.unique <- salmon.illumina.ref[,c(5,3)] 
# salmon.illumina.ref.unique$condition <- "unique"
# colnames(salmon.illumina.ref.unique)[2] <- "reads"
# 
# salmon.illumina.ref.plot <- rbind(rbind(salmon.illumina.ref.Unmapped, salmon.illumina.ref.ambigous), salmon.illumina.ref.unique)
# 
# salmon.illumina.ref.plot$condition <- factor(salmon.illumina.ref.plot$condition, levels = c("ambigous", "unique","Unmapped"))
# 
# 
# png(paste(out.folder, "salmon_IlluminaReference.png", sep = "/"), width = 900)
# ggplot(salmon.illumina.ref.plot, aes(fill = condition, y = reads, x = sample_2)) + 
#   geom_bar(stat="identity") + 
#   theme_classic() +
#   ggtitle("Reads alligned to transcriptome (with Illumina annotation)") + 
#   ylab("Salmon quantified reads") + 
#   scale_fill_manual(values = c("#56B4E9","#F0E442", "#F04442")) + 
#   theme(axis.text.x = element_text(angle = 290, hjust = 1, vjust = 0.5))
# dev.off()
# 
# png(paste(out.folder, "salmon_IlluminaReferencePercent.png", sep = "/"), width = 900)
# ggplot(salmon.illumina.ref.plot, aes(fill = condition, y = reads, x = sample_2)) + 
#   geom_bar(stat="identity", position = "fill") + 
#   theme_classic() +
#   ggtitle("Reads alligned to transcriptome (with Illumina annotation)") + 
#   ylab("Salmon quantified reads") + 
#   scale_fill_manual(values = c("#56B4E9","#F0E442", "#F04442")) + 
#   theme(axis.text.x = element_text(angle = 290, hjust = 1, vjust = 0.5))
# dev.off()
# 
# ###with new reference
# 
# salmon.hybrid.ref <- read.csv(file = "04_DE_pipeline/salmon_parsed.tsv", sep = "\t")
# salmon.hybrid.ref$sample_2 <- read_numbers[,1]
# salmon.hybrid.ref$Unmapped <- read_numbers[,5]
# 
# salmon.hybrid.ref.Unmapped <- data.frame("sample_2" = salmon.hybrid.ref$sample_2, "reads" = salmon.hybrid.ref$Unmapped - salmon.hybrid.ref$Aligned)
# salmon.hybrid.ref.Unmapped$condition <- "Unmapped"
# 
# salmon.hybrid.ref.ambigous <- salmon.hybrid.ref[,c(5,4)] 
# salmon.hybrid.ref.ambigous$condition <- "ambigous"
# colnames(salmon.hybrid.ref.ambigous)[2] <- "reads"
# 
# salmon.hybrid.ref.unique <- salmon.hybrid.ref[,c(5,3)] 
# salmon.hybrid.ref.unique$condition <- "unique"
# colnames(salmon.hybrid.ref.unique)[2] <- "reads"
# 
# salmon.hybrid.ref.plot <- rbind(rbind(salmon.hybrid.ref.Unmapped, salmon.hybrid.ref.ambigous), salmon.hybrid.ref.unique)
# 
# salmon.hybrid.ref.plot$condition <- factor(salmon.hybrid.ref.plot$condition, levels = c("ambigous", "unique","Unmapped"))
# 
# 
# png(paste(out.folder, "salmon_hybridReference.png", sep = "/"), width = 900)
# ggplot(salmon.hybrid.ref.plot, aes(fill = condition, y = reads, x = sample_2)) + 
#   geom_bar(stat="identity") + 
#   theme_classic() +
#   ggtitle("Reads alligned to transcriptome (with hybrid annotation)") + 
#   ylab("Salmon quantified reads") + 
#   scale_fill_manual(values = c("#56B4E9","#F0E442", "#F04442")) + 
#   theme(axis.text.x = element_text(angle = 290, hjust = 1, vjust = 0.5))
# dev.off()
# 
# png(paste(out.folder, "salmon_hybridReferencePercent.png", sep = "/"), width = 900)
# ggplot(salmon.hybrid.ref.plot, aes(fill = condition, y = reads, x = sample_2)) + 
#   geom_bar(stat="identity", position = "fill") + 
#   theme_classic() +
#   ggtitle("Reads alligned to transcriptome (with hybrid annotation)") + 
#   ylab("Salmon quantified reads") + 
#   scale_fill_manual(values = c("#56B4E9","#F0E442", "#F04442")) + 
#   theme(axis.text.x = element_text(angle = 290, hjust = 1, vjust = 0.5))
# dev.off()
