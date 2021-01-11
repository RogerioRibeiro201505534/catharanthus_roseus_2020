# this script reads .txt files with significantly DEG between different conditions and 
#draws Venn diagrams to find common genes
setwd("C:/Users/Faculdade/Desktop/Dissertação/results 2020/mapping_to_genome/07_DESeq2/featureCounts")

rm(list = ls())

library('VennDiagram')
library('gridExtra')
library('grid')
library('ggplot2')

### function ###

draw.venn.2samples <- function(names1, names2, conditions.venn, colors.venn){
  area1 <- length(names1)
  area2 <- length(names2)
  names12 <- names1[names1 %in% names2]
  area12 <- length(names12)
  
  grid.newpage()
  temp <- draw.pairwise.venn(area1 = area1,
                             area2 = area2, 
                             cross.area = area12, 
                             category = conditions.venn, 
                             lty = rep("blank", 2),
                             fill = colors.venn, 
                             alpha = rep(0.5, 2), 
                             cat.pos = c(202,158), 
                             cat.cex = 1.5,
                             cex = 1.5,
                             ind = FALSE)
  grid.draw(temp)
  
  return(names12)
}

draw.venn.3samples <- function(names1, names2, names3, conditions.venn, colors.venn){
  area1 <- length(names1)
  area2 <- length(names2)
  area3 <- length(names3)
  names12 <- names1[names1 %in% names2]
  area12 <- length(names12)
  names13 <- names1[names1 %in% names3]
  area13 <- length(names13)
  names23 <- names2[names2 %in% names3] 
  area23 <- length(names23)
  names123 <- names12[names12 %in% names3]
  area123 <- length(names123)
  
  grid.newpage()
  temp <- draw.triple.venn(area1 = area1,
                           area2 = area2,
                           area3 = area3,
                           n12 = area12,
                           n23 = area23,
                           n13 = area13,
                           n123 = area123,
                           category = conditions.venn, 
                           lty = rep("blank", 3), 
                           fill = colors.venn, 
                           alpha = rep(0.5, 3), 
                           cat.dist = rep(0.005, 3),
                           cex = 0.8,
                           ind = FALSE)
  grid.draw(temp)
  
  return(names123)
}


out.dir <- ("./Venn")
if (!(file.exists(out.dir))) {
  dir.create(out.dir)
}


#Upregulated genes
idio_up <- read.csv(file = "meso_vs_idio/results/_sigUP.csv", row.names = 1)
f1out_up <- read.csv(file = "F1in_vs_F1out/results/_sigUP.csv", row.names =  1)
f4out_up <- read.csv(file = "F4in_vs_F4out/results/_sigUP.csv", row.names =  1)

#####Get Venn Diagrams and lists of acession#####


#Between all 3 DeSeq2 results  
conditions <- c("meso_vs_idio_UP", "f1in_vs_f1out_Up", "f4in_vs_f4out_Up")
colours.Venn <- c("green", "blue", "red")

png(paste(out.dir, "Venn_Geral_UP.png", sep="/"), width=12,height=8, units='in', res=300)
Venn_geral <- draw.venn.3samples(idio_up$names.sig.UP, f1out_up$names.sig.UP, f4out_up$names.sig.UP, conditions, colours.Venn)
dev.off()


write.csv(Venn_geral, file = paste(out.dir, "Venn_geral_UP_list.csv", sep = "/"))


###Between idioUp and F1Up DeSeq2 results###

idio_F1_up <- intersect(idio_up$names.sig.UP, f1out_up$names.sig.UP)
write.csv(idio_F1_up, file = paste(out.dir, "Venn_Idio_F1Out_UP_list.csv", sep = "/"))

###Make diagram
#conditions <- c("idioUp", "f1out_Up")
#colours.Venn <- c("green", "blue")
#png(paste(out.dir, "Venn_Idio_F1Out_UP.png", sep="/"), width=12,height=8, units='in', res=300)
#Venn_Idio_F1Out <- draw.venn.2samples(idio_up$names.sig.UP, f1out_up$names.sig.UP, conditions, colours.Venn)
#dev.off()
#write.csv(Venn_Idio_F1Out, file = paste(out.dir, "Venn_Idio_F1Out_UP_list.csv", sep = "/"))



###Between idioUp and F4Up DeSeq2 results###

idio_F4_up <- intersect(idio_up$names.sig.UP, f4out_up$names.sig.UP)
write.csv(idio_F4_up, file = paste(out.dir, "Venn_Idio_F4Out_UP_list.csv", sep = "/"))

###Make diagram
#conditions <- c("idioUp", "f4out_Up")
#colours.Venn <- c("green", "red")
#png(paste(out.dir, "Venn_Idio_F4Out_UP.png", sep="/"), width=12,height=8, units='in', res=300)
#Venn_Idio_F4Out <- draw.venn.2samples(idio_up$names.sig.UP, f4out_up$names.sig.UP, conditions, colours.Venn)
#dev.off()
#write.csv(Venn_Idio_F4Out, file = paste(out.dir, "Venn_Idio_F4Out_UP_list.csv", sep = "/"))


###Between F1up and F4Up DeSeq2 results###
F1_F4_up <- intersect(f4out_up$names.sig.UP, f1out_up$names.sig.UP)
write.csv(F1_F4_up, file = paste(out.dir, "Venn_F1out_F4Out_UP_list.csv", sep = "/"))

###Make diagram
#conditions <- c("f1out_up", "f4out_Up")
#colours.Venn <- c("blue", "red")
#png(paste(out.dir, "Venn_F1OUT_F4Out_UP.png", sep="/"), width=12,height=8, units='in', res=300)
#Venn_F1OUT_F4Out <- draw.venn.2samples(f1out_up$names.sig.UP, f4out_up$names.sig.UP, conditions, colours.Venn)
#dev.off()
#write.csv(Venn_F1OUT_F4Out, file = paste(out.dir, "Venn_F1Out_F4Out_UP_list.csv", sep = "/"))
