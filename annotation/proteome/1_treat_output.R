#Set workng directory 
setwd("C:/Users/Faculdade/Desktop/Dissertação/results 2020/mapping_to_genome/annotation/02_proteome_annotation")
rm(list = ls())

#Library 
library(seqinr)
library(xlsx)


##Load caros2filtered_transcriptome dataset
#load the minimap2 file (minimap2 -c -splice merged_transcriptome.fasta --secondary=no Catro_mRNA_update.tfa_annotated_ORFs_longest_ORF.fasta > minimap_output.paf)
minimap_alignements <- read.table("minimap_output.paf")
colnames(minimap_alignements) <- c("Query", "Query Len", "Query start", "Query end", "Strand", "target", "Target len", "Target start", "Target end", "Nº matches", "align size", "mapping quality")


###Build a dataset ---- 
#FOr each row, calculate the identity and the percentage of the query sequence covered by the alin

minimap_alignements.metrics <- minimap_alignements
minimap_alignements.metrics$Identity <- minimap_alignements.metrics$`Nº matches` / minimap_alignements.metrics$`align size`
minimap_alignements.metrics$query_overlap <- minimap_alignements.metrics$`align size` / minimap_alignements.metrics$`Query Len`
minimap_alignements.metrics <- minimap_alignements.metrics[,c(1,2,5,6,7,10,11,24,25)]


###Filter the dataset, relating to indentity and query overlap metrics ----
indentity_treshold <- 0.95
qoverlap_treshold <- 0.90

#Filter for ind
minimap_alignements.filtered <- minimap_alignements.metrics[which(minimap_alignements.metrics$Identity >= indentity_treshold), ]
minimap_alignements.filtered <- minimap_alignements.filtered[which(minimap_alignements.filtered$query_overlap >= qoverlap_treshold), ]

#2635 matches (out of 2736 caros sequences) 

#Test if there are duplicates
length(unique(minimap_alignements.filtered$Query))
length(unique(minimap_alignements.filtered$target))
#there are 76 duplicates in the target transcripts



####Add to the dataset annotation of the samples where the proteins are present 
fasta_file <- read.fasta(file = "cathacyc queries/Catro_mRNA_update.tfa_annotated_ORFs_longest_ORF_v2.fasta")
annot_vector <- c()
i <- 0 
for (query_name in minimap_alignements.filtered$Query){
  i <- i + 1
  new_col_annot <- unlist(strsplit(query_name, "_"))[1]
  description <- attributes(fasta_file[[query_name]])$Annot
  if (grepl("VAC_1", description, fixed = TRUE)){
    new_col_annot <- paste(new_col_annot, " Vac_1")
  }
  if (grepl("VAC_2", description, fixed = TRUE)){
    new_col_annot <- paste(new_col_annot, " Vac_2")
  }
  if (grepl("Ton_1", description, fixed = TRUE)){
    new_col_annot <- paste(new_col_annot, " Ton_1")
  }
  if (grepl("Ton_2", description, fixed = TRUE)){
    new_col_annot <- paste(new_col_annot, " Ton_2")
  }
  annot_vector[i] <- new_col_annot
}

minimap.final <- minimap_alignements.filtered
minimap.final$annot <- annot_vector


##Load dammit annotation and the minimap association table
dammit_annotation <- read.csv(file = "dammit_annotation_with_go.tsv", sep = "\t")


#First, we have to convert transcripts ids to gene id 
#load the transcript to gene_id table and convert in to a named vector 

transcript2gene <- read.csv(file = "transcript2gene.tsv", sep = "\t", header =  FALSE)
gene_vector <- transcript2gene[,2]
names(gene_vector) <- transcript2gene[,1]

minimap_gene_vector <- c()
for (i in 1:length(minimap.final$target)){
  gene <- as.character(gene_vector[as.character(minimap.final$target[i])])
  minimap_gene_vector[i] <- gene
}

minimap.final$gene_id <- minimap_gene_vector
minimap.final.abr <- minimap.final[,c(1,4,10,11)]

#remove the _x from the caros acession number
new_caros_vector <- c()
for (i in 1:length(minimap.final$target)){
  new_caros_vector[i] <- unlist(strsplit(as.character(minimap.final.abr$Query[i]), "_"))[1]
}

minimap.final.abr$Query <- new_caros_vector

#Add the cathacyc annotation to this table
cathacyc_functional_annotation <- read.csv(file = "proteome_papper/orcae_catro.function.2020_06.csv", sep =  ";", header = F)
cathacyc_annot_vector <- cathacyc_functional_annotation[,5]
names(cathacyc_annot_vector) <- cathacyc_functional_annotation[,2]


annot_vector <- c()
for (i in 1:length(minimap.final.abr$Query)){
  gene_name <- as.character(cathacyc_annot_vector[as.character(minimap.final.abr$Query[i])])
  annot_vector[i] <- gene_name
}

minimap.final.abr$cathacyc <- annot_vector


#Finally pass the annotation to the dammit_annot_table
minimap.final_proteome_annotation <- minimap.final.abr$annot
minimap.final_cathacyc_annotation <-minimap.final.abr$cathacyc
names(minimap.final_proteome_annotation) <- minimap.final.abr$gene_id
names(minimap.final_cathacyc_annotation) <- minimap.final.abr$gene_id

proteome_annotation_vector <- rep(NA, length(dammit_annotation$gene.id))
cathacyc_annotation_vector <- rep(NA, length(dammit_annotation$gene.id))

for (i in 1:length(dammit_annotation$gene.id)){
  rows <- which(minimap.final.abr$gene_id == dammit_annotation$gene.id[i])
  if (length(rows) != 0){
    for (j in 1:length(rows)){
      if (j == 1){
        proteome_location <- minimap.final_proteome_annotation[rows[j]]
        cathacyc_functional_annotation <- minimap.final_cathacyc_annotation[rows[j]]
      }
      else {
        proteome_location <- paste(proteome_location, minimap.final_proteome_annotation[rows[j]], sep = "; ")
        cathacyc_functional_annotation <- paste(cathacyc_functional_annotation, minimap.final_cathacyc_annotation[rows[j]], sep = "; ")
      }
    }
    proteome_annotation_vector[i] <- proteome_location
    cathacyc_annotation_vector[i] <- cathacyc_functional_annotation
  }
}


dammit_annotation$proteome_location <- proteome_annotation_vector
dammit_annotation$cathacyc_annotation <- cathacyc_annotation_vector


length(unique((minimap.final.abr$Query)))
length(unique((minimap.final.abr$gene_id))) #2467
#168 duplicates 

write.table(minimap.final.abr, "minimap_final.tsv", sep = "\t", row.names = F, quote = F)
write.xlsx(minimap.final.abr, "minimap_final.xlsx", row.names = F)
write.table(dammit_annotation, "dammit_annotation_with_proteome.tsv", sep = "\t", quote = F, row.names = F)

