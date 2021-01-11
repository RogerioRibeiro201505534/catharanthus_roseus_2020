#!/bin/bash

#Use this script to map the caros transcript to the genome-guided transcriptome
#Note: the excel file with the caros sequences with the caros sequences may not be provided ... 


#extracts the caros with annotation from the supp file 
py extract_caros.py
#Get the ORF's from each caros 
getorf -sequence Catro_mRNA_update.tfa_annotated -outseq Catro_mRNA_update.tfa_annotated_ORF.txt -find 3 -minsize 100
#get the longest ORF in each Caros
py extract_nucleotidic_ORF.py 

minimap2 -c -splice merged_transcriptome_filtered.fasta --secondary=no Catro_mRNA_update.tfa_annotated_ORFs_longest_ORF.fasta > minimap_output.paf

#Use this output with the R script provided 
