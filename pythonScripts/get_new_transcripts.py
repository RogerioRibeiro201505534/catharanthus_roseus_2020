#use this script to get the new Genes from the StringTie transcriptome 
#New genes here are defined as genes id that do not have a transcript that appeard in the original
#Genome assembly V2 (CRO_T******). 
#during the StringTie running, genes with new transcripts, supported by the transcriptomic data are given 
#new id (in this case CATHA**** ). Therefore, new genes are CATHA*** genes_id without CRO_T*** transcripts 

import csv 
from Bio import SeqIO

input_file = "transcript2gene_id.csv"
fasta_file = "merged_transcriptome.fasta"
 

def make_map(input_file):  
    #This function makes a gene_2_transcript dictionary. The file has to be tsv and be in this order
    with open(input_file, "r") as fh: 
        map = {}
        text = csv.reader(fh, delimiter = "\t")
        for line in text:
            if line[1] in map.keys(): 
                map[line[1]].append(line[0])
            else: 
                map[line[1]] = [line[0]]
    return map

def get_new_genes(gene2transcript_map): 
    new_genes = []
    for key in gene2transcript_map.keys(): 
        if is_new(gene2transcript_map[key]): 
            new_genes.append(key)
    return new_genes
            

def is_new(transcript_list): 
    for transcript in transcript_list: 
        if "CRO_T" in transcript: 
            return False
    return True

def get_new_sequences(sequence_list, fasta_file):
    new_sequences_list = [] 
    merged_transcriptome = SeqIO.parse(fasta_file, "fasta")
    for transcript in merged_transcriptome: 
        if "CATHA" in transcript.id:
            gene_id = "CATHA.{}".format(transcript.id.split(".")[1])
        else: 
            gene_id = transcript.id
        
        if gene_id in sequence_list:
            new_sequences_list.append(transcript)

    SeqIO.write(new_sequences_list, "new_genes_transcripts.fasta", "fasta")
        

                
gene2transcrip_map = make_map(input_file)
new_gene_list = get_new_genes(gene2transcrip_map)

with open("new_sequences_id.txt", "a") as fh: 
    for gene in new_gene_list: 
        fh.write(gene)
        fh.write("\n")

get_new_sequences(new_gene_list, fasta_file)
