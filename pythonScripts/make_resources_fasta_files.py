#use this script to filter and annotate with the acession number the peptides and CDS fastas files 
#Edited 10/08/20


##Load libraries 
import csv 
import re
from Bio import SeqIO

#Load files
transcripts2gene = "transcript2gene.tsv"
name_map = "merged_transcriptome_modified.fasta.dammit.namemap.csv"
annotation_gff = "merged_transcriptome_modified.fasta.dammit.gff3"
cds_fasta_file = "merged_transcriptome_modified.fasta.transdecoder.cds"
pep_fasta_file = "merged_transcriptome_modified.fasta.transdecoder.pep"



def get_transcript_list(transcript2gene):
    #get the transcript_list to keep
    #in this case processed by reading the transcript2gene.tsv file 
    t_list = []
    with open(transcript2gene, "r") as fh: 
        t2g = csv.reader(fh, delimiter = "\t")
        for row in t2g: 
            t_list.append(row[0])
    return t_list

def get_name_map(name_map): 
    #get the association between dammit names and transcripts 
    dammit_id2transcript_map = {}
    with open(name_map, "r") as fh: 
        table = csv.reader(fh, delimiter = ",")
        next(table)
        for line in table: 
            dammit_id2transcript_map[line[1]] = line[0]
    return dammit_id2transcript_map


def parse_fasta_files(fasta_file, seq_type, t_list, dammit_id2transcript_map): 
    sequence_list = []
    with open(fasta_file, "r") as fh: 
        for record in SeqIO.parse(fh, "fasta"):    
            record_id = record.id.split(".")
            t_id = dammit_id2transcript_map[record_id[0]]
            if t_id in t_list: 
                record.id = "{}.{}".format(t_id, record_id[1])
                record.description = ""
                sequence_list.append(record)
    SeqIO.write(sequence_list, "{}_output.fasta".format(seq_type), "fasta")

if __name__ == "__main__":
    t_list = get_transcript_list(transcripts2gene)
    dammit_id2transcript_map = get_name_map(name_map)
    parse_fasta_files(cds_fasta_file, "cds", t_list, dammit_id2transcript_map)
    parse_fasta_files(pep_fasta_file, "pep", t_list, dammit_id2transcript_map)






#swissprot_mapping = "swissprot_mapping.tsv"


#def parse_last(annotation): 
#    #function to parse LAST lines in the annotation; returns protein name and Database of origin 
#    if annotation[1] != "LAST": #checks if the line parsing is a orthodb or a sprot line. If not skips 
#        return -1, -1
#    name_index = annotation[8]
#    db = re.search(";database=(OrthoDB|sprot)", name_index)
#    name = re.search("Name=(.*);Target=", name_index)
#   return db.group(1), name.group(1)

#def getsprot2name_map(swissprot_mapping): 
#    sprot2name_map = {}
#    with open(swissprot_mapping, "r") as fh: 
#       table = csv.reader(fh, delimiter = "\t")
#        next(table)
#        for row in table: 
#            sprot2name_map[row[0]] = row[1]
#    return sprot2name_map


#def parse_acc(annotation):
#    #get the acession for swissprot (could be donne with regular expression but this generates 
#    #anoying warning in the IDE) 
#    id = ""
#    flag = False
#    for i in annotation: 
#        if i == "|" : 
#            if flag == True: 
#                return id 
#            else: 
#                flag = True
#        if flag == True and i != "|": 
#           id += i
