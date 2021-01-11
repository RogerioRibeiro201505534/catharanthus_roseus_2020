import csv 
from Bio import SeqIO


input_file = "transcript2gene.tsv"
input_fasta = "merged_transcriptome_filtered.fasta"


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


def get_longest_transcript(map, input_fasta): 
    #Takes the map file and gets the longest transcript from each gene_id 
    fasta_file_seq = []
    sequence_map = {}
    fasta_sequences = SeqIO.parse(open(input_fasta), "fasta")
    for sequence in fasta_sequences: 
        sequence_len = len(sequence.seq)
        sequence_map[sequence.id] = (sequence, sequence_len)
    for items in map.items():
        transcript_list = []
        for transcript in items[1]: 
            transcript_list.append(sequence_map[transcript])
        sorted_list = sorted(transcript_list, key = lambda x:x[1], reverse = True)
        fasta_file_seq.append(sorted_list[0][0]) 
    SeqIO.write(fasta_file_seq, "merged_transcriptome_filtered_longest_transcript.fasta", "fasta")
        
    



map = make_map(input_file)
get_longest_transcript(map, input_fasta)
