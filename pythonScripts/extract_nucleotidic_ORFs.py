### Get the longest ORF from the output 

from Bio import SeqIO

input_file = "Catro_mRNA_update.tfa_annotated_ORF_v2.txt"
output_file = "Catro_mRNA_update.tfa_annotated_ORFs_longest_ORF_v2.fasta"

def parse_file_longest_transcript(input_file, output_file): 
    longest_sequence_map = {}
    fasta_file = SeqIO.parse(input_file, "fasta")
    for sequence in fasta_file: 
        sequence_id = sequence.id.split("_")[0]
        if sequence_id in longest_sequence_map.keys(): 
            if len(sequence.seq) >= len(longest_sequence_map[sequence_id].seq): 
                longest_sequence_map[sequence_id] = sequence
        else: 
            longest_sequence_map[sequence_id] = sequence

    sequence_list = []
    for value in longest_sequence_map.values(): 
        sequence_list.append(value)
    
    SeqIO.write(sequence_list, output_file, "fasta")



parse_file_longest_transcript(input_file, output_file)