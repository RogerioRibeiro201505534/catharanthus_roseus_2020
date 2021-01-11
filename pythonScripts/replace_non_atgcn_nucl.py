from Bio import SeqIO
from Bio.Seq import Seq 


input_file = "merged_transcriptome.fasta"
output_file = "merged_transcriptome_modified.fasta"

fasta_parser = SeqIO.parse(input_file, "fasta")
new_seq_list = []
non_atgcn_counter = 0
for SeqObj in fasta_parser: 
    new_seq = "" 
    for nucl in SeqObj.seq: 
        if nucl.upper() in ["A", "T", "C", "G", "N"]: 
            new_seq += nucl 
        else: 
            new_seq += "N"
            non_atgcn_counter += 1
    SeqObj.seq = Seq(new_seq)
    new_seq_list.append(SeqObj)
print(non_atgcn_counter)
SeqIO.write(new_seq_list, output_file, "fasta")
