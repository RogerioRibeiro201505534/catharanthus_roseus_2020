from Bio import SeqIO 

out_table_file = "transcripts_by_size.tsv" 
input_file = "merged_transcriptome.fasta" 

output_list = []
for record in SeqIO.parse(input_file, "fasta"): 
    output_list.append([record.id, str(len(record.seq))])

with open(out_table_file, "w") as outputfh: 
    for output in output_list:
        outputfh.write("\t".join(output))
        outputfh.write("\n") 