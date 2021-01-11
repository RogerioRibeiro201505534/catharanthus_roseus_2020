##This script takes the trinity transcripts fasta file as an input and renames the transcripts,
#keeping the cluster, gene, and isoform indentifier

from Bio import SeqIO
import re 

transcriptome = "experience_2_transcriptome.fasta"
output_file = "{}_renamed".format(transcriptome)
out_transcript_map = "{}_renamed.gene_trans_map".format(transcriptome)
prefix = "InOut"

def rename_transcripts(transcriptome, output_file, prefix, out_transcript_map): 
    gene2transcript_id = [] 
    output_SeqObj = []
    geneId2newIdmap = {}
    fasta_parser = SeqIO.parse(transcriptome, "fasta")
    gene_counter = 0
    for SeqObj in fasta_parser: 
        regex = re.search("(TRINITY_DN[0-9]+_c[0-9]+_g[0-9]+)_(i[0-9]+)", SeqObj.id)
        gene = regex.group(1)
        isoform_id = regex.group(2)
        if gene in geneId2newIdmap:
            gene_numb = geneId2newIdmap[gene]
        else: 
            gene_counter = gene_counter + 1 
            geneId2newIdmap[gene] = gene_counter
            gene_numb = gene_counter
        SeqObj.id = "{}{}_{}".format(prefix, gene_numb,isoform_id)
        output_SeqObj.append(SeqObj)
        gene2transcript_id.append("{}{}\t{}\n".format(prefix, gene_numb, SeqObj.id))
    SeqIO.write(output_SeqObj, output_file, "fasta")
    with open(out_transcript_map, "w") as fh:
        fh.writelines(gene2transcript_id) 
    


rename_transcripts(transcriptome, output_file, prefix, out_transcript_map)