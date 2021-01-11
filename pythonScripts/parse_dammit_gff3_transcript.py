### Use this scrip to parse the dammit swissprot output 


import csv 
import re 

dammit_gff = "merged_transcriptome_modified.fasta.dammit.gff3"
dammit_genemap_file = "merged_transcriptome_modified.fasta.dammit.namemap.csv"


def get_transcript2gene_map(file):
    #get transcript2gene dictionary
    transcript2gene_map = {}
    with open(file, "r") as fh: 
        table = csv.reader(fh, delimiter = ",")
        next(table)
        for line in table: 
            transcript2gene_map[line[1]] = line[0]
    return transcript2gene_map



def parse_gff(dammit_gff, dammit_genemap):
    transcript2swissprot = {}
    with open(dammit_gff, "r") as fh: 
        annotation = csv.reader(fh, delimiter = "\t") 
        next(annotation)    
        for line in annotation:
            if line[1] == 'LAST': 
                name_index = line[8]
                db = re.search(";database=(OrthoDB|sprot)", name_index)
                if db.group(1) == "sprot":
                    name = re.search("Name=(.*);Target=", name_index)
                    gene_name = name.group(1)
                    transcript2swissprot[dammit_genemap[line[0]]] = gene_name
    with open("swissprot_transcript_output.tsv", "w") as outfh: 
        for key, value in transcript2swissprot.items(): 
            outfh.write(key)
            outfh.write("\t")
            outfh.write(value)
            outfh.write("\n")



dammit_genemap = get_transcript2gene_map(dammit_genemap_file)
parse_gff(dammit_gff, dammit_genemap)





