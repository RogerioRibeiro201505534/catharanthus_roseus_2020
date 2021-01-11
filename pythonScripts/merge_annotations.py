#Use this script to merge annotations, per gene id 

import csv

gene_2_transcript_table = "transcript2gene_id.csv"
input_transcript = "00_Venn_geral_kallisto_uniref90.txt"
output_gene =  "00_Venn_geral_kallisto_unifref90_merged.txt"

def get_gene2transcript_map(file): 
    #get gene2transcript dictionary 
    gene_2_transcript_map = {}
    with open(file, "r") as fh: 
        table = csv.reader(fh, delimiter = "\t")
        for line in table: 
            if line[1] in gene_2_transcript_map.keys(): 
                gene_2_transcript_map[line[1]].append(line[0])
            else: 
                gene_2_transcript_map[line[1]] = [line[0]]
    return gene_2_transcript_map

def get_annot_map(input_transcript, header=True): 
    #get dictionary where the key are transcript_id´s 
    annot_map = {}
    with open(input_transcript, "r") as fh: 
        annotation = csv.reader(fh, delimiter = "\t")
        if header == True: 
            next(annotation)
        for line in annotation: 
            annot_map[line[0]] = line[0:13]
    return annot_map

def merge_annotation(gene_annot_list): 
    #merge annotations based on gene_id, return result as list 
    merged_annotation = [[],[],[],[],[],[],[],[],[],[],[],[],[]]
    for transcript in gene_annot_list: 
        i = 0 
        for annot in transcript: 
            if ";" in annot:# if there are several, semicollon separated values 
                annot = annot.split(";")
                for info in annot: 
                    if info not in merged_annotation[i]: 
                        merged_annotation[i].append(info)
            else: # if there aren´t several, semicollon separated vallues 
                if annot == '': 
                    pass
                elif annot not in merged_annotation[i]: 
                    merged_annotation[i].append(annot)
            i += 1 
    return merged_annotation



def main_function(gene_2_transcript_map,annot_map, output_gene):
    with open(output_gene, "w") as fh: 
        fh.write("ID\tTRANSCRIPTS\tGENENAME	DESCRIPTION	ENZYME	GO(P)ID	GO(P)NAME	GO(F)ID	GO(F)NAME	GO(C)ID	GO(C)NAME	KEYWORD	PATHWAY	GOSLIM\n")
    for gene in gene_2_transcript_map.keys(): 
        gene_annot_list = []
        for transcript in gene_2_transcript_map[gene]:
            if transcript in annot_map:
                gene_annot_list.append(annot_map[transcript])
        if gene_annot_list != []: 
            merged_annotation = merge_annotation(gene_annot_list)
            merged_annotation.insert(0, [gene])
            output_annotation = []
            for collum in merged_annotation: 
                output_annotation.append(";".join(collum))
            with open(output_gene, "a") as fh: 
                fh.write("\t".join(output_annotation))
                fh.write("\n")
                    

if __name__ == "__main__":
    gene_2_transcript_map = get_gene2transcript_map(gene_2_transcript_table)
    annot_map = get_annot_map(input_transcript, header = True)
    main_function(gene_2_transcript_map,annot_map, output_gene)