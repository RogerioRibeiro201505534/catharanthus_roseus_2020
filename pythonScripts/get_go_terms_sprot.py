##use this script to put together the go_terms annotation from swissprot into the table 


import csv
import re 

uniprot_file = "swissprot_mapping.tsv"
output_file = "dammit_annotation_with_go.tsv"
input_file = "dammit_annotation_output.tsv"

def read_uniprot_output(file): 
    #get a uniprot2goterms dictionary 
    #key: gene_acc 
    #values [per index] protein_name [0], GOBP [1], GOBP term [2] GOMF [1], GOMF term [2] GOCC [1], GOCC term [2] 
    #go_term_dictionary = {}   
    accesion_annotation = {} 
    with open(file, "r") as fh:
        table = csv.reader(fh, delimiter = "\t")
        next(table)
        for row in table:
            accesion_annotation[row[0]] = [[row[1]], [], [], [], [], [], []] 
            if row[2] != "": 
                go_terms_BP = row[2].split("; ")
                for go_term in go_terms_BP: 
                    go_term_name, go_term_acc = parse_go_term(go_term)
                    accesion_annotation[row[0]][1].append(go_term_name)
                    accesion_annotation[row[0]][2].append(go_term_acc)

            if row[4] != "": 
                go_terms_CC = row[4].split("; ")
                for go_term in go_terms_CC: 
                    go_term_name, go_term_acc = parse_go_term(go_term)
                    accesion_annotation[row[0]][3].append(go_term_name)
                    accesion_annotation[row[0]][4].append(go_term_acc)

            if row[3] != "": 
                go_terms_MF = row[3].split("; ")
                for go_term in go_terms_MF: 
                    go_term_name, go_term_acc = parse_go_term(go_term)
                    accesion_annotation[row[0]][5].append(go_term_name)
                    accesion_annotation[row[0]][6].append(go_term_acc)
                    
    return accesion_annotation

def parse_go_term(go_term):  
    info = re.search("(.+) \[(GO:.+)\]", go_term)
    return info.group(1), info.group(2)
            

def insert_go_notation(input_file, output_file, go_annot_dictionary): 
    outputfh = open(output_file, "w")
    outputfh.write("gene id\torthodb\tSwissprot\tSwissprot name\tpfam domain\tpfam domain description\tpfam acession\tinfernal\ttranscript\tGO(BP)\tGO(BP) terms\tGO(MF)\tGO(MF) terms\tGO(CC)\tGO(CC) terms\n")
    outputfh.close()
    with open(input_file) as inputfh: 
        input_table = csv.reader(inputfh, delimiter = "\t")
        next(input_table)
        for row in input_table: 
            output_row = row [:]
            GO_BP = []
            GO_BP_TERMS = []
            GO_MF = []
            GO_MF_TERMS = []
            GO_CC = []
            GO_CC_TERMS = []
            gene_name_list = []
            if row[2] != "": #skip these steps if there is no swissprot annotation 
                swissprot_annots = row[2].split(";")
                for swissprot_annot in swissprot_annots:
                    get_uniprot_id = re.search("\|(.+)\|", swissprot_annot)
                    uniprot_id = get_uniprot_id.group(1)
                    uniprot_id_annotation = go_annot_dictionary[uniprot_id]
                    gene_name_list.append(uniprot_id_annotation[0][0])
                    if uniprot_id_annotation[1] != []: 
                        for go_term, go_name in zip(uniprot_id_annotation[2], uniprot_id_annotation[1]):
                            if go_term not in GO_BP: 
                                GO_BP.append(go_term)
                                GO_BP_TERMS.append(go_name)
                    if uniprot_id[3] != []:
                        for go_term, go_name in zip(uniprot_id_annotation[4], uniprot_id_annotation[3]):
                            if go_term not in GO_CC: 
                                GO_CC.append(go_term)
                                GO_CC_TERMS.append(go_name)
                    if uniprot_id[5] != []:
                        for go_term, go_name in zip(uniprot_id_annotation[6], uniprot_id_annotation[5]):
                            if go_term not in GO_MF: 
                                GO_MF.append(go_term)
                                GO_MF_TERMS.append(go_name)
            #Append rows to output file row
            output_row.insert(3,";".join(gene_name_list))    
            output_row.append(";".join(GO_BP))
            output_row.append(";".join(GO_BP_TERMS))
            output_row.append(";".join(GO_MF))
            output_row.append(";".join(GO_MF_TERMS))
            output_row.append(";".join(GO_CC))
            output_row.append(";".join(GO_CC_TERMS))
            #Write each row to the output file
            outputfh = open(output_file, "a")
            outputfh.write("\t".join(output_row))
            outputfh.write("\n")
            outputfh.close()
            




go_annot_dictionary = read_uniprot_output(uniprot_file)
insert_go_notation(input_file, output_file, go_annot_dictionary)