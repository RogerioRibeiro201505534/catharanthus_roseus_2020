### Use this scrip to parse the dammit output and merge transcripts annotations, as well as get the swissprot 
#gene id to parse

import csv 
import re 

dammit_sum_gff = "merged_transcriptome_modified.fasta.dammit.gff3"
transcript2genefile = "transcript2gene.tsv"
dammit_genemap = "merged_transcriptome_modified.fasta.dammit.namemap.csv"
output_file = "dammit_annotation_output.tsv"
annotation_to_keep = ["HMMER", "LAST", "Infernal"]
swissprot_id_list_file = "swissprot_id_query.txt"

def get_transcript2gene_map(file): 
    #get transcript2gene dictionary
    transcript2gene_map = {}
    with open(file, "r") as fh: 
        table = csv.reader(fh, delimiter = "\t")
        for line in table: 
            transcript2gene_map[line[0]] = line[1]
    return transcript2gene_map

def get_dammit_id2transcript_map(file): 
    #dammit_id 2 transcript id dictionary 
    dammit_id2transcript_map = {}
    with open(file, "r") as fh: 
        table = csv.reader(fh, delimiter = ",")
        next(table)
        for line in table: 
            dammit_id2transcript_map[line[1]] = line[0]
    return dammit_id2transcript_map



def parse_gff(dammit_sum_gff, transcript2gene_map, dammit_id2transcript_map):
    #Parse relevant annotations to a dictionary where the keys are gene names 
    #and replace the transcripts name in the annotation with the real transcript names  
    gene2_annotation_map = {}
    with open(dammit_sum_gff, "r") as fh: 
        annotation = csv.reader(fh, delimiter = "\t") 
        next(annotation)    
        for line in annotation:
            transcript_id = dammit_id2transcript_map[line[0]]
            if transcript_id in transcript2gene_map.keys(): 
                gene_id = transcript2gene_map[transcript_id]
            else:
                continue
            line[0] = transcript_id
            if gene_id in gene2_annotation_map.keys():
                gene2_annotation_map[gene_id].append(line)
            else: 
                gene2_annotation_map[gene_id] = [line]
    #note that not all genes are annotated. Therefore genes, with no annotation have to be added 
    #in the map with empty values 
    for gene in transcript2gene_map.values(): 
        if gene not in gene2_annotation_map.keys(): 
            gene2_annotation_map[gene] = []
    return gene2_annotation_map


def parse_infernal(annotation):
    #function to parse infernal lines in the annotation; returns NAME of the hits 
    name_index = annotation[8]
    name = re.search("Name=(.*);Target=", name_index)
    return name.group(1)

def parse_last(annotation): 
    #function to parse LAST lines in the annotation; returns protein name and Database of origin 
    name_index = annotation[8]
    db = re.search(";database=(OrthoDB|sprot)", name_index)
    name = re.search("Name=(.*);Target=", name_index)
    return db.group(1), name.group(1)

def parse_hmmer(annotation): 
    #function to parse pfam lines in the annotation; returns name, description and pfam id 
    name_index = annotation[8]
    name = re.search("Name=(.*);Target=", name_index)
    note = re.search("Note=(.*);accuracy=", name_index)
    pfamref = re.search("Dbxref=\"Pfam:(.*)\"", name_index)
    return name.group(1), note.group(1), pfamref.group(1)

def append_in_output(annotation, output_files): 
    last_index = len(annotation)
    text_to_write = ""
    for collum_numb in range(len(annotation)): 
        text_to_write += ";".join(annotation[collum_numb])
        if collum_numb == last_index - 1: 
            text_to_write += "\n"
        else: 
            text_to_write += "\t"
    with open(output_file, "a") as fh: 
        fh.write(text_to_write)
        
def parse_acc(annotation):
    #get the acession for swissprot (could be donne with regular expression but this generates anoying warning) 
    #in the IDE
    id = ""
    flag = False
    for i in annotation: 
        if i == "|" : 
            if flag == True: 
                return id 
            else: 
                flag = True
        if flag == True and i != "|": 
            id += i



def output_annotation(gene2_annotation_map, output_file, annotation_to_keep): 
    #outputs the a file with the annotation merged by gene as well as a swissprot id, which has to be used in order to have the swissprot annotation name
    header = "gene id\torthodb\tSwissprot\tpfam domain\tpfam domain description\tpfam acession\tinfernal\ttranscript\n"
    i = 0 
    swiss_prot_id_list = [] # This list will keep unique swissprot_ids (will be used latter to get go terms)
    with open(output_file, "w") as fh: 
        fh.write(header)
    for gene in gene2_annotation_map.items():
        i += 1 
        annotation_line = [[],[],[],[],[],[],[],[]]
        annotation_line[0].append(gene[0])
        for transcript_annotation in gene[1]: 
            if transcript_annotation[0] not in annotation_line[7]:
                annotation_line[7].append(transcript_annotation[0]) #append transcript id in index 7
            if transcript_annotation[1] == "Infernal": #index 6 in annotation_line list (starting from 0)
                annotation = parse_infernal(transcript_annotation)
                if annotation not in annotation_line[6]: 
                    annotation_line[6].append(annotation)
            if transcript_annotation[1] == "LAST": 
                db, annotation = parse_last(transcript_annotation)
                if db == "OrthoDB": #index 1 in annotation_line list
                    if annotation not in annotation_line[1]: 
                        annotation_line[1].append(annotation)
                if db == "sprot": 
                    if annotation not in annotation_line[2]: #index 2 in annotation_line_list 
                        annotation_line[2].append(annotation)  
                    swiss_prot_id = parse_acc(annotation)
                    if swiss_prot_id not in swiss_prot_id_list:  # get Uniprot acession id for latter use 
                        swiss_prot_id_list.append(swiss_prot_id)       
            if transcript_annotation[1] == "HMMER": 
                domain, note, pfamref = parse_hmmer(transcript_annotation)
                if pfamref not in  annotation_line[5]:  #Pfam ref stored in the index 5 
                    annotation_line[5].append(pfamref)
                    annotation_line[3].append(domain) #pfam domain name stored in index 3 
                    annotation_line[4].append(note) #pfam description stored in index 4 
        append_in_output(annotation_line, output_file)
    
    return swiss_prot_id_list





if __name__ == "__main__":
    transcript2gene_map = get_transcript2gene_map(transcript2genefile)    
    dammit_id2transcript_map = get_dammit_id2transcript_map(dammit_genemap)
    gene2_annotation_map = parse_gff(dammit_sum_gff, transcript2gene_map, dammit_id2transcript_map)
    swiss_prot_id_list = output_annotation(gene2_annotation_map, output_file, annotation_to_keep)
    with open(swissprot_id_list_file, "w") as fh: 
        for swiss_id in swiss_prot_id_list: 
            fh.write("{}\n".format(swiss_id))

            
