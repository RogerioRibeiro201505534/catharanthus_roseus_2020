# Use this scrip to parse the trinotate output and merge transcripts annotations. 
# The different transcoder prediction will be regarded as different transcript for partical purposes
#  and as such will be merged in the same row 

import csv 
import re

trinotate_tsv_test = "Trinotate.tsv"
output = "trinotate_parsed.tsv"


class gene_annotation: 
    #annotation class, storing several annotation type
    #Todo method to write in file here? 
    def __init__(self):     
        self.transcripts_id_list = []
        self.swissprot_id_list = []
        self.gene_name_list = []
        self.rnammer_list = []
        self.prot_id_list = []
        self.pfam_acc_list = []
        self.pfam_name_list = []
        self.pfam_description_list = []
        self.thmmer_list = []
        self.SignalP_list = []
        self.eggNog_acc_list = []
        self.eggNog_name_list = []
        self.kegg_list = []
        self.GO_BP_list = []
        self.GO_BP_terms_list = []
        self.GO_MF_list = []
        self.GO_MF_terms_list = []
        self.GO_CC_list = []
        self.GO_CC_terms_list = [] 
        

    def set_gene(self, gene): 
        self.gene = gene 

    def write_gene2output(self, output):
        text_to_write = "\n"
        text_to_write += self.gene + "\t"
        text_to_write += ";".join(self.transcripts_id_list) + "\t"
        text_to_write += str(len(self.prot_id_list)) + "\t"
        text_to_write += ";".join(self.swissprot_id_list) + "\t"
        text_to_write += ";".join(self.gene_name_list) + "\t"
        text_to_write += ";".join(self.rnammer_list) + "\t"
        text_to_write += ";".join(self.pfam_acc_list) + "\t"
        text_to_write += ";".join(self.pfam_name_list) + "\t"
        text_to_write += ";".join(self.pfam_description_list) + "\t"
        text_to_write += ";".join(self.SignalP_list) + "\t"
        text_to_write += ";".join(self.thmmer_list) + "\t"
        text_to_write += ";".join(self.eggNog_acc_list) + "\t"
        text_to_write += ";".join(self.eggNog_name_list) + "\t"
        text_to_write += ";".join(self.kegg_list) + "\t"
        text_to_write += ";".join(self.GO_BP_list) + "\t"
        text_to_write += ";".join(self.GO_BP_terms_list) + "\t"      
        text_to_write += ";".join(self.GO_MF_list) + "\t"
        text_to_write += ";".join(self.GO_MF_terms_list) + "\t"      
        text_to_write += ";".join(self.GO_CC_list) + "\t"
        text_to_write += ";".join(self.GO_CC_terms_list) + "\t"
        with open(output, "a") as fh: 
            fh.write(text_to_write)

    

def get_gene2transcript_annot_map(trinotate_tsv): 
    #Get a dictionary with gene keys and transcript annotations as the value
    gene2transcript_annot = {}
    with open(trinotate_tsv, "r") as fh: 
        annotation = csv.reader(fh, delimiter = "\t")
        next(annotation)
        for row in annotation: 
            if row[0] in gene2transcript_annot.keys(): 
                gene2transcript_annot[row[0]].append(row[1:])
            else: 
                gene2transcript_annot[row[0]] = [row[1:]]
    return gene2transcript_annot

def parse_swissprot_collum(swissprot_id): 
    #Parses sprot_best_blastp_hit and sprot_best_blastx_hit 
    search_full_name = re.search("Full=(.*);", swissprot_id)
    inter_name_string = search_full_name.group(1)
    gene_name = "" 
    for i in inter_name_string: 
        if i == ";"or i == "{": 
            return swissprot_id.split("^")[0], gene_name 
        gene_name += i

def rnammer_parser(rnammer_annot): 
    return rnammer_annot.split("^")[0]

def parse_signalP(signal_annot): 
    return signal_annot.replace("^", ";")

def pfam_parser(pfam_annot): 
    pfam_accs = []
    pfam_names = []
    pfam_descriptions = []
    pfam_hit_list = pfam_annot.split("`")
    for hit in pfam_hit_list: 
        hit_info = hit.split("^")
        pfam_accs.append(hit_info[0])
        pfam_names.append(hit_info[1])
        pfam_descriptions.append(hit_info[2])
    return pfam_accs, pfam_names, pfam_descriptions

def parse_thmmer(thmmer_row): 
    return "({})".format(thmmer_row.replace("^", ";"))    

def  parse_eggnog(eggnog_annot): 
    eggnogs_accs = []
    eggnogs_names = []
    eggnog_hit_list = eggnog_annot.split("`")
    for hit in eggnog_hit_list: 
        hit_info = hit.split("^")
        eggnogs_accs.append(hit_info[0])
        eggnogs_names.append(hit_info[1])
    return eggnogs_accs, eggnogs_names

def parse_kegg(kegg_annot): 
    return kegg_annot.split("`")

def parse_go_terms(go_annot):
    go_BPs = []
    go_BP_terms = []
    go_MFs = []
    go_MF_terms = []
    go_CCs = []
    go_CC_terms = []
    go_hit_list = go_annot.split("`")
    for hit in go_hit_list: 
        hit_info = hit.split("^")
        if hit_info[1] == "biological_process": 
            if hit_info[0] not in go_BPs:
                go_BPs.append(hit_info[0])
                go_BP_terms.append(hit_info[2])
        if hit_info[1] == "molecular_function": 
            if hit_info[0] not in go_MFs:
                go_MFs.append(hit_info[0])
                go_MF_terms.append(hit_info[2])
        if hit_info[1] == "cellular_component": 
            if hit_info[0] not in go_CCs:
                go_CCs.append(hit_info[0])
                go_CC_terms.append(hit_info[2])
                
    return go_BPs, go_BP_terms, go_MFs, go_MF_terms, go_CCs, go_CC_terms

def merge_annotation_by_geneid(gene2transcript_annot, output): 
    fh = open(output, "w")
    fh.write("gene_id\ttranscript_id\tnº Transcoder peptides\tswissprot\tgene_name\tRNAMMER\tPfam\tpfam_name\tpfam_description\tSignalP\tTmHMM\teggnog\teggnog term\tKegg\tGO(BP)\tGO(BP)term\tGO(MF)\tGO(MF)term\tGO(CC)\tGO(CC)term")
    fh.close()
    for gene, transcripts_annot in gene2transcript_annot.items():
        #create object from the gene annotation class (class with a bunch of list and a write method)
        gene_annot = gene_annotation()
        gene_annot.set_gene(gene)
        for trinity_row in transcripts_annot: #For each row in annotation
            
            if trinity_row[0] not in gene_annot.transcripts_id_list:
                gene_annot.transcripts_id_list.append(trinity_row[0])
            
            if trinity_row[1] != ".":  
                swissprot_id, gene_name = parse_swissprot_collum(trinity_row[1])
                if swissprot_id not in gene_annot.swissprot_id_list:
                    gene_annot.swissprot_id_list.append(swissprot_id)
                    gene_annot.gene_name_list.append(gene_name)
            
            if trinity_row[2] != ".": 
                rna_annot = rnammer_parser(trinity_row[2])
                if rna_annot not in gene_annot.rnammer_list:
                    gene_annot.rnammer_list.append(rna_annot)
            
            if trinity_row[3] != ".": 
                if trinity_row[3] not in gene_annot.prot_id_list: 
                    gene_annot.prot_id_list.append(trinity_row[3])
            
            if trinity_row[5] != ".":  
                swissprot_id, gene_name = parse_swissprot_collum(trinity_row[5])
                if swissprot_id not in gene_annot.swissprot_id_list:
                    gene_annot.swissprot_id_list.append(swissprot_id)
                    gene_annot.gene_name_list.append(gene_name)

            if trinity_row[6] != ".": 
                pfam_accs, pfam_names, pfam_descriptions = pfam_parser(trinity_row[6])
                for i in range(len(pfam_accs)): 
                    if pfam_accs[i] not in gene_annot.pfam_acc_list: 
                        gene_annot.pfam_acc_list.append(pfam_accs[i])
                        gene_annot.pfam_name_list.append(pfam_names[i])
                        gene_annot.pfam_description_list.append(pfam_descriptions[i])

            if trinity_row[7] != ".": 
                signalp = parse_signalP(trinity_row[7])
                if signalp not in gene_annot.SignalP_list: 
                    gene_annot.SignalP_list.append(signalp)
            
            if trinity_row[8] != ".": 
                thmmer_annot = parse_thmmer(trinity_row[8])
                if thmmer_annot not in gene_annot.thmmer_list: 
                    gene_annot.thmmer_list.append(thmmer_annot)
        
            if trinity_row[9] != ".": 
                eggNog_accs, eggNog_names = parse_eggnog(trinity_row[9])
                for i in range(len(eggNog_accs)): 
                    if eggNog_accs[i] not in gene_annot.eggNog_acc_list: 
                        gene_annot.eggNog_acc_list.append(eggNog_accs[i])
                        gene_annot.eggNog_name_list.append(eggNog_names[i])

            if trinity_row[10] != ".": 
                keggs = parse_kegg(trinity_row[10])
                for keg in keggs: 
                    if keg not in gene_annot.kegg_list: 
                        gene_annot.kegg_list.append(keg)

            if trinity_row[11] != ".": 
                go_BPs, go_BP_terms, go_MFs, go_MF_terms, go_CCs, go_CC_terms = parse_go_terms(trinity_row[11])
                for i in range(len(go_BPs)):
                    if go_BPs[i] not in gene_annot.GO_BP_list: 
                        gene_annot.GO_BP_list.append(go_BPs[i])
                        gene_annot.GO_BP_terms_list.append(go_BP_terms[i])
                for i in range(len(go_MFs)): 
                    if go_MFs[i] not in gene_annot.GO_MF_list: 
                        gene_annot.GO_MF_list.append(go_MFs[i])
                        gene_annot.GO_MF_terms_list.append(go_MF_terms[i])
                for i in range(len(go_CCs)):
                    if go_CCs[i] not in gene_annot.GO_CC_list: 
                        gene_annot.GO_CC_list.append(go_CCs[i])
                        gene_annot.GO_CC_terms_list.append(go_CC_terms[i])
            
            if trinity_row[12] != ".":
                go_BPs, go_BP_terms, go_MFs, go_MF_terms, go_CCs, go_CC_terms = parse_go_terms(trinity_row[12])
                for i in range(len(go_BPs)):
                    if go_BPs[i] not in gene_annot.GO_BP_list: 
                        gene_annot.GO_BP_list.append(go_BPs[i])
                        gene_annot.GO_BP_terms_list.append(go_BP_terms[i])
                for i in range(len(go_MFs)): 
                    if go_MFs[i] not in gene_annot.GO_MF_list: 
                        gene_annot.GO_MF_list.append(go_MFs[i])
                        gene_annot.GO_MF_terms_list.append(go_MF_terms[i])
                for i in range(len(go_CCs)):
                    if go_CCs[i] not in gene_annot.GO_CC_list: 
                        gene_annot.GO_CC_list.append(go_CCs[i])
                        gene_annot.GO_CC_terms_list.append(go_CC_terms[i])

            if trinity_row[13] != ".":
                go_BPs, go_BP_terms, go_MFs, go_MF_terms, go_CCs, go_CC_terms = parse_go_terms(trinity_row[13])
                for i in range(len(go_BPs)):
                    if go_BPs[i] not in gene_annot.GO_BP_list: 
                        gene_annot.GO_BP_list.append(go_BPs[i])
                        gene_annot.GO_BP_terms_list.append(go_BP_terms[i])
                for i in range(len(go_MFs)): 
                    if go_MFs[i] not in gene_annot.GO_MF_list: 
                        gene_annot.GO_MF_list.append(go_MFs[i])
                        gene_annot.GO_MF_terms_list.append(go_MF_terms[i])
                for i in range(len(go_CCs)):
                    if go_CCs[i] not in gene_annot.GO_CC_list: 
                        gene_annot.GO_CC_list.append(go_CCs[i])
                        gene_annot.GO_CC_terms_list.append(go_CC_terms[i])
        
        gene_annot.write_gene2output(output)
        

if __name__ == "__main__":
    gene2transcript_annot = get_gene2transcript_annot_map(trinotate_tsv_test)
    merge_annotation_by_geneid(gene2transcript_annot, output)
    