#convert the dammit_annotation and annotation.mapped.gaf into a dammit_annotation with goslim 
#the input should be the final annotation already with the included proteome

import csv 
from get_id2term_name_obo import get_id2term_name_obo

plant_go_goslim = "goslim_plant.obo"
goslim_converted_file = "annotations.goslim.gaf"
dammit_annotation = "dammit_annotation_final.tsv"
dammit_annotion_goslim = "dammit_annotation_goslim.tsv"




class goslim_annot: 
    def __init__(self, gene_name): 
        #this class is used to store all the information related to the goslim annotation from a gene
        self.gene_name = gene_name
        self.goslimBP = []
        self.goslimMF = []
        self.goslimCC = [] 

    def add_annotation(self, annot, domain_type): 
        if domain_type == "BP":
            if annot not in self.goslimBP: 
                self.goslimBP.append(annot)
        
        if domain_type == "MF": 
            if annot not in self.goslimMF: 
                self.goslimMF.append(annot)
        
        if domain_type == "CC": 
            if annot not in self.goslimCC: 
                self.goslimCC.append(annot)

def make_goslim_dictionary(goslim_converted_file): 
    with open(goslim_converted_file, "r") as fh: 
        goslim_file = csv.reader(fh, delimiter = "\t")
        i = 0
        while i != 5: 
            next(goslim_file)
            i+= 1 
        gene2goslim_map = {}
        for row in goslim_file:
            if row[1] not in gene2goslim_map.keys(): 
                gene2goslim_map[row[1]] = goslim_annot(row[1])
            
            if row[8] == "BP": 
                gene2goslim_map[row[1]].add_annotation(row[4], "BP")
            if row[8] == "MF": 
                gene2goslim_map[row[1]].add_annotation(row[4], "MF")
            if row[8] == "CC": 
                gene2goslim_map[row[1]].add_annotation(row[4], "CC")
    return gene2goslim_map

def go_annotation_file(dammit_annotation, dammit_annotion_goslim, gene2goslim_map, go2term_map): 
    with open(dammit_annotation, "r") as fh: 
        dammit_annot = csv.reader(fh, delimiter = "\t")
        next(dammit_annot)
        
        header = ['gene id', 'orthodb', 'Swissprot', 'Swissprot.name', 'pfam domain', 'pfam domain description', 'pfam acession', 'infernal', 'transcript', "goslim BP", "goslim BP Term", "goslim MF", "goslim MF Term", "goslim CC", "goslim CC Term", "proteome_location", "cathacyc_annotation"]
        outputfh = open(dammit_annotion_goslim, "w")
        outputfh.write("\t".join(header))
        outputfh.write("\n")
        outputfh.close()
        
        for line in dammit_annot: 
            row = line[:9]
            proteome_annotation = line[15:]
            if row[0] in gene2goslim_map.keys():
                goslim_annot = gene2goslim_map[row[0]] 
                BP_term_names = []
                MF_term_names = []
                CC_term_names = []
                if goslim_annot.goslimBP != []: 
                    for BP in goslim_annot.goslimBP: 
                        BP_term_names.append(go2term_map[BP])
                if goslim_annot.goslimMF != []: 
                    for MF in goslim_annot.goslimMF: 
                        MF_term_names.append(go2term_map[MF])
                if goslim_annot.goslimCC != []: 
                    for CC in goslim_annot.goslimCC: 
                        CC_term_names.append(go2term_map[CC])

                row.append(";".join(goslim_annot.goslimBP))
                row.append(";".join(BP_term_names))
                row.append(";".join(goslim_annot.goslimMF))
                row.append(";".join(MF_term_names))
                row.append(";".join(goslim_annot.goslimCC))
                row.append(";".join(CC_term_names))
            else: 
                for i in range(6): 
                    row.append("")
            row.append(proteome_annotation[0])
            row.append(proteome_annotation[1])
            #finally write the row into a new file
            outputfh = open(dammit_annotion_goslim, "a")
            outputfh.write("\t".join(row))
            outputfh.write("\n")
            outputfh.close()

            
            


gene2goslim_map = make_goslim_dictionary(goslim_converted_file)
go2term_map = get_id2term_name_obo(plant_go_goslim)
go_annotation_file(dammit_annotation, dammit_annotion_goslim, gene2goslim_map, go2term_map)