import csv 
import re

busco_file = "full_table_from StringTie_transcriptome.tsv"
gtf_file = "merged_transcriptome.gtf"

def read_Busco_file(busco_file): 
    with open(busco_file, 'r') as bfh: 
        busco_file = csv.reader(bfh, delimiter = "\t")
        next(busco_file)
        next(busco_file)
        next(busco_file)
        next(busco_file)
        busco_dictionary = {}
        for line in busco_file:
            if line[1] != "Missing": 
                gene_name = "{}_{}_({})".format(line[2], line[0], line[1])
                #dictionary with gene_name and calssification 
                busco_dictionary[line[2]] = [gene_name, line[1]]
    return busco_dictionary

def parse_gtf_file(gtf_file, busco_dictionary): 
    important_lines = [] 
    with open(gtf_file, 'r') as gfh: 
        gtf_coord = csv.reader(gfh, delimiter = "\t")
        next(gtf_coord)
        next(gtf_coord)
        for line in gtf_coord: 
            if line[2] == 'transcript': 
                re_search = re.search("; transcript_id \"(.+)\";", line[8])
                transcript = re_search.group(1)
                if transcript in busco_dictionary.keys(): 
                    line[8] = busco_dictionary[transcript][0]
                    important_lines.append(line)
    with open("new_gff_file.gff", "w") as fhw: 
        for line in important_lines: 
            fhw.write("\t".join(line))
            fhw.write("\n")



if __name__ == "__main__":
    busco_dictionary = read_Busco_file(busco_file)
    parse_gtf_file(gtf_file, busco_dictionary)