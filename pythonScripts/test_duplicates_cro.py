#use this script to test how many "Cro locus are in the same gene id in the genome"

import csv 
import re 

gff_compare_file = "merged_transcriptome.annotated.gtf"

def count_cro_in_gene_ids(gff_compare_file): 
    with open(gff_compare_file) as fh: 
        reader = csv.reader(fh, delimiter = "\t")
        dictionary_by_gene_id = {}
        dictionary_by_gene_name = {}
        for line in reader: 
            if line[2] == "transcript": 
                gff_line_9 = line[8].split(";")
                if "gene_name" in gff_line_9[2]:
                    #process cases where there is a reference gene_name
                    re_transcript = re.search("transcript_id \"(.*)\"", gff_line_9[0])
                    transcript_id = re_transcript.group(1)
                    re_gene = re.search("gene_id \"(.*)\"", gff_line_9[1])
                    gene_id = re_gene.group(1)
                    re_gene_name = re.search("gene_name \"(.*)\"", gff_line_9[2])
                    gene_name = re_gene_name.group(1)

                else: 
                    re_transcript = re.search("transcript_id \"(.*)\"", gff_line_9[0])
                    transcript_id = re_transcript.group(1)
                    re_gene = re.search("gene_id \"(.*)\"", gff_line_9[1])
                    gene_id = re_gene.group(1)
                    gene_name = gene_id
                
                if gene_id not in dictionary_by_gene_id.keys(): 
                    dictionary_by_gene_id[gene_id] = [transcript_id]
                else: 
                    dictionary_by_gene_id[gene_id].append(transcript_id)
                
                if gene_name not in dictionary_by_gene_name.keys(): 
                    dictionary_by_gene_name[gene_name] = [transcript_id]
                else: 
                    dictionary_by_gene_name[gene_name].append(transcript_id)

        #After processing all the rows in the transcript file
        fh = open("duplicate_Cro_ids.txt", "w")
        dupe_cro_counter = 0 
        cro_counter_dic = {}
        for key in dictionary_by_gene_id.keys():
            cro_counter = 0
            for transcript in dictionary_by_gene_id[key]: 
                if "CRO_T" in transcript: 
                    cro_counter += 1 
            if cro_counter >= 2:
                dupe_cro_counter += 1 
                print(key)
                fh.write(key + "\n")
            if str(cro_counter) not in cro_counter_dic.keys(): 
                cro_counter_dic[str(cro_counter)] = 1
            else: 
                cro_counter_dic[str(cro_counter)] += 1
        print(cro_counter_dic)
        print(dupe_cro_counter)
        fh.close()

        with open("geneName2transcript_id.tsv", "w") as fh: 
            for key in dictionary_by_gene_name.keys(): 
                for transcript_id in dictionary_by_gene_name[key]: 
                    fh.write(key + "\t" + transcript_id)
                    fh.write("\n")
        

 
                     




count_cro_in_gene_ids(gff_compare_file)