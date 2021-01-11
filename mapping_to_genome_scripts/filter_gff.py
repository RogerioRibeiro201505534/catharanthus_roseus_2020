#this script filters the gffcompare file, only kepping certain tags -> List of tags to keep 
#and addin a gene_name tag where there is none (new genes)
#also generates a list of transcript to gene association as well as a list of new transcripts -> class code "u"
import csv 
import re



input_file = "mapping_to_genome/04_gffcompare/merged_transcriptome.annotated.gtf"
output_file ="mapping_to_genome/04_gffcompare/merged_transcriptome.filtered.gtf"

list_of_taggs_to_keep = ["j", "=", "k", "m", "n", "u"]

def helper_function(collum9, atribute): 
    #small helper function to retrieve the class code
    if atribute == "class_code": 
        code_search = re.search("; class_code \"([a-z]|=)\"", collum9)
        return code_search.group(1)
    if atribute == "gene_id": 
        atribute_list = collum9.split(";")
        gene_id = atribute_list[1].replace("\"", "")
        transcript_id = atribute_list[0].replace("\"", "")
        return gene_id.split(" ")[2], transcript_id.split(" ")[1]
    if atribute == "gene_name": 
        atribute_list = collum9.split(";")
        gene_name = atribute_list[2].replace("\"", "")
        transcript_id = atribute_list[0].replace("\"", "")
        return gene_name.split(" ")[2], transcript_id.split(" ")[1]

def add_gene_name_to_u(collum9, new_gene_id): 
    #use this function to add a gene_name to the rows where the genes are "new"
    gene_name = " gene_name \"{}\"".format(new_gene_id)
    collum9_list = collum9.split(";")
    collum9_list.insert(2, gene_name)
    collum9 = ";".join(collum9_list)
    return collum9



def filter_gff(input_file, output_file, list_of_taggs_to_keep): 
    with open(input_file, "r") as fh: 
        gene_name_list = []
        transcript_list = []
        new_genes_list = []
        new_transcript_list = []
        filter_flag = False
        list_of_rows_to_keep = []
        gffcompare_file = csv.reader(fh, delimiter = "\t")
        for row in gffcompare_file:
            if row[2] == "transcript" and helper_function(row[8], "class_code") in list_of_taggs_to_keep: 
                filter_flag = True
                if helper_function(row[8], "class_code") == "u": 
                    gene_id, transcript_id = helper_function(row[8], "gene_id")
                    new_genes_list.append(gene_id)
                    new_transcript_list.append(transcript_id)
                    row[8] = add_gene_name_to_u(row[8], gene_id)
                gene_name, transcript_id = helper_function(row[8], "gene_name")
                gene_name_list.append(gene_name)
                transcript_list.append(transcript_id)
                list_of_rows_to_keep.append(row)
                current_gene_name = gene_name
                continue
            if row[2] == "transcript" and helper_function(row[8], "class_code") not in list_of_taggs_to_keep:
                filter_flag = False
                continue
            if row[2] != "transcript" and filter_flag == True: 
                row[8] = add_gene_name_to_u(row[8], current_gene_name)
                list_of_rows_to_keep.append(row)
                continue
    with open(output_file, "w") as output_gtf: 
        for out_row in list_of_rows_to_keep: 
            output_gtf.write("\t".join(out_row)) 
            output_gtf.write("\n")

    with open("mapping_to_genome/04_gffcompare/new_transcript2geneList.txt", "w") as output_t2g_new:
        for new_gene, new_transcript in zip(new_genes_list, new_transcript_list): 
            output_t2g_new.write("{}\t{}".format(new_gene, new_transcript))
            output_t2g_new.write("\n")
    
    with open("mapping_to_genome/04_gffcompare/new_genes.txt", "w") as output_genes_new: 
        #this will change the order, which it does not matter a lot
        output_genes_new.write("\n".join(list(set(new_genes_list))))

    with open("mapping_to_genome/04_gffcompare/transcript2gene.tsv", "w") as t2gout: 
        for gene, transcript in zip(gene_name_list, transcript_list): 
            t2gout.write("{}\t{}\n".format(transcript, gene))


filter_gff(input_file, output_file, list_of_taggs_to_keep)


