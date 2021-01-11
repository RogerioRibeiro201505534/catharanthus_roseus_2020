import csv

blastx_file = "BLASTxout_Catha2arath.txt"
blastn_file = "tBLASTnout_Arath2catha.txt"

def make_blast_dictionary(blast_file): 
    query2result_dictionary = {}
    with open(blast_file, "r") as fh: 
        blastresult = csv.reader(fh, delimiter = "\t")
        for line in blastresult:
            if line[0] not in query2result_dictionary.keys(): 
                query2result_dictionary[line[0]] = line
            else: 
                if float(query2result_dictionary[line[0]][11]) < float(line[11]): 
                    query2result_dictionary[line[0]] = line
    
    output = "{}_best_hit.txt".format(blast_file)
    print(output)
    with open(output, "w") as ofh: 
        for values in query2result_dictionary.values(): 
            ofh.write("\t".join(values)) 
            ofh.write("\n")


make_blast_dictionary(blastx_file) 
make_blast_dictionary(blastn_file)