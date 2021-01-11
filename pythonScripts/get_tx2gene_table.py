import csv 

input_file = "transcript2gene_id_raw.csv" 
output_file = "transcript2gene_id.csv" 


with open(input_file, "r") as fh: 
    csv = csv.reader(fh)
    for line in csv: 
        string = line[0].replace("\"", "").split(";")
        tup = (string[0].replace("transcript_id", "").replace(" ", ""), string[1].replace("gene_id", "").replace(" ", ""))
        outstring = "{}\t{}\n".format(tup[0], tup[1])
        outfh = open(output_file, "a")
        outfh.write(outstring)
        outfh.close()