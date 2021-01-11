import csv 
annotation_file = "dammit_annotation_final.tsv"
output_file = "annotations.gaf"


def convert_gafe(annotation_file, output_file): 
    with open(annotation_file, "r") as fh: 
        dammit_annot = csv.reader(fh, delimiter = "\t")
        next(dammit_annot)
        for row in dammit_annot:
            #For BP
            if row[9] != "": 
                go_terms = row[9].split(";")
                for go_term in go_terms: 
                    gaf_row = "Uniprotkb\t{}\t482159\toptional\t{}\tpubmed\tunknown\toptional\t{}\toptional\toptional\tprotein\tunknown\t482159\tWB\toptional\toptional\n".format(row[0], go_term, "BP")
                    fh = open(output_file,"a")
                    fh.write(gaf_row)
                    fh.close()
            if row[11] != "": 
                go_terms = row[11].split(";")
                for go_term in go_terms: 
                    gaf_row = "Uniprotkb\t{}\t482159\toptional\t{}\tpubmed\tunknown\toptional\t{}\toptional\toptional\tprotein\tunknown\t482159\tWB\toptional\toptional\n".format(row[0], go_term, "MF")
                    fh = open(output_file,"a")
                    fh.write(gaf_row)
                    fh.close()
            if row[13] != "": 
                go_terms = row[13].split(";")
                for go_term in go_terms: 
                    gaf_row = "Uniprotkb\t{}\t482159\toptional\t{}\tpubmed\tunknown\toptional\t{}\toptional\toptional\tprotein\tunknown\t482159\tWB\toptional\toptional\n".format(row[0], go_term, "CC")
                    fh = open(output_file,"a")
                    fh.write(gaf_row)
                    fh.close()

convert_gafe(annotation_file, output_file)
