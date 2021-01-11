#use this script to get the subset of reads that did not aligned with featureCounts 

from Bio import SeqIO
import time
import csv

#input file 
input_file = "corefiles2fastqfiles.txt"
filter_class = "Unassigned_NoFeatures"

def get_class_list(core_file, filter_class): 
    interesting_read_set = set()
    other_read_set = set()
    with open(core_file, "r") as fh: 
        core_file_table = csv.reader(fh, delimiter = "\t")
        for row in core_file_table: 
            if row[1] == filter_class: 
                interesting_read_set.add(row[0])
            if row[1] != filter_class: 
                other_read_set.add(row[0])

    read_set = interesting_read_set - other_read_set
    print(len(read_set))
    return read_set

def print_sequences_to_file(fastq_file, read_set): 
    prefix = ""
    i = 0 
    while(fastq_file[i]) != ".": 
        prefix += fastq_file[i]
        i += 1 
    output_file = "{}_output_reads.fasta".format(prefix)
    
    for record in SeqIO.parse(fastq_file, "fastq"):
        if record.id in read_set: 
            with open(output_file, "a") as output:
                SeqIO.write(record, output, "fasta")



if __name__ == "__main__":
#read the input file and make a dictionary
    input_list = []
    inputfh = open(input_file, "r")
    input_table = csv.reader(inputfh, delimiter = "\t")
    for row in input_table:
        input_list.append((row[0], row[1]))
    inputfh.close()

    for CORE_file, fastq_file in input_list:  
        print("processing {} and {}".format(CORE_file, fastq_file))
        read_set = get_class_list(CORE_file, filter_class)
        print_sequences_to_file(fastq_file, read_set)


