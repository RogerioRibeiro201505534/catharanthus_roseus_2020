#Use this script to convert the gtf file to gff combatible with JBrowser
import csv 
import re 

swissprot_transcript_file = "uniprot_out.tsv"
transcript_to_uniprot = "swissprot_transcript_output.tsv"
gtf_file = "merged_transcriptome.gtf"
gff3_output = "merged_transcriptome_JBrows.gff3"

def make_swissprot_name_dictionary(swissprot_transcript_file):
    id2swissprotname_map = {} 
    with open(swissprot_transcript_file) as fh: 
        table = csv.reader(fh, delimiter = "\t")
        next(table)
        for t in table: 
            id2swissprotname_map[t[0]] = t[1]
    return id2swissprotname_map

def transcript2uniprot_map(transcript_to_uniprot): 
    transcript2uniprot_map = {}
    with open(transcript_to_uniprot) as fh: 
        table = csv.reader(fh, delimiter = "\t")
        for t in table: 
            swiss_id = t[1].split("|")[1]
            transcript2uniprot_map[t[0]] = swiss_id
    return transcript2uniprot_map


def conver_gtf_to_gff(gtf_file, gff3_output, id2swissprotname_map, trans_to_uniprot_map):
    with open(gtf_file) as fh: 
        input_gtf = csv.reader(fh, delimiter= "\t")
        next(input_gtf)
        next(input_gtf)
        i = 0 
        for line in input_gtf:
            i += 1  
            new_line = line[:]
            if line[2] == "transcript":
                transcript_id_re = re.search("transcript_id \"(.*)\"",new_line[8].split(";")[1])
                transcript_id = transcript_id_re.group(1)
                if transcript_id in trans_to_uniprot_map.keys(): 
                    swissprot_name = id2swissprotname_map[trans_to_uniprot_map[transcript_id]]
                else:
                    swissprot_name = "hypothetical protein"
                new_line[8] = "ID={};Name={};Alias={};Note={}".format(transcript_id, transcript_id, transcript_id, swissprot_name)

            if line[2] == "exon": 
                new_line[2] = "CDS"
                new_line[8] = "Parent={}".format(transcript_id)
                
            
            outputfh = open(gff3_output, "a")
            outputfh.write("\t".join(new_line))
            outputfh.write("\n")
            outputfh.close()
            if i%10000 == 0: 
                print(str(i) + " lines were processed")

 
#ID=CRO_119431;Name=CRO_119431;Alias=CRO_119431;Note=hypothetical%20protein


id2swissprotname_map = make_swissprot_name_dictionary(swissprot_transcript_file)
trans_to_uniprot_map = transcript2uniprot_map(transcript_to_uniprot)
conver_gtf_to_gff(gtf_file, gff3_output, id2swissprotname_map, trans_to_uniprot_map)