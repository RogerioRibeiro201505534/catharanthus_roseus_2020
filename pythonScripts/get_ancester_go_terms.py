#python script to generate a gene to geneset table 
#(with anceestor go terms with is_a relation only)
from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag
import csv
import os

file = "dammit_annotation_with_go.tsv" 


def make_gene2gs_map(file, ontology): 
  godag = GODag(obo_file= "go_ontology_basic_06_01.obo", optional_attrs={'relationship'})
  gene2gs_map = {}
  gene_processed = 0
  with open(file) as csv_file:
    dammit_annot = csv.reader(csv_file, delimiter='\t')
    next(dammit_annot)
    for row in dammit_annot: 
      gene_processed += 1 
      if row[ontology] != "":
        gene_id = row[0]
        go_terms = row[ontology].split(";")
        for term in go_terms:
          gosubdag_r0 = GoSubDag([term], godag, prt=None, relationships={'part_of'})
          try:
            go_terms_to_add = gosubdag_r0.rcntobj.go2parents[term]
          except KeyError: 
            print("The following go term was not processed {}".format(term))
            continue
          go_terms_to_add.add(term)
          for go_term_to_add in go_terms_to_add: 
            if gene_id not in gene2gs_map.keys():
              gene2gs_map[gene_id] = [go_term_to_add]
            else: 
              if go_term_to_add not in gene2gs_map[gene_id]: 
                gene2gs_map[gene_id].append(go_term_to_add)
      if gene_processed%1000 == 0: 
        print(gene_processed)
  return gene2gs_map

def make_gene2gs_table(gene2gs_map, ontology): 
  outfile = "gene_set_{}_with_ancestors.txt".format(ontology)
  with open(outfile, "w") as fh:  
    fh.write("gene\tgeset\n")
    for key, values in gene2gs_map.items():
      for value in values:
        row_to_write = "{}\t{}\n".format(key, value)
        fh.write(row_to_write)



          

  
gene2gs_map = make_gene2gs_map(file, ontology = 9)
make_gene2gs_table(gene2gs_map, "BP")
gene2gs_map = make_gene2gs_map(file, ontology = 11)
make_gene2gs_table(gene2gs_map, "MF")
gene2gs_map = make_gene2gs_map(file, ontology = 13)
make_gene2gs_table(gene2gs_map, "CC")



