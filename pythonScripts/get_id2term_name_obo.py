import csv 

input_file = "goslim_plant.obo"

def get_id2term_name_obo(input_file): 
    id2termname_map = {}
    with open(input_file, "r") as fh: 
        obo = csv.reader(fh, delimiter = ":")
        for i in range(20): 
            next(obo)
        for line in obo: 
            if line == []: 
                continue
            if line[0] == "[Term]": 
                flag = True
                continue
            if line[0] == "id" and flag == True: 
                new_id = "GO:{}".format(line[2])
                continue
            if line[0] == "name" and flag == True:
                term_name = line[1][1:] 
                id2termname_map[new_id] = term_name
                flag = False
                new_id = ""
    return id2termname_map
            
                


#get_id2term_name_obo(input_file, output_file)