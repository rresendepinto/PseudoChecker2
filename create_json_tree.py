import sys
from ete3 import Tree
import argparse as ap
import json
from Bio import SeqIO
import re

def intersect_strings(cds, utr_exon):
    #works only for first exon
    clean_exon = utr_exon

    while clean_exon != cds[0:len(clean_exon)] and len(clean_exon) >= 6:

        clean_exon = clean_exon[1:]

    if clean_exon != cds[0:len(clean_exon)]:

        return None

    else:

        return clean_exon

def intersect_strings_intronless(cds, utr_exon):

    clean_exon = utr_exon

    while cds != clean_exon[0:len(cds)]:

        clean_exon = clean_exon[1:]

    clean_exon = clean_exon[0:len(cds)]

    if clean_exon != cds[0:len(clean_exon)]:

        return None

    else:

        return clean_exon

def get_pos_exons(list_exs_cds, cds):
    num_exons = len([e for e in list_exs_cds if e.id != 'CDS'])
    pos_exons = {}
    last_end = 0

    for e in list_exs_cds:  # Aligns each exon against the genomic region
        if e.id != 'CDS':

            if e.description.endswith('Exon_1') and e.description.endswith('Exon_' + str(num_exons)):
                #for intronless genes
                exon_seq = intersect_strings_intronless(cds, e.seq)

                if len(exon_seq) == 0:
                    return None
                else:
                    pos_exons[e.id] = (1, len(exon_seq))
                    last_end = len(exon_seq)

            elif e.description.endswith('Exon_1'):

                exon_seq = intersect_strings(cds, e.seq)

                if len(exon_seq) == 0:
                    return None
                else:
                    pos_exons[e.id] = (1, len(exon_seq))
                    last_end = len(exon_seq)

            elif e.description.endswith('Exon_' + str(num_exons)):
                #this excludes 3' utrs 
                pos_exons[e.id] = (last_end + 1, min(len(cds), last_end + len(e.seq)))

            else:

                pos_exons[e.id] = (last_end + 1, last_end + len(e.seq))
                last_end = last_end + len(e.seq)

    return(pos_exons)   

def convert_exon_pos2global_pos(list_exs_cds, exon, exon_pos):

    cds = str(list_exs_cds[len(list_exs_cds) - 1].seq)

    exon_pos_dict = get_pos_exons(list_exs_cds, cds)

    exon = 'Exon_' + str(exon)

    global_pos = exon_pos + exon_pos_dict[exon][0] -1 

    return(global_pos)

def get_exon_from_globalpos(list_exs_cds, global_pos):

    cds = str(list_exs_cds[len(list_exs_cds) - 1].seq)

    exon_pos_dict = get_pos_exons(list_exs_cds , cds)

    exon_determined = 0

    for e in exon_pos_dict.keys():
        
        if global_pos >= exon_pos_dict[e][0] and global_pos <= exon_pos_dict[e][1]:
            
            exon_determined = e 

    if exon_determined != 0:

        return exon_determined


def add_functionalState_mutations(json):
    global target_pseudoindex
    global target_missing_exons
    global target_frameshifts
    global target_premature_stops
    global list_exs_cds
    global pos_exons

    pseudoindex =  target_pseudoindex[json['name']]
    list_mutations = []

    if json['name'] in target_frameshifts:

        for frameshift in target_frameshifts[json['name']]:
            

            list_mutations.append({"exon": frameshift[0], "position": convert_exon_pos2global_pos(list_exs_cds, int(frameshift[0]), int(frameshift[3])), "type": frameshift[1], "length": frameshift[2]})

    

    for pstop in target_premature_stops[json['name']]:
        if pstop != '':
            list_mutations.append({"exon": get_exon_from_globalpos(list_exs_cds, int(pstop)), "position": pstop, "type": 'Premature Stop', "length": 3})

    for missing_exon in target_missing_exons[json['name']]:

        list_mutations.append({"exon": missing_exon, "position": pos_exons['Exon_' + str(missing_exon)][0], "type": 'Missing Exon', "length": pos_exons['Exon_' + str(missing_exon)][1] - pos_exons['Exon_' + str(missing_exon)][0]})

    json['Pseudoindex'] = pseudoindex
    json['mutations'] = list_mutations
    


def get_json(node):
    # Read ETE tag for duplication or speciation events
    if not hasattr(node, 'evoltype'):
        dup = ''
    elif node.evoltype == "S":
        dup = "N"
    elif node.evoltype == "D":
        dup = "Y"

    node.name = node.name.replace("'", '')
        
    json = { "name": node.name, 
             "display_label": node.name,
             "duplication": dup,
             "branch_length": str(node.dist),
             "common_name": node.name,
             "seq_length": 0,
             "type": "node" if node.children else "leaf",
             "uniprot_name": "Unknown",
             }
    
    global target_pseudoindex
    if node.name in target_pseudoindex.keys():
        add_functionalState_mutations(json)

    if node.children:
        json["children"] = []
        for ch in node.children:
            json["children"].append(get_json(ch))

    return json

def compare_mutations(mut1, mut2):

    if mut1['exon'] == mut2['exon']:
        if mut1['position'] == mut2['position']:
            if mut1['type'] == mut2['type']:
                if mut1['length'] == mut2['length']:
                    return True
                
    return False

def compare_list_mutations(l1,l2):
    #print('comparing mutations lists')
    common_mutations = []
    for mut1 in l1:
        for mut2 in l2:
            #print(mut1, mut2)
            if mut1 == mut2:
                common_mutations.append(mut1)

    return common_mutations


def get_pseudoindex(node): #will iterate through the node (in json) children and add a pseudoindex value and mutations to this node
    #nodes are json here
    Pseudoindex = [] 
    

    if 'children' in node.keys():
        # Start with mutations of the first child
        list_mutations = get_pseudoindex(node['children'][0])['mutations']
        for ch in node['children']:
            #Pseudoindex.append(get_pseudoindex(ch)['Pseudoindex'])
            ch_node = get_pseudoindex(ch)
            #print(ch_node['name'])
            #print(ch_node['Pseudoindex'])

            Pseudoindex.append(ch_node['Pseudoindex'])
            list_mutations = compare_list_mutations(list_mutations, ch_node['mutations'])
            #print(f'{ch_node["name"]} : {ch_node["mutations"]}')
        #print('mutations remaining:')
        #print(list_mutations)
        node['mutations'] = list_mutations

        if len(set(Pseudoindex)) == 1:
            node['Pseudoindex'] = Pseudoindex[0]
            return(node)
        else:
            node['Pseudoindex'] = (min([int(x) for x in Pseudoindex]))
            return node
        
        

    else:
        return node

    

def change_target(text):
    pattern= r".+:\d+-\d+_"
    modified_text = re.sub(pattern, '', text)
    return modified_text

if __name__ == '__main__':

    parser = ap.ArgumentParser(description='Create a json file with a phylogeny to be displayed in PseudoViz')

    parser.add_argument("-f", "--results_folder", help="results folder from pseudochecker analysis")
    parser.add_argument("-r", "--reference", help="reference file used for pseudochecker analysis")
    parser.add_argument("-t", "--tree")

    args = parser.parse_args()

    if args.results_folder.endswith('/'):
        results_folder = args.results_folder
    else:
        results_folder = args.results_folder + '/'


    target_pseudoindex = {}
    target_missing_exons = {}

    list_exs_cds = list(SeqIO.parse(args.reference, 'fasta'))  # Parses the FASTA file containing the reference exons and cds

    n_exons = len(list_exs_cds)

    with open(results_folder + 'species_pseudo_stats.csv', 'r') as input:

        for line in input:
            fields = line.split(',')
            target = change_target(fields[1])
            pseudoindex = fields[2]

            target_pseudoindex[target] = pseudoindex

            exons = [int(x.strip()) for x in fields[6].split()]

            target_missing_exons[target] = []
          
            for i in range(1,n_exons):

                if i not in exons:

                    target_missing_exons[target].append(i)

    target_premature_stops = {}

    with open(results_folder + 'stop_codons.csv', 'r') as input:

        for line in input:

            fields = line.split(',')
            target = change_target(fields[0])
            
            target_premature_stops[target] =  [x.strip() for x in fields[1::]]

    target_frameshifts = {}

    with open(results_folder + 'frameshifts.csv') as input:

        for line in input:

            fields = line.split(',')
            
            target = change_target(fields[0])
            
            if target in target_frameshifts.keys():

                target_frameshifts[target].append([x.strip() for x in fields[1::]])

            else:

                target_frameshifts[target] = [ [x.strip() for x in fields[1::]] ]

    #pos_exons = get_pos_exons(list_exs_cds)
    
    tree = Tree(args.tree, format=1)
    #print(tree.write(format=0))
    #print(target_pseudoindex.keys())
    tree.prune(target_pseudoindex.keys())

   
    # TreeWidget seems to fail with simple quotes
    json_tree = get_json(tree)
    #print(json_tree)
    print(str(get_pseudoindex(json_tree)).replace("'", '"'))