import argparse as ap
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord
import sys
import re


def check_for_reverted_order(exon_dict, cds):

    pos_of_exons_in_cds = {}

    for e in exon_dict.keys():
        substr = str(exon_dict[e].seq)
        pos_of_exons_in_cds[exon_dict[e].id] = cds.find(substr)
    
    values_pos = [p for p in pos_of_exons_in_cds.values() if p > -1]
    
    if values_pos == sorted(values_pos, reverse=True):
        return True
    elif values_pos == sorted(values_pos):
        return False
    else:
        return None

        


def intersect_strings(cds, utr_exon):
    #works only for first exon
    clean_exon = utr_exon

    while clean_exon != cds[0:len(clean_exon)] and len(clean_exon) >= 6:

        clean_exon = clean_exon[1:]

    if clean_exon != cds[0:len(clean_exon)]:

        return None

    else:

        return clean_exon

def get_pos_exons(dict_exs_cds, cds):
    num_exons = len([e for e in list_exs_cds if e.id != 'CDS'])
    pos_exons = {}
    last_end = 0

    for e in list_exs_cds:  # Aligns each exon against the genomic region
        if e.id != 'CDS':
            if e.description.endswith('Exon_1'):

                exon_seq = intersect_strings(cds, e.seq)
                print('Exon_seq1: ' + exon_seq)

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

def correct_exons(exon_dict, starting_exon, ending_exon):
    #used to remove exons which do not make up part of the cds

    
    names = [key for key in exon_dict.keys()]
    exons_extract = names[starting_exon:ending_exon+1]
    print('names:' + ','.join(names))
    dict_subset = {key: exon_dict[key] for key in exons_extract}
    
    exon_n = 0
    dict_subset_new_names= {}

    for e in dict_subset.keys():
        print(e)
        dict_subset_new_names[names[exon_n]] = dict_subset[e]
        exon_n += 1

    return dict_subset_new_names

def correct_ending_exon(exon_dict, ending_exon):
    

    names = [key for key in exon_dict.keys()]
    print('Printing names:')
    exons_extract = names[:ending_exon + 1]
    print(exons_extract)
    print('correcting end exon')
    print('names:' + ','.join(names))
    dict_subset = {key: exon_dict[key] for key in exons_extract}
    
    exon_n = 0
    dict_subset_new_names= {}

    for e in dict_subset.keys():
        print(e)
        dict_subset_new_names[names[exon_n]] = dict_subset[e]
        exon_n += 1

    return dict_subset_new_names

if __name__ == '__main__':

    parser= ap.ArgumentParser(description='convert ensembl exons to pseudochecker input')

    parser.add_argument('-f', '--file', help='fasta file')
    parser.add_argument('-o', '--out', default=None, help='If filtered output wanted, specify if no stop is desired or if stop at end of sequence is OK')
    #parser.add_argument('-t', '--type_output', default='No_stop_at_end')

    args = parser.parse_args()

    file =  open(args.file, 'r')

    Exon_n = 1

    record_dict = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
    file.close()
    exon_dict = {}

    check_transcript = []

    for record in record_dict.keys():

        transcript = re.split('-|/ +/', record)[1]
         
        check_transcript.append(transcript)

        if(len(set(check_transcript))!=1):
            sys.exit("Exons/cds of different transcripts are present. Cannot be converted like this, please fix it!")
        
    for record in record_dict.keys():

        
        if record.startswith('exon'):
            
            exon_dict['Exon_' + record.split('-')[-1]] = record_dict[record]
            
        elif record.startswith('rna'):
            CDS=record_dict[record]
    
    exonscat = ''
        
    
    if 'Exon_1' in exon_dict:
        print('Key exists in the dictionary')
    else:
        print('Key does not exist in the dictionary')
    for i in range(1, len(exon_dict.keys()) + 1):
        
        exonscat += str(exon_dict['Exon_' + str(i)].seq)

    pos_start_cds = 0
    #print(exon_dict.keys())
    
    for i in range(len(exonscat)):
        
        if exonscat[i:i+len(CDS.seq)] == CDS.seq:
            
            pos_start_cds = i
            break

    for i in range(len(exonscat), -1, -1):
       
        if exonscat[i-len(CDS.seq):i] == CDS.seq:

            pos_end_cds = i
            break


    #account for index starting on 0
    pos_start_cds += 1 

    
    tmp_exon_dict = {}

    exon_dict_numbers = [int(x.split('_')[1]) for x in exon_dict.keys()]
    exon_dict_numbers.sort()

    for e in exon_dict_numbers:

        tmp_exon_dict['Exon_' + str(e)] = exon_dict['Exon_' + str(e)]

    exon_dict = tmp_exon_dict
    

    list_exs_cds = [exon_dict[e] for e in exon_dict.keys()]
    exons_pos_dict = get_pos_exons(list_exs_cds, CDS)
    
    print(exons_pos_dict)

    cur_startexon_n = 0

    exon_names = [key for key in exons_pos_dict.keys()]

    #remove exons before start of cds
    if pos_start_cds >= len(exon_dict['Exon_1'].seq):
        
        #for when the coding sequence does not start in the first exon
        
 
        while pos_start_cds > exons_pos_dict[exon_names[cur_startexon_n]][1]:
            cur_startexon_n += 1


    #remove exons after end of cds
    cur_endexon_n = 0
    while pos_end_cds > exons_pos_dict[exon_names[cur_endexon_n]][1]:
        
        cur_endexon_n += 1

    exon_dict = correct_exons(exon_dict, starting_exon=cur_startexon_n, ending_exon=cur_endexon_n)

    """
    if cur_exon_n < len(list_exs_cds) - 1:
        print(cur_exon_n)
        print('Need to correct end exon')
        exon_dict = correct_ending_exon(exon_dict, ending_exon=cur_exon_n)
    """
    
    

    with open(args.out, 'wt') as out:
        
        exon_n = 1

        while exon_n <= len(exon_dict.keys()):

            exon_name = 'Exon_' + str(exon_n)

            new_record = SeqRecord(seq= exon_dict[exon_name].seq, id= exon_name, description='') 

            SeqIO.write(new_record, out, 'fasta')

            exon_n += 1

        cds_record = SeqRecord(seq= CDS.seq, id='CDS', description='')         
        SeqIO.write(cds_record, out, 'fasta')










        
  