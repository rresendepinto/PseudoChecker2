#extra functions for pseudochecker
import os, sys, time, re, subprocess, shutil
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import Align
from operator import itemgetter, length_hint
import argparse
import multiprocessing as mp
from multiprocessing.dummy import Pool as ThreadPool   
from contextlib import closing
import time
import pandas as pd
import json


class Mutation():

    def __init__(self, start, end, type_mut):   

        self.start = start
        self.end = end
        self.type_mut = type_mut

    def get_mut_length(self):

        return(self.end + 1 - self.start)

def intersect_strings(cds, utr_exon):
    #works only for first exon
    clean_exon = utr_exon

    #possible issues with very short exonss
    while clean_exon != cds[0:len(clean_exon)]:

        clean_exon = clean_exon[1:]

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
            if e.description.endswith('Exon_1'):

                exon_seq = intersect_strings(cds, e.seq)
                
                if len(exon_seq) == 0 or exon_seq == None:
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

def find_stop_codon(seq): #deprecated
    stopcodons = {'TAA', 'TAG', 'TGA'}
    stop_presence=False
    stop_list = []
    for p in range(0, len(seq), 3):
        if p != 0 and p != len(seq) - 3:
            if seq[p : p + 3] in stopcodons:
                stop_presence = True
                stop_list.append(p + 1)

    return(stop_list)

def find_stop_codon2(seq1, refseq2):

    stopcodons = {'TAA', 'TAG', 'TGA'}
    
    stop_list = []
    #find where reference sequence actually starts, i.e. after any possible leading gaps
    start_of_ref = 0
    while refseq2[start_of_ref] == '-':
        start_of_ref += 1  
    #get indices of gaps in the reference sequence
    indices = [i for i, x in enumerate(refseq2) if x == "-" and i > start_of_ref]
    
    #now search the target sequence with gaps for stop codons
    for p in range(0, len(seq1), 3):

        if p != len(seq1) - 3:

            if seq1[p : p + 3] in stopcodons:   
                #add the stop codons to the list but subtract the number of gaps in the reference sequence before the stop codon
                stop_list.append(p + 1 - len([x for x in indices if x < p + 1])) 

    return(stop_list)
                       
        
def convert_exon_pos2global_pos(list_exs_cds, cds, exon, exon_pos):

    exon_pos_dict = get_pos_exons(list_exs_cds, cds)

    exon = 'Exon_' + str(exon)

    global_pos = exon_pos + exon_pos_dict[exon][0] -1 

    return(global_pos)

    

def check_exon_mutations(ref, target):

    mutations = []
    
    #gap opened
    gap = None

    #where is the gap
    gap_which = ''

    for i in range(len(ref)):
        #if difference between sequences
        
        if ref[i] != target[i]:
            
            if gap == None:

                if ref[i] == '-':

                    gap=[i, i]

                    gap_which = 'Reference'

                elif target[i] == '-':

                    gap=[i, i]

                    gap_which = 'Target'

            else:

                if gap_which == 'Reference':

                    if ref[i] == '-':
                        gap[1] = i


                    else:

                        x = Mutation(gap[0], gap[1], 'insertion')
                        mutations.append(x)                        

                        gap=None

                        gap_which = ''

                        if target[i] == '-':
                            
                            gap=[i, i]

                            gap_which = 'Target'

                elif gap_which == 'Target':

                    
                    if target[i] == '-':
                        gap[1] = i


                    else:

                        x = Mutation(gap[0], gap[1], 'deletion')
                        mutations.append(x)

                        gap=None

                        gap_which = ''

                        if ref[i] == '-':
                            
                            gap=[i, i]

                            gap_which = 'Reference'


        elif gap:

            if gap_which == 'Reference':

                x = Mutation(gap[0] + 1, gap[1] + 1, 'insertion')
                mutations.append(x)
                gap = None

            elif gap_which == 'Target':

                x = Mutation(gap[0] + 1, gap[1] + 1, 'deletion')
                mutations.append(x)
                gap = None

        

    
    return mutations
   

       


def create_mutations_table(g, path):
    #this function takes a genomic target class and returns a table with annotated mutation which will also be outputted to a file
    #it iterates over the exon_alns attribute of this class

    exon_dict = {}
    
    for exon in g.exon_alns.keys():
        
        record_ref = SeqRecord(Seq(g.exon_alns[exon][1].upper()),  id="Reference_species")
        
        record_target = SeqRecord(Seq(g.exon_alns[exon][0].upper()),  id=g.id)

        exon_mutations = check_exon_mutations(record_ref.seq, record_target.seq)
        
        exon_dict[exon] = exon_mutations

    list_frameshifts = []
    
    with open(path, 'at') as handle:

        for exon in exon_dict.keys():

            if exon_dict[exon]:

                for mut in exon_dict[exon]:

                    #if mut.get_mut_length % 3 != 0:

                    output = f'{g.id},{exon},{mut.type_mut},{mut.get_mut_length()},{mut.start},{mut.end}\n'
                    list_frameshifts.append([g.id,exon,mut.type_mut, mut.get_mut_length(), mut.start, mut.end])
                    handle.write(output)
    
    table_frameshifts = pd.DataFrame(list_frameshifts, columns=['id', 'exon', 'type','len', 'start', 'end'])
    return(table_frameshifts)
          
def add_stops_to_mutation_table(t, path , list_exs_cds):
    #this functions creates a table with stop codons in the job's main path
    
    t_aligned_exons = [x for x in t.exon_alns.keys()]
    
    #to find stop codons first we need a sequence with the exons joined but with gaps representing the missing exons
    #so the lack of an exon doesnt alter the reading frame like a frameshift would
    exonscat_withgaps = ''
    for e in list_exs_cds:
        if e.id != 'CDS':
            
            exon_n = int(e.description.split('_')[-1])
            
            if exon_n in t_aligned_exons :
                exonscat_withgaps += t.exon_alns[exon_n][0]
            else:
                exonscat_withgaps += '-' * len(e.seq)
    #returns a list of positions containing premature stop codons
    pstops = find_stop_codon(exonscat_withgaps)
    
    with open(path, 'at') as handle:

        for pstop in pstops:

            output = f'{t.id},{"stop_codon"},{pstop}\n'

            handle.write(output)



def shifted_codons(table_mut, cds, list_exs_cds):
    #outputs the percentage of the cds that is shifted    
    
    #how shifted the frame is
    frameshift_score = 0

    #length of shifted frame
    total_shifted = 0

    #save position of previous mutation
    previous_mut = 0



    for index, row in table_mut.iterrows():
       
        #get position of mutation in cds instead of in exon
        mut_pos = convert_exon_pos2global_pos(list_exs_cds, cds, row['exon'], row['start'])

        #if frame is shifted add the length since previous mutation or since start of cds
        if frameshift_score % 3 != 0:

            total_shifted += mut_pos - previous_mut

        #change the score according to shiftedness of frame
        if row['type'] == "insertion":

            frameshift_score += row['len']

        elif row['type'] == "deletion":

            frameshift_score -= row['len']

        previous_mut = mut_pos
       
        
    #at the end add the lenght of shifted cds since last mutation (if the frame is still shifted)
    if frameshift_score % 3 != 0:
        
        total_shifted += len(cds) - previous_mut
    
    return total_shifted/len(cds) * 100

def truncated_codons(pstop_list, cds):

    first_pstop = min(pstop_list)

    return (len(cds) - first_pstop ) / len(cds) * 100

def absent_exons(genomic_target, list_exs_cds, cds):
    #determines the percentage of cds lost by missing exons
    list_exs = [e for e in list_exs_cds if e.id != 'CDS']
    
    present_exons = [x for x in genomic_target.dic_exons.keys()]

    exon_pos_dict = get_pos_exons(list_exs, cds)

    total_present = 0
    print(exon_pos_dict)
    for e in present_exons:

        total_present += exon_pos_dict['Exon_' + str(e)][1] - (exon_pos_dict['Exon_' + str(e)][0] - 1 )
    print(total_present)
    return (1 - total_present/len(cds)) * 100

def alt_pseudoindex(shifted_codons, truncated_codons, absent_exons):
    #pseudoindex created from the first step of the software , the alignment of the exons
    #since the original pseudoindex is created from the MACSE alignment
    pseudoindex = 0

    #first lets check the percentage of shifted codons from frameshifts
    if shifted_codons < 10:
        pseudoindex = max(0, pseudoindex)
    elif shifted_codons > 10 and shifted_codons < 15:
        pseudoindex = max(1, pseudoindex)
    elif shifted_codons > 15 and shifted_codons < 20:
        pseudoindex = max(2, pseudoindex)
    elif shifted_codons > 20 and shifted_codons < 25:
        pseudoindex = max(3, pseudoindex)
    elif shifted_codons > 25 and shifted_codons < 30:
        pseudoindex = max(4, pseudoindex)
    elif shifted_codons > 30:
        pseudoindex = 5
    
    #now lets check the percentage of absent cds from missing exons
    if absent_exons < 10:
        pseudoindex = max(0, pseudoindex)
    elif absent_exons > 10 and absent_exons < 15:
        pseudoindex = max(1, pseudoindex)
    elif absent_exons > 15 and absent_exons < 20:
        absent_exons = max(2, pseudoindex)
    elif absent_exons > 20 and absent_exons < 25:
        pseudoindex = max(3, pseudoindex)
    elif absent_exons > 25 and absent_exons < 30:
        pseudoindex = max(4, pseudoindex)
    elif absent_exons > 30:
        pseudoindex = 5


    #now lets check the percentage of truncated cds
    if truncated_codons < 10:
        pseudoindex = max(0, pseudoindex)
    elif truncated_codons > 10 and truncated_codons < 15:
        pseudoindex = max(1, pseudoindex)
    elif truncated_codons > 15 and truncated_codons < 20:
        pseudoindex = max(2, pseudoindex)
    elif truncated_codons > 20 and truncated_codons < 25:
        pseudoindex = max(3, pseudoindex)
    elif truncated_codons > 25 and truncated_codons < 30:
        pseudoindex = max(4, pseudoindex)
    elif truncated_codons > 30:
        pseudoindex = 5

    return pseudoindex

def map_numbers_to_ref_exon(e_start, exon):
    n = int(e_start) 
    n_spaces = 0
    out_string = ''
    while n < (e_start + len(exon)):
        if (n_spaces == 0):
            out_string += '|' + str(n)
            n_spaces = 7
            n+=4
            
        else:
            out_string += ' '
            n_spaces -= 1
            n += 1

        

    return(out_string)






def exon_alns_2_json(genomic_targets_list, mainpath):

    targets_for_json = {}

    for t in genomic_targets_list:
    
        targets_for_json[t.id] = t.exon_alns

    jsonStr = json.dumps(targets_for_json)

    with open(mainpath + "exon_alns.json", "w") as outfile:
        outfile.write(jsonStr)


def needle_alignment_emboss(s1, s2, match, mismatch, target_id):
    import subprocess
    from Bio.Emboss.Applications import NeedleCommandline
    from Bio import AlignIO
    
    if match and mismatch:
        

        cline = NeedleCommandline(auto=True, stdout=True, gapopen=8, gapextend=1, datafile=os.path.join(os.path.dirname(__file__), f"scoring_matrices/scoring_matrix_{match}_{mismatch}.txt"))
    else:
        cline = NeedleCommandline(auto=True, stdout=True, gapopen=8, gapextend=1)
    
   
    cline.asequence = "asis:" + s1
    cline.bsequence = "asis:" + s2

    try:
        process = subprocess.Popen(str(cline), shell=True, stdout=subprocess.PIPE, universal_newlines=True)
    except:
        print("Failed running process")

        print(target_id)
        print('match: ' + str(match))
        print('mismatch: ' + str(mismatch))

        print('asequence: ' + str(s1))
        print('bsequence: ' + str(s2))

    return AlignIO.read(process.stdout, "emboss")
    
    

def pairwise_aligner(seq1, seq2, match, mismatch):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match = match
    aligner.mismatch = mismatch
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    alignments = aligner.align(seq1, seq2)
    return alignments

def get_flanking_exons(aligned_exons, exon):
   
    upstream = max((x for x in aligned_exons if x < exon), default=None)
    downstream = min((x for x in aligned_exons if x > exon), default=None)
    return upstream, downstream


def get_exon_flanking_coordinates(g_target, upstream_exon, downstream_exon):
    
    if upstream_exon and downstream_exon:
        start_pos = g_target.exon_alns[downstream_exon][-1]
        end_pos = g_target.exon_alns[upstream_exon][-2]

    elif upstream_exon:

        start_pos = g_target.exon_alns[downstream_exon][-1]
        end_pos = len(g_target.genomicregion)
    
    elif downstream_exon:

        start_pos = 0
        end_pos = g_target.exon_alns[upstream_exon][-2]

    else:

        return None, None


    return start_pos, end_pos
    
    

def recover_exons(g_target, list_exs_cds):

    for i in range(1, len(list_exs_cds) - 1):
        e = list_exs_cds[i]
        if e.id != 'CDS':
            t_aligned_exons = [x for x in g_target.exon_alns.keys()]
            if i not in t_aligned_exons:
                g_target.curr_exon = i
                exon = str(e.seq)
                exon = exon.replace('-', '')
                upstream, downstream = get_flanking_exons(t_aligned_exons, i)
                if upstream and downstream:
                    print(f'Finding exon {i} of {g_target.id} using exons {upstream} and {downstream}')
                    g_target.align_exon_to_genomic(exon, recovering_exon=True, flanking_exons = (upstream, downstream))
                    print(f'Success in aligning missing exon {i} of {g_target.id} !!!')
                else:
                    print(f'Was not able to successfully retrieve coordinates of flanking exons for exon {i} of {g_target.id}')

def correct_index_for_gaps(seq1, index, type):
    #this function will convert the index in a gapless version of seq1 that are meant to represent a subsequence
    #to the corresponding coords in a gap-containing version
    gapless_seq1 = seq1.replace('-', '')

    #end means subsequence is from start of seq1 to index, start means subseq is from index to end of seq
    assert type in ['start', 'end']

    if type =='end':
        seq2 = gapless_seq1[index:]

        for i, char in enumerate(seq1):
            if seq1[i:].replace('-', '') == seq2:
                return i

 
    else:
        seq2 = gapless_seq1[:index]
        
        for i, char in enumerate(seq1):
            if seq1[:i].replace('-', '') == seq2:
                return i
    
    


