# PseudoChecker - V1.0 - Core
# Imports all the necessary dependencies and external python scripts to run PseudoChecker
import os, sys, time, re, subprocess
from pseudochecker_analysis import *
from pseudochecker_utils import *
from Bio import SeqIO, pairwise2
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from operator import itemgetter
import argparse
from multiprocessing.dummy import Pool as ThreadPool   
from contextlib import closing
import time
import dill
import json 

def get_local_time():
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    print(current_time)

def trim_utr(exonscat, cds):
    '''Determines the UTR from the 1st and last exon (if exist)'''
    
    term = None
    global exon1_no_utr_start 
    global last_exon_no_utr_end
    for i in range(len(exonscat) + 1):
        if cds == exonscat[i : i + len(cds)]:
            exon1_no_utr_start = i
            break
    exonscat = exonscat[exon1_no_utr_start:]
    
    for i in range(len(exonscat) + 1):
        if exonscat[:i] == cds:
            
            term = i
            break
    '''the following asserts cds corresponds to the exons and specifies the error in the input data'''
    try:  
        assert term != None
    except AssertionError as msg:

        print("Please review your input. There are discrepancies between exons and cds that are not caused by UTR:")

        for pos, i in enumerate(exonscat):
                if i != cds[pos]:
                    print('Discrepancy at pos ' + str(pos) + ' of the CDS')

        sys.exit()

    if numtotalexons > 1:
        last_exon_no_utr_end = len(lastexon) - (len(exonscat) - term)
    else:
        last_exon_no_utr_end = term

def reverse_complement(seq):

    record = Seq(seq)
    return str(record.reverse_complement())



class GenomicTarget():

    def __init__(self, i):

        self.genomicregion = str(i.seq).upper()
        self.genomicregion = self.genomicregion.replace('-', '')
        self.pred_cds = ''
        self.curr_exon = 0 
        self.speciesname = str(i.id)
        self.dic_exons = ({})  # Main data structure, Key = exon ID (number), Values = [predexonwogaps,startofalign,endofalign,splicesite3,splicesite5,identity,nogapspredexon]
        self.dic_exon_cds_position = ({})  # Records the coordinates occupied by each exon in the predicted (CDS), Key = exon ID, Values = (inpos,finalpos)
        self.nr = 0
        self.pstop = 0
        self.pos = 0  # accumulates the last value of (endofalign) when partioning the entire genomic target regionzsa
        self.cdspos = 0  # accumulates the last value in the cds length occupied by the last exon
        self.preliminar_gaps_size = set()
        self.id = i.id
        self.exon_alns = {}
        

    def heuristicsearch(self, exon):
        '''Returns the most likely region in the genomic target for performing deterministic exon alignments'''
        hits = []

        def searchmatches(self, exon):
            for i in range(len(self.genomicregion) - len(exon) - 1):
                genomicregionwindow = self.genomicregion[i : i + len(exon)]
                matches = 0
                mismatches = 0
                for j in range(len(genomicregionwindow)):
                    if genomicregionwindow[j] == exon[j]:
                        matches += 1
                    else:
                        mismatches += 1
                    if mismatches > 2:
                        break
                if mismatches <= 2:
                    hits.append([matches, mismatches, i, i + len(exon) - 1])

        p = 18
        if len(exon) < p:
            p = len(exon)
        for q in range(p, 0, -1):
            if len(hits) >= 1:
                break
            wordsize = q
            for i in range(len(exon) - wordsize + 1):
                if len(hits) >= 1:
                    break
                else:
                    searchmatches(self, exon[i : i + wordsize])
        
        hits = sorted(hits, reverse=True)
        hits2 = []
        
        if len(hits) == 0:
            if args.verbose:
                print('No seeded alignment detected.')
            return None
        else:
            maxmatches = hits[0][0]
            for i in hits:
                if i[0] == maxmatches:
                    hits2.append(i)
            hits2 = sorted(hits2, key=itemgetter(2))
            return [hits2[0][2] - int(5 * len(exon)), hits2[0][3] + int(5 * len(exon))]
            

    def check_gaps_length(self, flanked_predexon, flanked_refexon):
        '''Captures the size of preliminar gaps'''
        for j in re.finditer('[-]+', flanked_predexon):
            self.preliminar_gaps_size.add(j.end() - j.start())
        for j in re.finditer('[-]+', flanked_refexon):
            self.preliminar_gaps_size.add(j.end() - j.start())

    

        

    
    def align_exon_to_genomic(self, exon, recovering_exon = False, flanking_exons=None):
        '''Alignment of each reference exon against the target genomic region'''
        usingheuristics = False
        if args.verbose:
            print('Exon ' + str(self.curr_exon) + ' of ' + str(self.id) + '\n')

        utr_alignment=False

        if utr_status == 0:  # If UTR's are not flanked
            
            if self.curr_exon == 1:
                if  len(exon[exon1_no_utr_start:]) > 12:
                    exon = exon[exon1_no_utr_start:]  # Trims 5' UTR
                else:
                    utr_alignment = True
                    
            if self.curr_exon == numtotalexons:
                if len(exon[:last_exon_no_utr_end]) > 12:  # Trims 3' UTR
                    exon = exon[:last_exon_no_utr_end]
                else:
                    utr_alignment = True

        if len(self.genomicregion) == 0:
            if args.verbose:
                print('No more alignments possible.')
            if self.speciesname not in partialsequences:
                partialsequences[self.speciesname] = [self.curr_exon]
            else:
                partialsequences[self.speciesname] = partialsequences[self.speciesname] + [self.curr_exon]
            return 0
        

        #print('Using heuristcs:' +  str(args.heuristics))

        if len(exon) * len(self.genomicregion) >= 200000 and recovering_exon == False and args.heuristics:
        #    print('Using Heuristics!')
            usingheuristics = True
            limits = self.heuristicsearch(exon)
            if limits == None:
                if self.speciesname not in partialsequences:
                    partialsequences[self.speciesname] = [self.curr_exon]
                else:
                    partialsequences[self.speciesname] = partialsequences[self.speciesname] + [self.curr_exon]
                return 0
            lowerbound = limits[0]
            upperbound = limits[1]
            if lowerbound <= 0:
                lowerbound = 0
            if upperbound >= len(self.genomicregion) - 1:
                upperbound = len(self.genomicregion) - 1
            genomicregionshort = self.genomicregion[lowerbound:upperbound]

            if args.debug:
                with open(mainpath + 'genomic_region_short.fa', 'at') as out:
                    new_record = SeqRecord(Seq(genomicregionshort ), id=self.speciesname+ '_' + str(self.curr_exon), description=self.speciesname)
                    SeqIO.write(new_record, out, 'fasta')

        if recovering_exon:

            #create short genomic region between already aligned exons to find missing exon
            upstream_pos, downstream_pos = get_exon_flanking_coordinates(self, flanking_exons[0], flanking_exons[1])

            if upstream_pos and downstream_pos:

                genomicregionshort = self.genomicregion[upstream_pos:downstream_pos]
            
            else: 
                if args.verbose:
                    print('No previously aligned exons found, not possible to continue search for missing exon!')
                return None

                
        global match
        global mismatch
        global gap
        global extend

        if bestfit_status == 0:  # Pre selected scoring scheme
            if args.verbose:
                print('Pre selected scoring scheme')
            if usingheuristics is False:
                if needle:
                    alignment = needle_alignment_emboss(self.genomicregion, exon)
                else:
                    alignment = pairwise2.align.globalms(
                        self.genomicregion,
                        exon,
                        match,
                        mismatch,
                        gap,
                        extend,
                        one_alignment_only=True,
                        penalize_end_gaps=(True, False),
                    )  # Global alignment with free end gaps
            else:
                if needle:
                    alignment = needle_alignment_emboss(genomicregionshort, exon)
                else:
                    alignment = pairwise2.align.globalms(
                        genomicregionshort,
                        exon,
                        match,
                        mismatch,
                        gap,
                        extend,
                        one_alignment_only=True,
                        penalize_end_gaps=(True, False),
                    )  # Global alignment with free end gaps
           
            if needle:
                predexon = str(alignment[0].seq)  # String containing the aligned genomic target sequence
                refexon = str(alignment[1].seq)  # String containing the aligned reference exon sequence    

        
            else:   
                
                predexon = alignment[0][0]  # String containing the aligned genomic target sequence
                refexon = alignment[0][1]  # String containing the aligned reference exon sequence
            
            for i in range(len(refexon)):
                if refexon[i] != '-':
                    startofalign = i
                    break
            for i in range(len(refexon) - 1, -1, -1):
                if refexon[i] != '-':
                    endofalign = i + 1
                    break
                else:
                    endofalign = len(refexon)
            flanked_predexon = predexon[startofalign:endofalign]  # Flanks the right and left end gaps of the alignment
            flanked_refexon = refexon[startofalign:endofalign]  # Flanks the right and left end gaps of the alignment

            if utr_status == 0 and utr_alignment:  # If UTR's are not flanked
        
                if self.curr_exon == 1:

                    exon_noutr = exon[exon1_no_utr_start: ]  # Trims 5' UTR
                    
                    corrected_utr_start = correct_index_for_gaps(flanked_refexon, exon1_no_utr_start, type='end')

                    flanked_predexon = flanked_predexon[corrected_utr_start:]
                    flanked_refexon = flanked_refexon[corrected_utr_start:]
                    
                    exonoccupancy[self.curr_exon] = len(exon_noutr) / len(cds) 

                if self.curr_exon == numtotalexons:

                    exon_noutr = exon[:last_exon_no_utr_end]  # Trims 3' UTR

                    corrected_utr_end = correct_index_for_gaps(flanked_refexon, last_exon_no_utr_end, type='start')

                    flanked_predexon = flanked_predexon[:corrected_utr_end]
                    flanked_refexon = flanked_refexon[:corrected_utr_end]
                    
                    exonoccupancy[self.curr_exon] = len(exon_noutr) / len(cds) 

                

                
            else:
                exonoccupancy[self.curr_exon] = len(exon) / len(cds) 



            flanked_predexon_nogaps = flanked_predexon.replace('-', '')
            splicesite5 = predexon[endofalign : endofalign + 2]  # 5' splice site (right splice site)
            splicesite3 = predexon[startofalign - 2 : startofalign]  # 3' splice site (left splice site)
            ident = identity(flanked_predexon, flanked_refexon)

        if bestfit_status == 1:  # Tests three different alignment scoring schemes
            
            testedalignments = []

            matchreward = [1, 2, 3, 4, 5, 6, 7, 8]
            mismatchreward = [-1, -2, -3, -4, -5]
            c = 0
            opt = False
            
            for s in matchreward:
                if opt == True:
                    break
                for j in mismatchreward:
                    c += 1

                    if usingheuristics is False:
                        if needle:

                            alignment = needle_alignment_emboss(self.genomicregion, exon, s, j, self.id)
                            if alignment == None:
                                if args.verbose:
                                    print('alignmentNone')
                                break
                        else:
                            alignment = pairwise2.align.globalms(
                                self.genomicregion, exon, s, j, gap, extend, one_alignment_only=True, penalize_end_gaps=(True, False)
                            )  # Global alignment with free end gaps

                            if len(alignment) == 0:
                                break
                    else:
                        if needle:
                            alignment = needle_alignment_emboss(genomicregionshort, exon, s, j, self.id)
                            if alignment == None:
                                if args.verbose:
                                    print('alignmentNone')
                                break
                        else:
                            alignment = pairwise2.align.globalms(
                                genomicregionshort,
                                exon,
                                s,
                                j,
                                gap,
                                extend,
                                one_alignment_only=True,
                                penalize_end_gaps=(True, False),
                            )  # Global alignment with free end gaps 

                            if len(alignment) == 0:
                                break
                    
                    if needle:
                        predexon = str(alignment[0].seq)  # String containing the aligned genomic target sequence
                        refexon = str(alignment[1].seq)  # String containing the aligned reference exon sequence     

                       
                    else:   
                        predexon = alignment[0][0]  # String containing the aligned genomic target sequence
                        refexon = alignment[0][1]  # String containing the aligned reference exon sequence
                    
                    for i in range(len(refexon)):  # Flanks the left end gaps of the alignment
                        if refexon[i] != '-':
                            
                            startofalign = i
                            break
                    for i in range(len(refexon) - 1, -1, -1):  # Flanks the right end gaps of the alignment
                        if refexon[i] != '-':
                            endofalign = i + 1
                            break
                        else:
                            endofalign = len(refexon)
                    
                    flanked_predexon = predexon[
                        startofalign:endofalign
                    ]  # Flanks the right and left end gaps of the alignment
                    flanked_refexon = refexon[
                        startofalign:endofalign
                    ]  # Flanks the right and left end gaps of the alignment
                    
                   
                    alignment_identity = identity(flanked_predexon, flanked_refexon)


                    if utr_status == 0 and utr_alignment:  # If UTR's are not flanked
        
                        if self.curr_exon == 1:

                            exon_noutr = exon[exon1_no_utr_start: ]  # Trims 5' UTR
                            
                            corrected_utr_start = correct_index_for_gaps(flanked_refexon, exon1_no_utr_start, type='end')

                            flanked_predexon = flanked_predexon[corrected_utr_start:]
                            flanked_refexon = flanked_refexon[corrected_utr_start:]
                            
                            exonoccupancy[self.curr_exon] = len(exon_noutr) / len(cds) 

                        if self.curr_exon == numtotalexons:

                            exon_noutr = exon[:last_exon_no_utr_end]  # Trims 3' UTR

                            corrected_utr_end = correct_index_for_gaps(flanked_refexon, last_exon_no_utr_end, type='start')

                            flanked_predexon = flanked_predexon[:corrected_utr_end]
                            flanked_refexon = flanked_refexon[:corrected_utr_end]
                            
                            exonoccupancy[self.curr_exon] = len(exon_noutr) / len(cds) 

                       

                        
                    else:
                        exonoccupancy[self.curr_exon] = len(exon) / len(cds) 

                    
                   

                    if self.curr_exon != numtotalexons:
                        splicesite5 = predexon[endofalign : endofalign + 2]  # 5' splice site (right splice site)
                    else:
                        splicesite5='NA'
                    if self.curr_exon != 1:
                        splicesite3 = predexon[startofalign - 2 : startofalign]  # 3' splice site (left splice site)
                    else:
                        splicesite3='NA'
                    gaps = int(flanked_predexon.count('-') + flanked_refexon.count('-'))
                    reliable = None
                    if (flanked_predexon.count('-') - flanked_refexon.count('-')) % 3 == 0:
                        reliable = True
                    else:
                        reliable = False
                  
                    testedalignments.append(
                        [
                            alignment_identity,
                            flanked_predexon,
                            flanked_refexon,
                            startofalign,
                            endofalign,
                            splicesite5,
                            splicesite3,
                            gaps,
                            reliable,
                        ]
                    )
                    if numtotalexons > 1:
                        if self.curr_exon == 1:
                            if splicesite5 == 'GT' and (flanked_predexon.count('-') - flanked_refexon.count('-')) % 3 == 0:
                                opt = True
                                break
                        if self.curr_exon == numtotalexons:
                            if splicesite3 == 'AG' and (flanked_predexon.count('-') - flanked_refexon.count('-')) % 3 == 0:
                                opt = True
                                break
                        else:
                            if (
                                splicesite5 == 'GT'
                                and splicesite3 == 'AG'
                                and (flanked_predexon.count('-') - flanked_refexon.count('-')) % 3 == 0
                            ):
                                opt = True
                                break
                    else:
                        if (flanked_predexon.count('-') - flanked_refexon.count('-')) % 3 == 0:
                            opt = True
                            break
            splicefunctionalalignments = []
            
            if args.debug:
                jsonStr = json.dumps(testedalignments)

                with open(mainpath + "tested_alignments.json", "at") as outfile:
                    outfile.write(jsonStr)

            if numtotalexons > 1:
                if self.curr_exon == 1:
                    for i in testedalignments:
                        if i[5] == 'GT' or i[5] == 'GC':
                            splicefunctionalalignments.append(i)
                elif self.curr_exon == numtotalexons:
                    for i in testedalignments:
                        if i[6] == 'AG':
                            splicefunctionalalignments.append(i)
                else:
                    for i in testedalignments:
                        if (i[5] == 'GT' or i[5] == 'GC') and i[6] == 'AG':
                            splicefunctionalalignments.append(i)

            if len(splicefunctionalalignments) > 0:  # At least an alignment has functional splice sites
                if args.verbose:
                    print('End splicefunctionalalignments')

                splicefunctionalreliable = []
                for i in splicefunctionalalignments:
                    if i[8] == True:
                        splicefunctionalreliable.append(i)
                splicefunctionalreliable = sorted(
                    splicefunctionalreliable, key=itemgetter(7, 0)
                )  # Sorts by the number of gaps
               
                if len(splicefunctionalreliable) > 0:
                    flanked_predexon = splicefunctionalreliable[0][1]
                    flanked_refexon = splicefunctionalreliable[0][2]
                    flanked_predexon_nogaps = flanked_predexon.replace('-', '')
                    startofalign = splicefunctionalreliable[0][3]
                    endofalign = splicefunctionalreliable[0][4]
                    splicesite5 = splicefunctionalreliable[0][5]
                    splicesite3 = splicefunctionalreliable[0][6]
                    ident = alignment_identity
                else:
                    splicefunctionalalignments = sorted(
                        splicefunctionalalignments, key=itemgetter(7, 0)
                    )  # Sorts by the number of gaps
                    flanked_predexon = splicefunctionalalignments[0][1]
                    flanked_refexon = splicefunctionalalignments[0][2]
                    flanked_predexon_nogaps = flanked_predexon.replace('-', '')
                    startofalign = splicefunctionalalignments[0][3]
                    endofalign = splicefunctionalalignments[0][4]
                    splicesite5 = splicefunctionalalignments[0][5]
                    splicesite3 = splicefunctionalalignments[0][6]
                    ident = alignment_identity

            elif (
                len(splicefunctionalalignments) == 0 or numtotalexons == 1
            ):  # No alignment with functional splice sites or single exon gene
                
                nosplicingsitesbutnoframeshifts = []
                #assert len(testedalignments) > 0
                for i in testedalignments:
                    if i[8] == True:
                        nosplicingsitesbutnoframeshifts.append(i)
                if len(nosplicingsitesbutnoframeshifts) > 0:
                    
                    nosplicingsitesbutnoframeshifts = sorted(
                        nosplicingsitesbutnoframeshifts, key=itemgetter(7, 0)
                    )  # Sorts by the number of gaps
                    
                    flanked_predexon = nosplicingsitesbutnoframeshifts[0][1]
                    flanked_refexon = nosplicingsitesbutnoframeshifts[0][2]
                    flanked_predexon_nogaps = flanked_predexon.replace('-', '')
                    startofalign = nosplicingsitesbutnoframeshifts[0][3]
                    endofalign = nosplicingsitesbutnoframeshifts[0][4]
                    splicesite5 = nosplicingsitesbutnoframeshifts[0][5]
                    splicesite3 = nosplicingsitesbutnoframeshifts[0][6]
                    ident = alignment_identity
                elif len(testedalignments) > 0:
                    testedalignments = sorted(testedalignments, key=itemgetter(7, 0))  # Sorts by the number of gaps
                    flanked_predexon = testedalignments[0][1]
                    flanked_refexon = testedalignments[0][2]
                    flanked_predexon_nogaps = flanked_predexon.replace('-', '')
                    startofalign = testedalignments[0][3]
                    endofalign = testedalignments[0][4]
                    splicesite5 = testedalignments[0][5]
                    splicesite3 = testedalignments[0][6]
                    ident = alignment_identity
                else:
                    return(None)

        '''Finds downstream stop codon in the last exon'''
        r = 0
        if self.curr_exon == numtotalexons and findalternativestop == 1 and numtotalexons > 1:
            for i in range(0, len(flanked_predexon), 3):
                codontargetexon = flanked_predexon[len(flanked_predexon) - i - 3 : len(flanked_predexon) - i]
                codonrefexon = flanked_refexon[len(flanked_refexon) - i - 3 : len(flanked_refexon) - i]
                if codonrefexon.count('-') % 3 != 0 or codontargetexon.count('-') % 3 != 0:
                    r = 1
                    break
                if codontargetexon in stopcodons:
                    r = 1
                    break
            if r == 0:
                nogaps = flanked_predexon.count('-')
                for i in range(0, 15, 3):
                    codontarget = self.genomicregion[endofalign + nogaps + i : endofalign + nogaps + i + 3]
                    
                    if codontarget in stopcodons:
                        
                        flanked_predexon = (
                            flanked_predexon + self.genomicregion[endofalign + nogaps : endofalign + nogaps + i + 3]
                        )
                        flanked_refexon = flanked_refexon + '#' * (len(flanked_predexon) - len(flanked_refexon))
                        flanked_predexon_nogaps = flanked_predexon.replace('-', '')
                        break

        if usingheuristics == True:
            startofalign += lowerbound
            endofalign += lowerbound
        if ident < min_exon_ident:  # Avoids unspecific alignments
            if args.verbose:
                print( 
                    'Alignment Identity: '
                    + str(ident)
                    + '. Identity below the defined threshold by the user. ('
                    + str(min_exon_ident)
                    + '%)'
                )
                print('This exon either is lost in the target species or very eroded so that any similarity is destroyed.\n')
            if self.speciesname not in partialsequences:
                partialsequences[self.speciesname] = [self.curr_exon]
            else:
                partialsequences[self.speciesname] = partialsequences[self.speciesname] + [self.curr_exon]
            
        else:
            if len(splicesite3) < 2:
                splicesite3 = 'XX'  # If splice site does not exist (precise intron deletion or genomic region extremities)
            if len(splicesite5) < 2:
                splicesite5 = 'XX'  # If splice site does not exist (precise intron deletion or genomic region extremities)
            gapspredexon = flanked_predexon.count('-')
            self.dic_exon_cds_position[self.curr_exon] = [self.cdspos, self.cdspos + len(flanked_predexon_nogaps) - 1]
            self.cdspos += len(flanked_predexon_nogaps)
            
            if args.verbose:
                print('Alignment identity (%): ' + str(ident))
                print(
                    '\nAlignment in the target genomic region: from '
                    + str(startofalign + self.pos)
                    + ' to '
                    + str(endofalign - 1 + self.pos)
                    + '\n'
                )
            self.pred_cds += flanked_predexon_nogaps
            self.dic_exons[int(self.curr_exon)] = [
                flanked_predexon_nogaps,
                startofalign + self.pos,
                endofalign + self.pos - 1,
                splicesite3,
                splicesite5,
                ident,
                gapspredexon,
            ]  # Stores information about each predicted exon of the target species - Main data structure
            self.check_gaps_length(flanked_predexon, flanked_refexon)  # Checks for preliminar frameshift gaps
            
                
            self.exon_alns[self.curr_exon] = [flanked_predexon,
                flanked_refexon,
                self.curr_exon,
                splicesite3,
                splicesite5,
                numtotalexons,
                ident,
                startofalign + self.pos,
                endofalign - 1 + self.pos]
                
            self.genomicregion = self.genomicregion[
                endofalign:
            ]  # Forces the following exon to align in a region immediately after the location where the previous exon aligned (increases speed and memory saving)
            self.pos = self.pos + endofalign  # Accumulates the last end of align position
            self.exons_file.write('>Exon_' + str(self.curr_exon) + '\n' + flanked_predexon_nogaps + '\n')
            
            
    def process_exons(self, list_exs_cds):
        for e in list_exs_cds:  # Aligns each exon against the genomic region
            if e.id != 'CDS':
                self.curr_exon += 1
                exon = str(e.seq).upper()
                exon = exon.replace('-', '')
                self.align_exon_to_genomic(exon)
                
    


def identity(flanked_predexon, flanked_refexon):
    '''Returns the alignment identity'''
    identity = 0
    for i in range(len(flanked_refexon)):
        if flanked_predexon[i] == flanked_refexon[i]:
            identity += 1
    identity = identity / len(flanked_predexon) * 100
    identity = round(identity, 2)
    return float(identity)



def main_loop_wrapper(g_target):
    global list_exs_cds
    g_target.process_exons(list_exs_cds) 





if __name__ == '__main__':
    '''Input arguments for running an analysis'''

    ####################
    ###Get the options

    parser = argparse.ArgumentParser(description='Identify very divergent potentially non-homologous windows in a protein multiple sequence alignment.')

    parser.add_argument('--file_exs_cds', help='Relative path to file containing the reference species exons and CDS of the in-study gene')
    parser.add_argument("--file_genomic", help='Relative path to file containing the target genomic sequence per target species')
    parser.add_argument("--match", type=int, default=1, help='Match Reward')
    parser.add_argument("--mismatch", type=float, default=-3, help='Mismatch reward')
    parser.add_argument("--bestfit_status", type=int, default=1, help='Defines if the user selects the best fit alignment scoring scheme (0 = no, 1 = yes)')
    parser.add_argument("--file_additional_cds", default=None, help='Relative path to file containing the additional (optional) predetermined coding sequences for macse alignment')
    parser.add_argument("--find_alternative_stop", type=int, default=0, help='Defines if the user pretends to find a downstream final stop codon in the last exon (0 = no, 1 = yes)')
    parser.add_argument("--min_exon_ident", type=int, default=65, help='Defines the minimum exon alignment identity for it to be considered as positive alignment')
    parser.add_argument("--utr_status", type=int, default=0, help='Defines if the reference species 1st and/or last exon are already UTR flanked (UTRs were removed) or not (0 = not flanked, 1 = flanked)')
    parser.add_argument("--fs_lr", type=str, default='17', help='MACSE frameshift cost for less reliable sequences')
    parser.add_argument("--fs_lr_term", type=str, default='10', help='MACSE terminal frameshift cost for less reliable sequences')
    parser.add_argument("--stop_lr", type=str, default='10', help='MACSE premature stop codon cost for less reliable sequences')
    parser.add_argument("--fs", type=str, default='100', help='MACSE frameshift cost for reliable sequences')
    parser.add_argument("--fs_term", type=str, default='10', help='MACSE terminal frameshift cost for reliable sequences')
    parser.add_argument("--stop", type=str, default='50', help='MACSE premature stop codon cost for reliable sequences')
    parser.add_argument("--gap_ext", type=str, default='1', help='MACSE gap extension cost')
    parser.add_argument("--gap_ext_term", type=str, default='0.9', help='MACSE terminal gap extension cost')
    parser.add_argument("--gap_op", type=str, default='7', help='MACSE gap opening cost')
    parser.add_argument("--gap_op_term", type=str, default='6.3', help='MACSE terminal gap opening cost')
    parser.add_argument("--analysis_name",  help='Name of the analysis (job title)')
    parser.add_argument("--main_path", help='Main path to the job folder')
    parser.add_argument("-c", "--cpu",type=int, default=1, help="number of threads")
    parser.add_argument("--skip_MACSE", dest='skip_MACSE', action='store_true')
    parser.add_argument("--recover_exons", dest='recover_exons', action='store_true')
    parser.add_argument("--needle", dest='needle', default=False, action='store_true', help='uses biopython implementation of emboss needle (faster but still in testing)')
    parser.add_argument("--debug", dest='debug', default=False, action='store_true', help='outputs genomicregionshort as fasta and tested alignments as json')
    parser.add_argument("--disable_heuristics", dest='heuristics',  default=True, action='store_false', help = 'Stops heuristics step from first selecting a likely genomic region to align each exon')
    parser.add_argument("--verbose", dest='verbose', default=False, action='store_true')
    #parser.add_argument("--PairwiseAligner", dest='PairwiseAligner', default=False, action='store_true', help='uses bioalign PairwiseAligner')
    parser.set_defaults(skip_MACSE=False)

    

    args = parser.parse_args()

    ####################
    match = args.match  # Match reward
    mismatch = args.mismatch  # Mismatch reward
    
    if args.needle:
        needle=args.needle
    else:
        needle = False
        
    bestfit_status = args.bestfit_status  # Defines if the user selects the best fit alignment scoring scheme (0 = no, 1 = yes)
    findalternativestop = args.find_alternative_stop  # Defines if the user pretends to find a downstream final stop codon in the last exon (0 = no, 1 = yes)
    file_exs_cds = args.file_exs_cds  # Relative path to file containing the reference species' exons and CDS of the in-study gene
    file_genomic = args.file_genomic  # Relative path to file containing the target genomic sequence per target species
    file_additional_cds = args.file_additional_cds  # Relative path to file containing the additional (optional) predetermined coding sequences
    min_exon_ident = args.min_exon_ident # Defines the minimum exon alignment identity for it to be considered as positive alignment
    utr_status = args.utr_status  # Defines if the reference species' 1st and/or last exon are already UTR flanked or not (0 = not flanked, 1 = flanked)
    fs_lr = args.fs_lr  # MACSE frameshift cost for less reliable sequences
    fs_lr_term = args.fs_lr_term  # MACSE terminal frameshift cost for less reliable sequences
    stop_lr = args.stop_lr # MACSE premature stop codon cost for less reliable sequences
    fs = args.fs  # MACSE frameshift cost for reliable sequences
    fs_term = args.fs_term # MACSE terminal frameshift cost for reliable sequences
    stop = args.stop  # MACSE stop codon cost for reliable sequences
    gap_ext = args.gap_ext  # MACSE gap extension cost
    gap_ext_term = args.gap_ext_term  # MACSE terminal gap extension cost
    gap_op = args.gap_op  # MACSE gap opening cost
    gap_op_term = args.gap_op_term  # MACSE terminal gap opening cost
    analysisname = args.analysis_name  # Name of the analysis (job title)

    
    if args.main_path[-1] != '/':
        mainpath = args.main_path + '/'
    else:
        mainpath = args.main_path  # Main path to the job folder
    if os.path.isdir(mainpath) == False:
        os.mkdir(mainpath)

    if min_exon_ident < 30:
        min_exon_ident = 30

    get_local_time()

    '''Global Variables and Small Loops'''
    list_genomic = list(SeqIO.parse(file_genomic, 'fasta'))  # Parses the FASTA file containing genomic target region/targets region
    list_exs_cds = list(SeqIO.parse(file_exs_cds, 'fasta'))  # Parses the FASTA file containing the reference exons and cds
    #need to change code to accept exon ids that are not just Exon_n
    if file_additional_cds:
        list_additional_cds = list(
            SeqIO.parse(file_additional_cds, 'fasta')
        )  # Parses the FASTA file containing the additional (optional) pre-determined coding sequences

        no_additional_cds = len(list_additional_cds)
    else:
        no_additional_cds = 0


    lastexon = str(list_exs_cds[len(list_exs_cds) - 2].seq)
    lastexon = lastexon.replace('-', '')  # Needed for trimming the UTR's

    '''UTR trimming variables'''
    exon1_no_utr_start = 0
    last_exon_no_utr_end = 0

    '''General variables'''
    gap = -8
    extend = -1

    nrglobal = 0
    numtotalexons = 0
    exonscat = ''  # Concatenation of all input reference exons for UTR trimming
    cds = str(list_exs_cds[len(list_exs_cds) - 1].seq)
    cds = cds.replace('-', '')  # Reference species coding sequence
    cdsposglobal = {}  # Data structure that saves the position of each predicted exon in the corresponding predicted cds

    macsereliablesequences = []
    macselessreliablesequences = []
    excludedsequences = []
    additionalcds = []
    partialsequences = {}
    exonoccupancy = {}
    stopcodons = {'TAA', 'TAG', 'TGA'}

    #analysis parameters will be saved to a file which will be needed to visualize results
    config_params = {"Match reward": match, "Mismatch reward": mismatch, "bestfit_status": bestfit_status, "findalternativestop": findalternativestop, 
                    "min_exon_ident": min_exon_ident, "utr_status": utr_status, "MACSE frameshift cost for less reliable sequences":fs_lr,
                    "MACSE terminal frameshift cost for less reliable sequences": fs_lr_term, "MACSE premature stop codon cost for less reliable sequences": stop_lr,
                    "MACSE frameshift cost for reliable sequences": fs, "MACSE terminal frameshift cost for reliable sequences": fs_term,
                    "MACSE stop codon cost for reliable sequences": stop, "MACSE gap extension cost": gap_ext, "MACSE terminal gap extension cost": gap_ext_term,
                    "MACSE gap opening cost": gap_op, "MACSE terminal gap opening cost": gap_op_term, "job_title": analysisname}

    json_config = json.dumps(config_params)

    with open(mainpath + "config_params.json", "w") as outfile:
        outfile.write(json_config)

    fasta_predicted_cds = open(
        mainpath + 'predictedcds.fasta', 'w'
    )  # FASTA file for saving all the predicted coding sequences during the analysis
    fasta_lessreliable = open(
        mainpath + 'lessreliable.fasta', 'w'
    )  # FASTA file for saving the MACSE less reliable predicted coding sequences
    fasta_reliable = open(
        mainpath + 'reliable.fasta', 'w'
    )  # FASTA file for saving the MACSE reliable predicted coding sequences
    partial_sequences = open(
        mainpath + 'partialsequences.txt', 'w'
    )  # FASTA file for saving the predicted partial sequences

    fasta_reliable.write(
        '>Reference_Species' + '\n' + cds + '\n'
    )  # Defines the reference coding sequence as a MACSE reliable sequence

    for i in list_exs_cds:  # Concatenates the exons for UTR trimming and captures the total no. of exons of the reference species gene
        if i.id != 'CDS':
            exon = str(i.seq)
            exon = exon.replace('-', '')
            exonscat += exon
            numtotalexons += 1

    if utr_status == 0:
        #if exonscat != cds:   #added this so it does not run trim_utr unecessarily which seems to cause an error when there is no utr
            
        trim_utr(exonscat, cds)  # Trims the UTR from the 1st and last exon if UTR status is set to non flanked UTRs
        #else:
        #    No_need_for_UTR_trimming = True
    #else:
    #    try:
    #        assert exonscat == cds
    #    except AssertionError:
    #        "Merged exons and CDS do not correspond! Are you sure the reference exons do not need UTR trimming?"
    #        sys.exit('''The program cannot continue until input is checked''')
    
    

    genomic_targets_list = []
    jobs = []

    '''Main Loop'''
    for t in list_genomic:  # Iterates over each target species genomic region

        i = GenomicTarget(t)
        
        i.exons_file = open(mainpath + i.speciesname + '_exons.fasta', 'w')
       
        if args.verbose:
            print('\nPredicting the CDS (coding sequence) in the target Species: ' + str(i.id))
            print('Aligning each exon against the genomic reference ...\n')
        
        genomic_targets_list.append(i)

  
    with closing(ThreadPool(processes=args.cpu)) as pool:

        pool.map(main_loop_wrapper, genomic_targets_list)
        pool.terminate()
    
    if args.recover_exons:

        def intermediate_func(g_target):
            #needed because pool.map does not accept more than one argument
            global list_exs_cds
            recover_exons(g_target, list_exs_cds)

        with closing(ThreadPool(processes=args.cpu)) as pool:

            pool.map(intermediate_func, genomic_targets_list)
            pool.terminate()


    
    if numtotalexons > 1:
        exon_pos_dict = get_pos_exons(list_exs_cds, cds)
    else:
        exon_pos_dict = {}
        exon_pos_dict[list_exs_cds[0].id] = [exon1_no_utr_start, last_exon_no_utr_end]

    exon_data = {}

        
    
    for i in genomic_targets_list:

        #needs reviewing or to go inside a function
        t_aligned_exons = [x for x in i.exon_alns.keys()]
        exonscat_withgaps = ''
        reference_exonscat_withgaps = ''

        exon_data[i.id] = i.exon_alns

        table_frameshifts = create_mutations_table(i, mainpath + 'frameshifts.csv')

        if table_frameshifts.empty == False:
            table_frameshifts_minexon = table_frameshifts[table_frameshifts.exon == table_frameshifts.exon.min()]
            table_frameshifts_min = table_frameshifts_minexon.start.values.tolist()[0]

        for e in list_exs_cds:
            if e.id != 'CDS':
                
                exon_n = int(e.id.split('_')[-1])
                
                if exon_n in t_aligned_exons :
                    exonscat_withgaps += i.exon_alns[exon_n][0]
                    reference_exonscat_withgaps += i.exon_alns[exon_n][1]
                else:
                    exon_len = (exon_pos_dict[e.id][1] - exon_pos_dict[e.id][0] + 1)
                    exonscat_withgaps += '-' * exon_len
                    reference_exonscat_withgaps += '*' * exon_len
        
        def reverse_modulo(n, x):

            mod = n%x
            rev_mod = x-mod

            return(rev_mod)
        
        
        #essential for correct search of stop codons after frameshift insertions
        if table_frameshifts.empty == False:
            for index, row in table_frameshifts.iterrows():
                if (row[3] % 3 != 0):
                    if row[2] == 'insertion':
                        cds=str(list_exs_cds[len(list_exs_cds) - 1].seq)
                        global_pos_frameshift = convert_exon_pos2global_pos(list_exs_cds, cds, row[1], row[4])
                        #like MACSE, it will add a gap to the insertion to keep the sequence in frame
                        exonscat_withgaps = exonscat_withgaps[0:global_pos_frameshift] + reverse_modulo(row[3], 3) * '-' + exonscat_withgaps[global_pos_frameshift:]
                        #to keep in frame the reference sequence also needs to be modified
                        reference_exonscat_withgaps = reference_exonscat_withgaps[0:global_pos_frameshift] + reverse_modulo(row[3], 3) * '-' + reference_exonscat_withgaps[global_pos_frameshift:]
        
        t_aligned_exons = [str(x) for x in t_aligned_exons]
        
        list_pstops = find_stop_codon2(exonscat_withgaps, reference_exonscat_withgaps)
        
        #count out insertions
            #now done in find_stop_codon2
        #list_pstops =  [x - reference_exonscat_withgaps[0:x].count('-') for x in list_pstops]

        with open(mainpath + 'stop_codons.csv', 'at') as stop_out:
            list_pstops_str = [str(i) for i in list_pstops]
            stop_out.write(i.id + ',' +  ",".join(list_pstops_str) + '\n')

        

        if len(list_pstops) > 0:
            truncated_percentage = truncated_codons(list_pstops, cds)
        else:
            truncated_percentage = 0

        shifted_percentage = shifted_codons(table_frameshifts, cds, list_exs_cds)

        if numtotalexons > 1:
            absent_percentage = absent_exons(i, list_exs_cds, cds)
        else:
            absent_percentage = 100-len(t_aligned_exons)*100

        i_pseudoindex = alt_pseudoindex(shifted_percentage, truncated_percentage, absent_percentage)

        with open(mainpath + 'species_pseudo_stats.csv', 'at') as handle:

            line = f'{args.analysis_name},{i.id},{i_pseudoindex},{shifted_percentage},{truncated_percentage},{absent_percentage}, {" ".join(t_aligned_exons)}'

            handle.write(line + '\n')
        
        

        if len(i.pred_cds) == 0:
            excludedsequences.append(i.speciesname)
        else:
            cdsposglobal[i.speciesname] = (i.dic_exon_cds_position)  # Adds the dictionary of the  analysed species into the general dictionary (covering the totality of the species) containing the pos occupied per exon in the CDS

            for p in i.preliminar_gaps_size:
                if p % 3 != 0:  # non multiple of 3 gap length
                    fasta_lessreliable.write(
                        '>' + i.speciesname + '\n' + i.pred_cds + '\n'
                    )  # Assings the predicted coding sequence as a MACSE less reliable sequence
                    i.nr = i.nr + 1
                    nrglobal += 1
                    if args.verbose:
                        print('Assigning this sequence as a MACSE less reliable coding sequence.\n')
                    macselessreliablesequences.append(i.speciesname)
                    break
            if i.nr == 0:
                for p in range(0, len(i.pred_cds), 3):
                    if p != 0 and p != len(i.pred_cds) - 3:
                        if i.pred_cds[p : p + 3] in stopcodons:
                            i.pstop = 1
                            nrglobal += 1
                            if args.verbose:
                                print('Assigning this sequence as a MACSE less reliable coding sequence.\n')
                            fasta_lessreliable.write(
                                '>' + i.speciesname + '\n' + i.pred_cds + '\n'
                            )  # Assings the predicted coding sequence as a MACSE less reliable sequence
                            macselessreliablesequences.append(i.speciesname)
                            break
                if i.pstop == 0:
                    fasta_reliable.write(
                        '>' + i.speciesname + '\n' + i.pred_cds + '\n'
                    )  # Assings the predicted coding sequence as a MACSE reliable sequence
                    macsereliablesequences.append(i.speciesname)
            fasta_predicted_cds.write('>' + i.speciesname + '\n' + i.pred_cds + '\n')
            
            
        i.exons_file.close()
    
    
    '''End of main loop'''
    general_stats_dict= general_stats(exon_data, list_exs_cds)
    general_stats_table = pd.DataFrame.from_dict(general_stats_dict).T.reset_index()
    print(general_stats_table)
    general_stats_table.columns = ['Target', "Avg. Exon Align. Pid", 'Avg. Exon Size', 'Min Exon Size', 'Larg. Exon size', 'Splice site integrity', 'Aligning exons']
    general_stats_table.to_csv(mainpath + 'general_stats.csv', sep=',', encoding='utf-8', index=False, header=True)
    
    get_local_time()


    if file_additional_cds:  # Adds the additional pre computed CDS to macse reliable sequences fasta file
        for i in list_additional_cds:
            additionalcds.append(str(i.id))
            fasta_reliable.write('>' + str(i.id) + '\n' + str(i.seq) + '\n')

    for i in partialsequences:
        partial_sequences.write(i + '\n')

    fasta_reliable.close()
    fasta_lessreliable.close()
    fasta_predicted_cds.close()
    partial_sequences.close()
    

    fasta_predicted_cds = open(mainpath + 'predictedcds.fasta', 'r')
    predicted_cdsfile = SeqIO.parse(fasta_predicted_cds, 'fasta')
    fasta_predicted_cds_translated = open(
        mainpath + 'predictedcdstranslated.fasta', 'w'
    )  # FASTA file for saving the predicted coding sequences during the analysis
    for i in predicted_cdsfile:  # Translates all the predicted coding sequences
        translatedseq = str(Seq(str(i.seq), generic_dna).translate())
        fasta_predicted_cds_translated.write('>' + str(i.id) + '\n')
        fasta_predicted_cds_translated.write(translatedseq + '\n')
    fasta_predicted_cds.close()
    fasta_predicted_cds_translated.close()

    if len(macselessreliablesequences) == 0 and len(macsereliablesequences) == 0:
        print('None of each reference exons aligned in any target species genomic region.')
        print('Analysis stopped.')
        noalignpage = open(mainpath + 'stopped.txt', 'w')
        noalignpage.close()
        sys.exit()

        #saving data for macse
    data_needed_for_macse = [cdsposglobal, partialsequences, excludedsequences, no_additional_cds, numtotalexons, exonoccupancy, additionalcds, nrglobal]

    with open(mainpath + 'first_step_data_dill.pkl', 'wb') as f:
        dill.dump(data_needed_for_macse, f)

    list_exs = [e for e in list_exs_cds if e.id != 'CDS']
    
    #save exon alns to json file
    exon_alns_2_json(genomic_targets_list, mainpath)
    #this will get the directory where the program is
    work_dir = os.getcwd() + '/'
    file_path = os.path.abspath(__file__)
    file_dir = file_path.replace('pseudochecker.py', '')

    

    if args.skip_MACSE == False:
        if args.verbose:
            print('Analysing the results from MACSE...\n')
            print('MACSE reliable sequences: ', macsereliablesequences)
            print('MACSE less reliable sequences: ', macselessreliablesequences)

        
        os.chdir(file_dir)

    
        if nrglobal >= 1:
            
            os.popen(
                'java -Xms1G -Xmx20G -jar '
                + 'macse_v2.03.jar -prog alignSequences -seq '
                + work_dir + '/' + mainpath
                + 'reliable.fasta -seq_lr '
                + work_dir + '/' + mainpath
                + 'lessreliable.fasta -out_NT '
                + work_dir + '/' + mainpath
                + 'macseanalysis_nt.fasta -out_AA '
                + work_dir + '/' + mainpath
                + 'macseanalysis_aa.fasta -fs '
                + str(fs)
                + ' -fs_term '
                + str(fs_term)
                + ' -fs_lr '
                + str(fs_lr)
                + ' -fs_lr_term '
                + str(fs_lr_term)
                + ' -stop_lr '
                + str(stop_lr)
                + ' -stop '
                + str(stop)
                + ' -gap_ext '
                + str(gap_ext)
                + ' -gap_ext_term '
                + str(gap_ext_term)
                + ' -gap_op '
                + str(gap_op)
                + ' -gap_op_term '
                + str(gap_op_term)
                + ' -local_realign_init 1 -local_realign_dec 1'
            ).read()
        else:
            os.popen(
                'java -Xms1G -Xmx20G -jar ' 
                + 'macse_v2.03.jar -prog alignSequences -seq '
                + work_dir + '/' + mainpath
                + 'reliable.fasta -out_NT '
                + work_dir + '/' + mainpath
                + 'macseanalysis_nt.fasta -out_AA '
                + work_dir + '/' + mainpath
                + 'macseanalysis_aa.fasta -fs '
                + str(fs)
                + ' -fs_term '
                + str(fs_term)
                + ' -fs_lr '
                + str(fs_lr)
                + ' -fs_lr_term '
                + str(fs_lr_term)
                + ' -stop_lr '
                + str(stop_lr)
                + ' -stop '
                + str(stop)
                + ' -gap_ext '
                + str(gap_ext)
                + ' -gap_ext_term '
                + str(gap_ext_term)
                + ' -gap_op '
                + str(gap_op)
                + ' -gap_op_term '
                + str(gap_op_term)
                + ' -local_realign_init 1 -local_realign_dec 1'
            ).read()


        os.chdir(work_dir)
        print('\nFinished MACSE analysis!\n' + os.getcwd())

       
        order_alignment(mainpath)  # Auxiliary Function
        
        macse_analysis(mainpath)  # Auxiliary Functionh
        
        
        
        

        if args.verbose:
            print(excludedsequences)
        excluded = excludedsequences[:]
        
       

        successful_confirmation_file = open(mainpath + 'success.txt', 'w')
        successful_confirmation_file.close()
        
    elif args.skip_MACSE:

        print('MACSE alignment will be skipped!')

    print('End of the analysis!')
    
    
    
    get_local_time()
    
    
    if numtotalexons == 1:
        single_exon_confirm = open(mainpath + 'single_exon.txt', 'w')
        single_exon_confirm.close()

