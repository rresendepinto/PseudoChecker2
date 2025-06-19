from asyncio import FastChildWatcher
from distutils.command.build_scripts import first_line_re
import os, sys, time, re, subprocess
from pseudochecker_analysis import *
from pseudochecker_utils import *
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from operator import itemgetter
import argparse
import dill
import subprocess as sp

def prune_seqs(file, ids, out):

    with open(ids, 'r') as handle:
        first_line = handle.readline()
        print(first_line)
        try:
            assert first_line.strip() == 'Reference_Species'
        except AssertionError:
            sys.exit("Reference_Species must be in first line of ids to use for MACSE analysis")

    process_seqtk = f'seqtk subseq {file} {ids} > {out}' 
    process1 = subprocess.Popen(process_seqtk, stdout=subprocess.PIPE, shell=True)
    process1.wait()



def run_MACSE(mainpath, fs, fs_term, fs_lr, fs_lr_term, stop_lr, stop, gap_ext, gap_ext_term, gap_op, gap_op_term, nrglobal):

    if nrglobal >= 1:
        print('java -Xms1G -Xmx20G -jar '
            + 'macse_v2.03.jar -prog alignSequences -seq '
            + mainpath
            + 'reliable.fasta -seq_lr '
            + mainpath
            + 'lessreliable.fasta -out_NT '
            + mainpath
            + 'macseanalysis_nt.fasta -out_AA '
            + mainpath
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
            + ' -local_realign_init 1 -local_realign_dec 1')
        os.popen(
            'java -Xms1G -Xmx20G -jar '
            + 'macse_v2.03.jar -prog alignSequences -seq '
            + mainpath
            + 'reliable.fasta -seq_lr '
            + mainpath
            + 'lessreliable.fasta -out_NT '
            + mainpath
            + 'macseanalysis_nt.fasta -out_AA '
            + mainpath
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
        print('java -Xms1G -Xmx20G -jar ' 
            + 'macse_v2.03.jar -prog alignSequences -seq '
            + mainpath
            + 'reliable.fasta -out_NT '
            +  mainpath
            + 'macseanalysis_nt.fasta -out_AA '
            +  mainpath
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
            + ' -local_realign_init 1 -local_realign_dec 1')
        os.popen(
            'java -Xms1G -Xmx20G -jar ' 
            + 'macse_v2.03.jar -prog alignSequences -seq '
            + mainpath
            + 'reliable.fasta -out_NT '
            +  mainpath
            + 'macseanalysis_nt.fasta -out_AA '
            +  mainpath
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

        print('\nFinished MACSE analysis!\n' + os.getcwd())

if __name__ == '__main__':
    '''Input arguments for running an analysis'''

    ####################
    ###Get the options

    parser = argparse.ArgumentParser(description='Identify very divergent potentially non-homologous windows in a protein multiple sequence alignment.')

    #parser.add_argument('--file_exs_cds', help='Relative path to file containing the reference species exons and CDS of the in-study gene')
    #parser.add_argument("--file_genomic", help='Relative path to file containing the target genomic sequence per target species')
    #parser.add_argument("--file_additional_cds", default=None, help='Relative path to file containing the additional (optional) predetermined coding sequences')
    parser.add_argument("--match", type=int, default=1, help='Match Reward')
    parser.add_argument("--mismatch", type=float, default=-3, help='Mismatch reward')
    parser.add_argument("--bestfit_status", type=int, default=1, help='Defines if the user selects the best fit alignment scoring scheme (0 = no, 1 = yes)')
    parser.add_argument("--find_alternative_stop", type=int, default=0, help='Defines if the user pretends to find a downstream final stop codon in the last exon (0 = no, 1 = yes)')
    parser.add_argument("--utr_status", type=int, default=0, help='Defines if the reference species 1st and/or last exon are already UTR flanked or not (0 = not flanked, 1 = flanked)')
    parser.add_argument("--fs_lr", type=str, default='17', help='MACSE frameshift cost for less reliable sequences')
    parser.add_argument("--fs_lr_term", type=str, default='10', help='MACSE terminal frameshift cost for less reliable sequences')
    parser.add_argument("--stop_lr", type=str, default='10', help='MACSE premature stop codon cost for less reliable sequences')
    parser.add_argument("--fs", type=str, default='30', help='MACSE frameshift cost for reliable sequences')
    parser.add_argument("--fs_term", type=str, default='10', help='MACSE terminal frameshift cost for reliable sequences')
    parser.add_argument("--stop", type=str, default='50', help='MACSE premature stop codon cost for reliable sequences')
    parser.add_argument("--gap_ext", type=str, default='1', help='MACSE gap extension cost')
    parser.add_argument("--gap_ext_term", type=str, default='0.9', help='MACSE terminal gap extension cost')
    parser.add_argument("--gap_op", type=str, default='7', help='MACSE gap opening cost')
    parser.add_argument("--gap_op_term", type=str, default='6.3', help='MACSE terminal gap opening cost')
    #parser.add_argument("--analysis_name",  help='Name of the analysis (job tittle)')
    parser.add_argument("--species", default=None, help="text file with a list of species on which to perform MACSE alignment")
    #parser.add_argument("--exon_data", help="file with the exon information obtained through first step of pseudochecker")
    parser.add_argument("--main_path", help="path to results directory")
    parser.add_argument("--min_exon_ident", type=int, default=65, help='Defines the minimum exon alignment identity for it to be considered as positive alignment')
    parser.add_argument("--macse_analysis_name", default=None, help='name specifically for MACSE analysis')
    

    args = parser.parse_args()

    ####################
    match = args.match  # Match reward
    mismatch = args.mismatch  # Mismatch reward
    bestfit_status = args.bestfit_status  # Defines if the user selects the best fit alignment scoring scheme (0 = no, 1 = yes)
    findalternativestop = args.find_alternative_stop  # Defines if the user pretends to find a downstream final stop codon in the last exon (0 = no, 1 = yes)
    #file_exs_cds = args.file_exs_cds  # Relative path to file containing the reference species' exons and CDS of the in-study gene
    #file_genomic = args.file_genomic  # Relative path to file containing the target genomic sequence per target species
    #file_additional_cds = args.file_additional_cds  # Relative path to file containing the additional (optional) predetermined coding sequences
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
    #analysisname = args.analysis_name  # Name of the analysis (job tittle)
    min_exon_ident = args.min_exon_ident # Defines the minimum exon alignment identity for it to be considered as positive alignment

    if args.main_path[-1] != '/': #path to results directory
        mainpath = args.main_path + '/'
    else:
        mainpath = args.main_path

    with open(mainpath + 'first_step_data_dill.pkl', 'rb') as f:
        data_needed_for_macse =  dill.load(f)

    cdsposglobal = data_needed_for_macse[0]
    partialsequences = data_needed_for_macse[1]
    excludedsequences = data_needed_for_macse[2]
    no_additional_cds = data_needed_for_macse[3]
    numtotalexons = data_needed_for_macse[4]
    exonoccupancy  = data_needed_for_macse[5]
    additionalcds = data_needed_for_macse[6]
    nrglobal = data_needed_for_macse[7]

    work_dir = os.getcwd() + '/'
    file_path = os.path.abspath(__file__)
    file_dir = file_path.replace('pseudochecker_exclusive_macse.py', '')
     
    
    
    #prune reliable and less reliable fastas to have only sequences wanted for the analysis
    if args.species:
        new_mainpath = mainpath + args.macse_analysis_name + '/'
        if os.path.exists(new_mainpath) == False:
            os.mkdir(new_mainpath)
        prune_seqs(mainpath + 'reliable.fasta', ids=args.species, out=new_mainpath +  'reliable.fasta')
        prune_seqs(mainpath + 'lessreliable.fasta', ids=args.species, out=new_mainpath  + 'lessreliable.fasta')
        num_lines = sum(1 for line in open(new_mainpath +  'lessreliable.fasta'))
        if num_lines == 0:
            nrglobal = 0
        os.chdir(file_dir)
        run_MACSE(new_mainpath, fs, fs_term, fs_lr, fs_lr_term, stop_lr, stop, gap_ext, gap_ext_term, gap_op, gap_op_term, nrglobal)
        
    else:
        print(mainpath)
        new_mainpath = mainpath
        os.chdir(file_dir)
        run_MACSE(new_mainpath, fs, fs_term, fs_lr, fs_lr_term, stop_lr, stop, gap_ext, gap_ext_term, gap_op, gap_op_term, nrglobal)
        
    os.chdir(work_dir)

    print('STEP 1 DONE')
    #codingsequencepredictionparametersframe(new_mainpath, bestfit_status, match, mismatch, utr_status, findalternativestop, min_exon_ident )
    print('STEP 2 DONE')
    order_alignment(new_mainpath)  # Auxiliary Function

    print('STEP 3 DONE')
    #macse_analysis(new_mainpath, cdsposglobal, partialsequences)  # Auxiliary Function

    print('STEP 4 DONE')
    #alignment_table(new_mainpath)  # Auxiliary Function
    print('STEP 5 DONE')
    #additional_alignment_metrics(new_mainpath)
    print('STEP 6 DONE')
    
    
    excluded = excludedsequences[:]
    
    #generalinfo(new_mainpath, excluded, no_additional_cds, partialsequences)  # Auxiliary Function
    
    #computepseudoindex(new_mainpath, partialsequences, numtotalexons, exonoccupancy, additionalcds, excluded)
    print('STEP 7 DONE')

    successful_confirmation_file = open(new_mainpath + 'success.txt', 'w')
    successful_confirmation_file.close()

    
