# Methods related with the MACSE analysis
from Bio import SeqIO


def order_alignment(mainpath):
    '''Forces the Reference Species CDS to appear in the top of the alignment'''
    results_nt = open(mainpath + 'macseanalysis_nt.fasta', 'r')
    results_nt_ordered = open(mainpath + 'macseanalysis_nt_ordered.fasta', 'w')
    results_nt = list(SeqIO.parse(results_nt, 'fasta'))  # Parses the MACSE Alignment in Nucleotides
    results_aa = open(mainpath + 'macseanalysis_aa.fasta', 'r')
    results_aa_ordered = open(mainpath + 'macseanalysis_aa_ordered.fasta', 'w')
    results_aa = list(SeqIO.parse(results_aa, 'fasta'))  # Parses the MACSE Alignment in AminoAcids
    for i in results_nt:
        if i.id == 'Reference_Species':
            results_nt_ordered.write('>Reference_Species\n' + str(i.seq) + '\n')
    for i in results_nt:
        if i.id != 'Reference_Species':
            results_nt_ordered.write('>' + str(i.id) + '\n' + str(i.seq) + '\n')
    for i in results_aa:
        if i.id == 'Reference_Species':
            results_aa_ordered.write('>Reference_Species\n' + str(i.seq) + '\n')
    for i in results_aa:
        if i.id != 'Reference_Species':
            results_aa_ordered.write('>' + str(i.id) + '\n' + str(i.seq) + '\n')
    results_nt_ordered.close()
    results_aa_ordered.close()


def macse_analysis(mainpath):
    
    results_nt = open(mainpath + 'macseanalysis_nt_ordered.fasta', 'r')
    results_nt = list(SeqIO.parse(results_nt, 'fasta'))  # Parses the MACSE Alignment in Nucleotides
    #results_aa = open(mainpath + 'macseanalysis_aa_ordered.fasta', 'r')
    #results_aa = list(SeqIO.parse(results_aa, 'fasta'))  # Parses the MACSE Alignment in AminoAcids
    
    for i in results_nt:
        if i.id == "Reference_Species":
            if '!' in str(i.seq):
                invalid_macse_costs = open(mainpath + 'invalid_macse_costs.txt', 'w')
                invalid_macse_costs.close()
    

def get_cds_to_exon_fs(sequence, index, id, cdsposglobal):
    '''Returns the MACSE frameshift detected mutation corresponding exon'''
    t = -1
    for i in range(index + 1):
        if sequence[i] not in {'-', '!'}:
            t += 1
    a = cdsposglobal[id]
    if t < 0:
        t = 0
    for i in a:
        if t >= a[i][0] and t <= a[i][1]:
            return i


def get_cds_to_exon_stop(sequence, index, id, cdsposglobal):
    '''Returns the MACSE premature stop codon detected mutation corresponding exon'''
    t = 0
    for i in range(index):
        if sequence[i] not in {'-', '!'}:
            t += 1
    a = cdsposglobal[id]
    for i in a:
        if t * 3 >= a[i][0] and t * 3 <= a[i][1]:
            return i
