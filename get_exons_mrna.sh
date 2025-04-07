#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

exons_fasta="${1}"_F_transcripts.fasta
cds_fasta="${1}"_F_cds.fasta

GN=$2

out=$3

grep -i "=gene-${GN} " ${exons_fasta} > Headers_exons.tmp
sed -i 's|[<>,]||g' Headers_exons.tmp

seqtk subseq ${exons_fasta} Headers_exons.tmp > ${out}

rm Headers_exons.tmp

grep -i "=gene-${GN} " ${cds_fasta} > Header_cds.tmp
sed -i 's|[<>,]||g' Header_cds.tmp

seqtk subseq ${cds_fasta} Header_cds.tmp >> ${out}

rm Header_cds.tmp

#convert to pseudochecker reference input
#python ${parent_path}/pseudochecker_converter.py -f ${out}_exons_mrna.fa -o ${out}_reference.fa