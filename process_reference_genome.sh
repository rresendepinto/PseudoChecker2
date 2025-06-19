GFF=$1
FNA=$2
NAME=$3
DIR=$4

# extract the largest isoform of each protein
 AGAT/bin/agat_convert_sp_gxf2gxf.pl -g "$DIR""$GFF" -o "$DIR"B_"$NAME".gff3 #General Corrections of GFF3
  AGAT/bin/agat_sp_statistics.pl --gff "$DIR"B_"$NAME".gff3 --gs "$DIR""$FNA" -d -o "$DIR"B_"$NAME"_stats.txt #Stats of GFF3 corrected
   AGAT/bin/agat_sp_fix_overlaping_genes.pl -f "$DIR"B_"$NAME".gff3 -o "$DIR"C_"$NAME".gff3 #Correct Overllaping genes
   AGAT/bin/agat_sp_add_start_and_stop.pl --gff "$DIR"C_"$NAME".gff3 --fasta "$DIR""$FNA" --out "$DIR"D_"$NAME".gff3 #Add missing Start/Stop codons
     AGAT/bin/agat_sp_filter_incomplete_gene_coding_models.pl -v -gff "$DIR"D_"$NAME".gff3 --fasta "$DIR""$FNA" -o "$DIR"E_"$NAME".gff3 #Remove Incomplet gene models
      AGAT/bin/agat_sp_statistics.pl --gff "$DIR"E_"$NAME".gff3 --gs "$DIR""$FNA" -d -o "$DIR"E_"$NAME".stats.txt #Stats of GFF3 corrected
       AGAT/bin/agat_sp_keep_longest_isoform.pl --gff "$DIR"E_"$NAME".gff3 -o "$DIR"F_"$NAME".gff3
        AGAT/bin/agat_sp_statistics.pl --gff "$DIR"F_"$NAME".gff3 --gs "$DIR""$FNA" -d -o "$DIR"F_"$NAME".stats.txt #Stats of GFF3 corrected
         AGAT/bin/agat_sp_extract_sequences.pl -g "$DIR"F_"$NAME".gff3 -f "$DIR""$FNA" -p protein -o "$DIR""$NAME"_F_proteins.fasta
         AGAT/bin/agat_sp_extract_sequences.pl -g "$DIR"F_"$NAME".gff3 -f "$DIR""$FNA" -t cds -o "$DIR""$NAME"_F_cds.fasta
          AGAT/bin/agat_sp_extract_sequences.pl -g "$DIR"F_"$NAME".gff3 -f "$DIR""$FNA" -t mrna -o "$DIR""$NAME"_F_mrna.fasta
           AGAT/bin/agat_sp_extract_sequences.pl -g "$DIR"F_"$NAME".gff3 -f "$DIR""$FNA" -t gene -o "$DIR""$NAME"_F_genes.fasta
            AGAT/bin/agat_sp_extract_sequences.pl -g "$DIR"F_"$NAME".gff3 -f "$DIR""$FNA" -t exon -o "$DIR""$NAME"_F_transcripts.fasta