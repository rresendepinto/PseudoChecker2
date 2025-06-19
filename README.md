## Overview

PseudoChecker consists of an integrated computational pipeline able to infer the coding status of a given eukaryotic nuclear protein-coding gene in single or multiple species of interest by taking advantage of existent genomic data.
It requires the genomic region where the gene is located from each target species and the exons and coding sequence from a reference species

## Installation

First, you will need to clone the repository into a local folder. Then, you can use on of three methods to install dependencies in order to be able to run the tool: conda, docker and singularity.

### Using conda


Install the current requirements through conda:


	conda update -n base -c defaults conda 
    conda install python=3.6 
    conda install anaconda-client 
    conda install pandas 
    conda install -c conda-forge biopython=1.77 
    conda install -c conda-forge sortedcontainers 
    conda install -c bioconda emboss 
    conda install -c bioconda macse 
    conda install -c bioconda java-jdk 
    conda install -c bioconda seqtk 
    conda install dill
    
### Using singularity

Pull the image like this:

	singularity pull pseudochecker2.sif docker://rrpinto/pseudochecker2.0:latest
	
Then run any pseudochecker command like this:
	
	singularity exec pseudochecker2.sif bash path_to_pseudochecker_folder/process_reference_genome.sh 
	singularity exec pseudochecker2.sif python path_to_pseudochecker_folder/pseudochecker.py --help
	

In addition, you will also need [AGAT](https://github.com/NBISweden/AGAT) to process the annotation files in order to create the input reference with exons and coding sequence for pseudochecker 

## How to run

### Reference creation

First, you will want to create the reference for your desired gene, containing the exons and the coding sequence. 

	>Exon_1    
	CATACTGGAGGGTGGGTAAAAGACTGTAAGGCCTCTGGTATTAGAAATTCCAAGATGAGTGACCTCTGGGTTGTGGCCACCCTCAAGAATGGATGCAGGCTCTTGAGGCATTTCCTCACAATGTCTGCATCAAGAACCCTGACCAAAGTGGAGCACCCTGTCCAACCCTGGACCTAAACCTCAGAATCCTCAGGATACTCAGACTTCCTGAATATCCACAAAACAGACTCCAGGGAAGGACCTCTGAGCCTTAACTTGGACAGGAGCACGAGGAGAAGCTCGAGTTGTTCAGCCTGGGAGGGAATTAGTCCCTATTCTGGGAGCTGATGGTGCCTTCTGATTGGTGCAGCAGGAAGATGCCCAGGGAAGTGTGAGTCCCTGTCCCAAGAGTAGGCCGGCCCCACAGAGGGTAGGTAGCTGGTCACAGAGCACAAGGTTGCAGACGGCATCCCTGAAACAACCAGATAAGGGCTGGATTATTGGGTTATGACTCAGTTTACATTTCACATAGAGAGCAGGACAGGAGGGTTGAAATCTCCAGACCGAGAGTGGGAAATTCTTCAAGGTGCCAAGGGCACTAAACTTAGCAGCTTGCCCCCACCCTTTGCTGGTCCCAGCTTCCTGGAAAGGACAAGAGGAGGAAGAGGAGGCATAGTTTTGCATTGGGGTGGAGGGATTAAGGGGAAGAGGAGTTAAGGACCGGATACAGGATTAAAAGAGAAGGGACAGGAAGGAGATGGAAGCTTGGGGTTTGAGCAAGAGTATGGGGTGCAGAGGGTGTAAGGGTTGGTTGTGTTGGGGTTGGGGGGCGTCTTCCAGGCCAACAATAGAGAGAAACAGGAATTTGAATAGAATAAAGGAGAGATAATTGAGCTCTGGAATCTGATGTCTCAAAAAGTCAAGCACTGGAAGCCCCAGGCCCCTGCTGGTAGGAGGCTGGCATTTCTGGAACAGGACTTGGCAAAAACAGGAAGGTAAGATTTCTTGGGTCTTTTATCTGTTTGACAGCAGCAACCTGGCACCAAAGAACTTTGATCTTCCCTCGTGTGCACTTATCTACCTCATCCTTCCTCTGTCTCATCTTCCTGGAAAAGCAACCCCCACCTCAGTTATCCCCAAGGCACCGCATTACCACCCAAAAGGGACTTCTCCCAGCTTCTCACTCTGAATGTGACTGGCCCAGGAGAGACGACGACTCAGTCTGCTTTTGGGACCCATTCCCATCAATGTCTGCTGATAAGACAAGCAGGCTAAGACAGGAAACAGGGAGGGGAAAACGAAAGCAAGCTGGTGCCAGCCTGCCCTGTGCTCGCACCCTACCAACAGTGCCCCCATCAGAATTTGGGTGTGCTAGTGTGATTGTGTGGTGGGATGACTGTGAAAAAACAACTGTGTGGGAGTGGGACTTTGAGTACAACTGGGGCTATATGATTCTGTCTGGGGTTGTCTGACCATGTGAGGTTGATCCAGTGACAGGTGTGACTGTGAATGTGTAACTACAGGTATGAGTGTGAATGTCCTCACTGGATGCCACTGTGGATCCTCTGCCCCCTTTACTGTACCATCAGCCCCCATGTGCATAACCCATCTCCCCAGCTTCTCTGCCCAGGACTTAGGAGTCCATGGAAAACGGGAGTCCTAAAAAAGGAGGGGTCTGTGGTAATTAGAGCCCAGGAATGTCACAAGATGAGCAATACATCTTCACCTGCGGCTGACTCACCTGGTGCCCCAGGTATAAAAGGCACCCAGATTAAAGGTAGGGAAGGAGGAGGAGAGAGAAGGAAGTATCCAGGCTGAGCATCATGAAGGGGCCCTCACCCACCAGCAGCCTCCTGCTGCTCCTGTTGTTCCTGAGCCCAGGCCCTGGAGGAG    
	>Exon_2    
	CATTCCCACTGGCACTCAGCACTGCATGCTGTACTCAGCTCTACCGACAGCCACTCCCAAACAAGCTACTGAGGAGGGTCATCCGAGTGGAACTTCAGGAAGCTGACGGGGACTGTCACCTCCAAGCCTTCGT    
	>Exon_3    
	GCTTCACCTGTCTCAACGCAGGGTCTGCATCCACCCCCAGAACCGCAGCCTGATTCGGTGGTTTGAACGCCAAGGGAAGATGCTCCAGGGAACTCAGCCCAACCAGAGTTTGGAGCTCAAAGGGAAAATGGGCTGGGGCCCCCAGAAGCCAAAGTAATAAAGCAGTGATGCATAATAGTCTCTGA    
	>CDS    
	ATGAAGGGGCCCTCACCCACCAGCAGCCTCCTGCTGCTCCTGTTGTTCCTGAGCCCAGGCCCTGGAGGAGCATTCCCACTGGCACTCAGCACTGCATGCTGTACTCAGCTCTACCGACAGCCACTCCCAAACAAGCTACTGAGGAGGGTCATCCGAGTGGAACTTCAGGAAGCTGACGGGGACTGTCACCTCCAAGCCTTCGTGCTTCACCTGTCTCAACGCAGGGTCTGCATCCACCCCCAGAACCGCAGCCTGATTCGGTGGTTTGAACGCCAAGGGAAGATGCTCCAGGGAACTCAGCCCAACCAGAGTTTGGAGCTCAAAGGGAAAATGGGCTGGGGCCCCCAGAAGCCAAAGTAA

Be sure that there are no different nucleotides between the exons and the coding sequence. Additionally, the coding sequence must not contain any UTRs (untranslated regions) and all exons must be found within the coding sequence (i.e. there must not be exons that contain only an UTR).

If you wish to create references for many genes from a single genome, you can generate this automatically from the annotation but first you need to process the genome you intend to use as reference with AGAT:

    
	process_reference_genome.sh your-annotation-file.gff your-genome-file.fna genome-prefix workdir
	
This will create a series of files, the ones needed to create references are the \{genome-prefix\}_F_transcripts.fasta and \{genome-prefix\}_F_cds.fasta. This reference can also be used for the online version of [PseudoChecker](http://pseudochecker.ciimar.up.pt/pseudochecker/index.html) .
This script assumes that the installed folder of Agat is in the working directory and should be altered to use a different path.

Then, to extract the desired gene and convert to pseudochecker format, the gene's exons must be extracted from the resulting transcripts file (${genome-prefix}_F_transcripts.fasta) and the coding sequence from the cds file (${genome-prefix}_F_cds.fasta).
For this, we suggest using the following scripts:

	get_exons_mrna.sh genome-prefix gene gene-exons-rna.fa
	python pseudochecker_converter.py -f gene-exons-rna.fa -o out_prefix_reference.fa

This will create two files: one with the exons and mrna extracted for that gene from the agat output files and the other with the prepared reference for pseudochecker.  
In this case, we use the gene symbol to extract the exons and coding sequence. To use another identifier, such as NCBI gene id, you can edit the provided get_exons_mrna.sh bash script.

It is also possible to provide your own reference to pseudochecker, just make sure it is in the same format as the reference provided in example_data (for example, exon header should be >Exon_n). You could also provide your own fasta file with the gene's exons and cds to pseudochecker_converter.py.
The next step will be to run pseudochecker:

	python pseudochecker.py --file_exs_cds  out_prefix_reference.fa --file_genomic genomic_region.fna   --analysis_name prefix_name --main_path output_folder -c number_cores --skip_MACSE 

For large datasets, we recommend using the "--skip_MACSE" flag. With this flag active, the software does not create a multiple sequence alignment for all targets, which can become very slow and less reliable with many targets. 
By default, PseudoChecker2 uses the same algorithm for exon alignments as PseudoChecker (Biopython's implementation of Needleman-Wunsch alignment: pairwise2.align.globalms). However, for large datasets, we recommend using the "--needle" flag, which uses Emboss's implementation of Needleman-Wunsch.

If you intend to create a multiple sequence alignment for a subset of targets, you can do so with:

	pseudochecker_exclusive_macse.py --main_path path_used_in_previous_step_for_analysis_output --species list_of_targets_to_align --macse_analysis_name name_for_output_folder 

If you had previously used "--skip_MACSE" flag and now wish to create a multiple sequence alignment for all targets, you can use this command without the "--targets" flag and without providing a targets file.

## Output

In the output folder, you will have:


-Fasta files with the exons for each target species


-frameshifts.csv - list of the frameshift mutations found and several details about them such as the exon in which they are in, the position in the frame, the number of nucleotides which it consists and whether it is an insertion or a deletion


-stop_codons.csv - list of premature stop codons per target


-species_pseudo_stats.csv - a file with the pseudoindex, a % of the shifted in case of a frameshift  mutation, a % of the truncated frame in case of a premature stop codon, missing % of the frame due to absent exons and a list of the aligning exons


-config_params.json - parameters of the analysis in json format


-predictedcds.fasta - predicted coding sequences for each target species in a nucleotide fasta


-predictedcdstranslated.fasta - predicted coding sequences for each target species in an aminoacid fasta


-exon_alns.json - file in json format with information on exon alignments (necessary for PseudoViz)


-first_step_data_dill.pkl - file in dill format with information on the analysis (necessary for running pseudochecker_exclusive_macse.py and for displaying MACSE alignment in PseudoViz)


For easier visualization of these results you can use the PseudoViz web-tool.

## Pseudoviz

[Pseudoviz](https://bitbucket.org/rresendepinto/pseudoviz/src/master/)  allows you to display the output of pseudochecker2.0 in a more intuitive way, showing: 


* a set of statistics that infer the likelihood of a gene being pseudogenized/lost; 
* the visualization of each exon; 
* the mapping of each exon on the genomic region;
* and alignments of the coding sequences for all or a subset of the target species. 

Just upload these files to PseudoViz:


* exon_alns.json
* reference-file-used-for-pseudochecker.fa
* fasta-with-genomic-regions-used-for-pseudochecker.fa (necessary for showing the exons within the region)

You can also install PseudoViz locally:

### Install and use with docker


* docker pull rrpinto/pseudoviz
* docker run -p 8888:5000 rrpinto/pseudoviz

Go to this link on your browser http://localhost:8888 and the app should be live

### Install and use with singularity


First pull the image:

* singularity pull pseudoviz.sif docker://rrpinto/pseudoviz

Then run it to go inside the image:

* singularity run pseudoviz.sif

Inside the image cd into PseudoViz folder and run:

* python flask_website.py

You can now go to http://193.136.51.225:5000/ on your browser and the app should be running.

## Dendrogram display

To get a phylogenetic overview of the mutations and the inferred PseudoIndex scores, you can create a json file from the results of PseudoChecker2, the reference used and a phylogenetic tree :

* python create_json_tree.py -r out_prefix_reference.fa -t phylogenetic_tree.treefile -f results_folder/ > results_dendrogram.json

Note that it is ok to have species in the tree provided that had not been analysed by PseudoChecker2. These simply won't show up in the dendrogram.

To then view the dendrogram, go into the respective page in PseudoViz and upload the created json. The tree should then be displayed with a color code representing PseudoIndex and the mutations when you hover above a specific node.





## Parameters:
 	
	usage: pseudochecker.py [-h] [--file_exs_cds FILE_EXS_CDS]
                        [--file_genomic FILE_GENOMIC]
                        [--file_additional_cds FILE_ADDITIONAL_CDS]
                        [--match MATCH] [--mismatch MISMATCH]
                        [--bestfit_status BESTFIT_STATUS]
                        [--find_alternative_stop FIND_ALTERNATIVE_STOP]
                        [--min_exon_ident MIN_EXON_IDENT]
                        [--utr_status UTR_STATUS] [--fs_lr FS_LR]
                        [--fs_lr_term FS_LR_TERM] [--stop_lr STOP_LR]
                        [--fs FS] [--fs_term FS_TERM] [--stop STOP]
                        [--gap_ext GAP_EXT] [--gap_ext_term GAP_EXT_TERM]
                        [--gap_op GAP_OP] [--gap_op_term GAP_OP_TERM]
                        [--analysis_name ANALYSIS_NAME]
                        [--main_path MAIN_PATH] [-c CPU] [--skip_MACSE]
                        [--needle] 
						
	optional arguments:
	  -h, --help            show this help message and exit
	  --file_exs_cds FILE_EXS_CDS
							Relative path to file containing the reference species
							exons and CDS of the in-study gene
	  --file_genomic FILE_GENOMIC
							Relative path to file containing the target genomic
							sequence per target species
	  --file_additional_cds FILE_ADDITIONAL_CDS
							Relative path to file containing the additional
							(optional) predetermined coding sequences
	  --match MATCH         Match Reward
	  --mismatch MISMATCH   Mismatch reward
	  --bestfit_status BESTFIT_STATUS
							Defines if the user selects the best fit alignment
							scoring scheme (0 = no, 1 = yes)
	  --find_alternative_stop FIND_ALTERNATIVE_STOP
							Defines if the user pretends to find a downstream
							final stop codon in the last exon (0 = no, 1 = yes)
	  --min_exon_ident MIN_EXON_IDENT
							Defines the minimum exon alignment identity for it to
							be considered as positive alignment
	  --utr_status UTR_STATUS
							Defines if the reference species 1st and/or last exon
							are already UTR flanked or not (0 = not flanked, 1 =
							flanked)
	  --fs_lr FS_LR         MACSE frameshift cost for less reliable sequences
	  --fs_lr_term FS_LR_TERM
							MACSE terminal frameshift cost for less reliable
							sequences
	  --stop_lr STOP_LR     MACSE premature stop codon cost for less reliable
							sequences
	  --fs FS               MACSE frameshift cost for reliable sequences
	  --fs_term FS_TERM     MACSE terminal frameshift cost for reliable sequences
	  --stop STOP           MACSE premature stop codon cost for reliable sequences
	  --gap_ext GAP_EXT     MACSE gap extension cost
	  --gap_ext_term GAP_EXT_TERM
							MACSE terminal gap extension cost
	  --gap_op GAP_OP       MACSE gap opening cost
	  --gap_op_term GAP_OP_TERM
							MACSE terminal gap opening cost
	  --analysis_name ANALYSIS_NAME
							Name of the analysis (job title)
	  --main_path MAIN_PATH
							Main path to the job folder
	  -c CPU, --cpu CPU     number of threads
	  --skip_MACSE
	  --needle              uses biopython implementation of emboss needle (faster 
							but still in testing)
	  --debug               outputs genomicregionshort as fasta and tested
							alignments as json

	pseudochecker_exclusive_MACSE.py
		usage: pseudochecker_exclusive_macse.py [-h] [--match MATCH]
												[--mismatch MISMATCH]
												[--bestfit_status BESTFIT_STATUS]
												[--find_alternative_stop FIND_ALTERNATIVE_STOP]
												[--utr_status UTR_STATUS]
												[--fs_lr FS_LR]
												[--fs_lr_term FS_LR_TERM]
												[--stop_lr STOP_LR] [--fs FS]
												[--fs_term FS_TERM] [--stop STOP]
												[--gap_ext GAP_EXT]
												[--gap_ext_term GAP_EXT_TERM]
												[--gap_op GAP_OP]
												[--gap_op_term GAP_OP_TERM]
												[--species SPECIES]
												[--main_path MAIN_PATH]
												[--min_exon_ident MIN_EXON_IDENT]
												[--macse_analysis_name MACSE_ANALYSIS_NAME]
		optional arguments:
		  -h, --help            show this help message and exit
		  --match MATCH         Match Reward
		  --mismatch MISMATCH   Mismatch reward
		  --bestfit_status BESTFIT_STATUS
								Defines if the user selects the best fit alignment
								scoring scheme (0 = no, 1 = yes)
		  --find_alternative_stop FIND_ALTERNATIVE_STOP
								Defines if the user pretends to find a downstream
								final stop codon in the last exon (0 = no, 1 = yes)
		  --utr_status UTR_STATUS
								Defines if the reference species 1st and/or last exon
								are already UTR flanked or not (0 = not flanked, 1 =
								flanked)
		  --fs_lr FS_LR         MACSE frameshift cost for less reliable sequences
		  --fs_lr_term FS_LR_TERM
								MACSE terminal frameshift cost for less reliable
								sequences
		  --stop_lr STOP_LR     MACSE premature stop codon cost for less reliable
								sequences
		  --fs FS               MACSE frameshift cost for reliable sequences
		  --fs_term FS_TERM     MACSE terminal frameshift cost for reliable sequences
		  --stop STOP           MACSE premature stop codon cost for reliable sequences
		  --gap_ext GAP_EXT     MACSE gap extension cost
		  --gap_ext_term GAP_EXT_TERM
								MACSE terminal gap extension cost
		  --gap_op GAP_OP       MACSE gap opening cost
		  --gap_op_term GAP_OP_TERM
								MACSE terminal gap opening cost
		  --species SPECIES     text file with a list of species on which to perform
								MACSE alignment
		  --main_path MAIN_PATH
								path to results directory
		  --min_exon_ident MIN_EXON_IDENT
								Defines the minimum exon alignment identity for it to
								be considered as positive alignment
		  --macse_analysis_name MACSE_ANALYSIS_NAME
								name specifically for MACSE analysis


