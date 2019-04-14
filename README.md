# Antibody Subgroup Assignment

This repository contains the code created for my MSci Biochemistry project at UCL. Created with Andrew Martin. 

Core analysis code consists of data_extraction.py, blast.py, derive_profiles.py, and core_analysis.py. For each distinct organism each of these scripts were repeated. 



## Data Extraction 

The **xml_parser** function is responsible for extracting sequence data from raw abYsis files. It checks each record within the 
*input_file* provided and extracts all records that belong to the *query_organism*. The *output_file* is a CSV file with entries stored in a vertical format, shown below. The first column contains the residue itself; the second column, the number of this residue within the entire sequence; the third column the Chothia (or Kabat when Chothia is unavailable) numbering of this residue. 

|Col1|Col2|Col3|
|----|----|----|
|>| ID| Accession |
|H|1|H1|
|L|2|H2|
|>| ID| Accession| 
H|3|H1|
V|4|H2|
L|5|H3|
|*|*|*|

The **count_abysis** function takes the output of the **xml_parser** function as it's *in_file* and determines how many records were 
extracted by counting the number of '>' characters it finds. 

The **original_make_fasta** takes the output of the **xml_parser** function as its input and converts ALL sequences into FASTA formatted 
sequences. It outputs a FASTA formatted file. 

The **make_fasta** function is a more sophisticated version of the **original_make_fasta** function. It takes the output of the **xml_parser** function as input and selects which sequences should be converted into FASTA formatted sequences based on user input. If 
the *version* variable is set to 'EA' then all sequences that begin Chothia numbering on values other than 1 are excluded; all other sequences are excluded. If the *version* variable is set to 'T6' then sequences that start Chothia numbering on a maximum of 7 are included however placeholder residues are inserted. A sequence that begins Chothia numbering on H5, would thus have four placeholder residues inserted before its sequence begins. 

The **remove_spaces** function takes the output of the **make_fasta** or **original_make_fasta** functions as an input. It removes all 
spaces from the input file. This is to prevent inconsistencies in either the sequence name or accession code causing issues downstream. It was noticed, for example, that some accession codes had spaces within them e.g. 0 001, which disrupted analysis. 

The **remove_short_sequences** function takes the output of **remove_spaces** and removes FASTA sequences of less than 21 residues.

The **seqkit_clean** function takes a FASTA formatted file and removes duplicate sequences by both ID and sequence. The *os* variable is provided to distinguish between Linux and Windows operating systems. The function runs a commmand line process utilisng the SeqKit...


