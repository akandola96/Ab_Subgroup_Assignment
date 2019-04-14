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
|*|*|*|
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

The **remove_short_sequences** function takes the output of **remove_spaces** and removes FASTA sequences of less than 21 residues. Outputs to a new file. 

The **seqkit_clean** function takes a FASTA formatted file and removes duplicate sequences by both ID and sequence. The *os* variable is provided to distinguish between Linux and Windows operating systems. The function runs a commmand line process utilisng the SeqKit library. 

The **convert_seqkit** function takes the output of **seqkit_clean** as an input. After running SeqKit the FASTA file will be formatted in a manner that is difficult to manually look at. This function simply re-formats the file so it can be viewed manually and understood. 

The **count_num_queries** function simply counts the number of FASTA sequences present in an *in_file* and returns output to the console

## BLAST 

The **extract_ref_data** function is used to extract reference query data from a FASTA formatted file. In the case of this project, this data derives from the IMGT. *Organism* represents the name of the query organism, this must be as found in IMGT data i.e. 'Mus musculus' rather than 'Mouse'. *Region* refers to V, D, J or C genes; in this project, V-region genes were exclusively used. The function extracts reference sequences that belong to the organism in question, provided they are functional genes. 

**make_ref_blastdb** utilises the data produced from **extract_ref_date** to construct a BLAST+ database. *-dbtype* is nucleotide and Seq IDs are parsed. Reference database should be named something appropriate e.g. m_musculus_ref_db. Six files will be produced from this step, all of which contribute to the BLAST database. Function works by running a command line process via Python's subprocess module. 

**tBLASTn_full** carries out a commmand line process to run tBLASTn. The *queries* input variable should be the FASTA formatted sequence file derived from the **Data Extraction** steps. The *db* variable derives from the **make_ref_blastdb** function, taking care not to specify any single file from this function but rather the name of the overall db e.g. m_musculus_ref_db rather than m_musculus_ref_db.nin. *-outfmt 5* specifies that the output of this step should be in XML format, this resulting output file is often very large and should not be manually viewed. The alignments in the XML file are ordered by E-value, with the best E-value being the first alignment. the *soft_masking false* argument is passed to turn off low complexity filtering. 

The **count_xml_blast_records** function simply counts the number of BLAST records produced. This number should agree with the number of sequences counted by the **count_num_queries** function, and is used solely to check this. 

**blast_output_xml2csv**


