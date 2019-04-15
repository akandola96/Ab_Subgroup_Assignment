# Antibody Subgroup Assignment

This repository contains the code created for my MSci Biochemistry project at UCL. Created with Andrew Martin. 

Core analysis code consists of data_extraction.py, blast.py, derive_profiles.py, and core_analysis.py. For each distinct organism each of these scripts were repeated. 



## Data Extraction 

### Aims to extract raw data from abysis, culminating in a FASTA formatted file containing sequences to be analysed

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

### Runs BLAST on all query sequences. Culminates in a csv file (referred to as blout queries file) that contains Query ID, Alignment Description, and 21 N-terminal query sequence residues for each query sequence.

The **extract_ref_data** function is used to extract reference query data from a FASTA formatted file. In the case of this project, this data derives from the IMGT. *Organism* represents the name of the query organism, this must be as found in IMGT data i.e. 'Mus musculus' rather than 'Mouse'. *Region* refers to V, D, J or C genes; in this project, V-region genes were exclusively used. The function extracts reference sequences that belong to the organism in question, provided they are functional genes. 

**make_ref_blastdb** utilises the data produced from **extract_ref_date** to construct a BLAST+ database. *-dbtype* is nucleotide and Seq IDs are parsed. Reference database should be named something appropriate e.g. m_musculus_ref_db. Six files will be produced from this step, all of which contribute to the BLAST database. Function works by running a command line process via Python's subprocess module. 

**tBLASTn_full** carries out a commmand line process to run tBLASTn. The *queries* input variable should be the FASTA formatted sequence file derived from the **Data Extraction** steps. The *db* variable derives from the **make_ref_blastdb** function, taking care not to specify any single file from this function but rather the name of the overall db e.g. m_musculus_ref_db rather than m_musculus_ref_db.nin. *-outfmt 5* specifies that the output of this step should be in XML format, this resulting output file is often very large and should not be manually viewed. The alignments in the XML file are ordered by E-value, with the best E-value being the first alignment. the *soft_masking false* argument is passed to turn off low complexity filtering. 

The **count_xml_blast_records** function simply counts the number of BLAST records produced. This number should agree with the number of sequences counted by the **count_num_queries** function, and is used solely to check this. 

**blast_output_xml2csv** converts the output of **tBLASTn_full** from XML to CSV. Extracts only the top alignment (by E-value) for each BLAST record. The query ID and alignment description are outputted to a new CSV file.

**convert_queries2csv** converts FASTA formatted query sequences into a csv file, extracting only the first 21 residues of each sequence. This CSV file is then joined to the BLAST output CSV file to create the blout_queries file. 

**join_queries2blout**: Attaches query sequences (residues) to the output of the file produced by **blast_output_xml2csv**. Resultant file is of the format shown below. Query sequence residues up to the 21st residue are included. This file is referred to as the **blout_queries** for the remainder of this ReadMe

|Query ID|Alignment Description|Query Sequence Residues|
|-----|-----|-----|
|Anti-A\00001|refseq5115\mus musculus\IGHV1-21*04|QVQLVQ...|
|Anti-B\00501|refseq3113\mus musculus\IGKV1-39*02|EVLKST...|

## Derive Profiles 

### Derives profiles for all subgroups. Culminates in a .txt file containing all subgroup profiles. Functions also used in later analysis.

**phrases** is a simple RegEx function that was designed to provide a way of differentiating IGHV1 from IGHV11 (as found in BLAST alignment description). Python's 'in' cannot distinguish between these. Function used throuhgout project from hereon in. 

**determine_subgroups**. Master function for the determination of subgroups. Uses other functions to create a dictionary in the form shown below containing all subgroups.

#### {Mus musculus Heavy Chain 1:IGHV1, Mus musculus Heavy Chain 2:IGHV2, Mus musculus Kappa Chain 1:IGKV1, Mus musculus Lambda Chain 1:IGLV1}

This dictionary is required to manage the difference in outputs of BLAST (which assigns subgroups as 'IGHV1') and hsubgroup (which assigns subgroups as 'Mus musculus Heavy Chain 1'). 

**get_sentences**. Converts the 'numeric' subgroup codes derived from the **get_numerics** function into sentences e.g. IGHV1 to Mus musculus Heavy Chain 1.

**get_numerics**: Root of all functions used to derive subgroup profiles. BLAST alignment description contains subgroup info in the form of IGHV1, IGHV2 etc. **get_numerics** creates strings e.g. IGHV1, IGHV2, IGHV3, IGHV4 ... IGHV30 and determines whether this string can be found in any of the BLAST alignment descriptions. If it is found, this subgroup is added to a list of confirmed subgroups, from which other functions work. *in_file* is the **blout queries file** shown above. 

  **get_profiles**: Master function to generate subgroup profiles. Must be run for each distinct chain e.g. Kappa, Lambda or Heavy ; this is done by changing the *locus* variable. *in_file* is the blout_queries file produced by earlier steps. *freq_type* can be 'log' or 'normal' and dictates whether the frequency of the residues should be count/number of sequences or log((count/number of sequences) + 1.003). 1.003 is present to ensure that ordinary frequencies are converted to non-zero numbers, log(1.003) = 0.001, which is the lowest possible value that can be stored as a frequency in the profiles. *matrix_type* can be 'full' or '2line'. This function calls one of two functions depending on the *matrix_type* variable, detailed below. 
  
  **derive_profiles_2line**:generates a 2line profile for a subgroup (1st and 2nd most common residues only). For each row in the blout_queries file, if the record alignment is found to contain the subgroup in question the query residues are added to distinct list. The first residue is added to the list pos1, the second residue to the list pos2, the third residue to the list pos3 and so on. From each of these lists the most commonly occuring and the second most commonly occuring residues are calculated as well as their frequencies. Passing *freq_type* as 'log' modifies the frequency such that f1 = log(f0 + 1.003), where f0 is the original frequency and f1 is the new frequency. 
  
  **derive_profiles_full**: generates a full profile for a subgroup (frequency of every residue at every position). For each row in the blout_querie. file, if the record alignment is found to contain the subgroup in question the query residues are added to distinct list. The first residue is added to the list pos1, the second residue to the list pos2, the third residue to the list pos3 and so on. From these lists, the frequency of every residue at every position is calculated. 
  
  ## Core Analysis
  
  ### Calculates the MCC for each subgroup. Culminates in a CSV file containing the MCC for each subgroup as well as other performance measures. 
  
  **fasta2pir**: hsubgroup requires that query sequences be in PIR format. This function takes a FASTA formatted queries file as *infile* and converts applies PIR formatting 
  
  **run_hsubgroup**. Takes the .txt file containing all subgroup profiles, as well as the PIR formatted queries file produced by **fasta2pir** as *matrix_file* and *query_file* respectively. Various parameters modify the action of hsubgroup: if the the profiles being used are 2line then *matrix_type* should be set to '2line', otherwise 'full'. If score calculation is desired to be based on logarithmic profiles then *score_type* should be set to 'product', otherwise 'sum'. For eeach query sequence two scores are generated, producing a CSV file of the form shown below:
  
  |Primary subgroup assignment|Score|Secondary subgroup assignment|Score|
  |----|----|----|----|
  |Mus musculus Kappa Chain 1|98.65|Mus musculus Kappa Chain 2|67.02|
  |Mus musculus Lambda Chain 3|59.35|Mus musculus Heavy Chain 14|32.92|
  
  This is referred to as the *scores_file*
  
  **attach_scores_to_queries**. As hsubgroup does not include SeqID information, it must be attached to the scores by means of this function. This function uses the FASTA formatted queries file to extract SeqIDs for each sequene and uses a zip function to write to a csv file of the format shown below:
  
  |SeqID|Primary subgroup assignment|Score|Secondary subgroup assignment|Score2|
  |----|----|----|----|----|
  |MOPC'ACL\00001|Mus musculus Kappa Chain 1|98.65|Mus musculus Kappa Chain 2|67.02|
  |Anti phos CD12\02059|Mus musculus Lambda Chain 3|59.35|Mus musculus Heavy Chain 14|32.92|
  
  
