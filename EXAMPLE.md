# Example code run

#### This file contains a simple code run through based on the data files present in the data folder. 

## Core Analysis

### data_extraction.py 
1) Run the *xml_parser* function on an input XML file: `xml_parser('kabat.xml','mus musculus','abysis.csv')` 

This will produce a CSV file containing extracted sequence data. 
If multiple XML files are present, each should be parsed individually and the resultant CSV files concatenated on the command line. 

2) Run the *data_extractor* function on the CSV file produced from step 1 to extract FASTA formatted query sequences: `data_extractor('abysis.csv','mus musculus','T6','windows','queries.fasta')`
The 'T6' refers to use of the placeholder approach. 
This will produce a FASTA formatted queries file 

3) Optional - Run the *tidy_up* function to remove redundant files: `tidy_up()`

### blast.py
4) Run the *blast_steps* function using the query sequence file created in step 2 and the IMGT reference data found in the *data* folder to perform a tBLASTn search: `blast_steps('queries.fasta','imgtrefseqs.fasta','mus musculs','ref_db')`
Extracts reference sequences from IMGT data  
Makes a BLAST database based on this 
Runs tBLASTN
This step will output a reference database in addition to an XML file containing BLAST results.   
This step can take some time depending on the number of queries being evaluated.  

5) Run the *blast_output_formatting* function to organise the results of step 4 and create the blout_queries file: `blast_output_formatting('queries.fasta')`
Parses the XML formatted BLAST output and converts it to CSV.  
Extracts the first 21 residues of query sequences and stores in a CSV.  
Combines these two CSV files.  
For each query sequence the top hit will be extracted.  
The result of these steps is the *blout_queries* file, the format of which can be seen in the ReadMe file.   

### derive_profiles.py

6) Run the master_derive_profiles function to generate subgroup profiles: `master_derive_profiles('Mus musculus','normal','2line')`
'normal' refers to normal residue frequencies i.e. not 'log'
'2line' refers to the profile style; the top 2 most commonly occuring residues at each position rather than a 'full' matrix

### core_analysis.py

7) Run the *master_core_analysis* function to produce an output file containing subgroup MCC values: `master_core_analysis('queries.fasta','2line','Mus musculus')`
Converts queries to PIR  
Runs hsubgroup   
Joins hsubgroup scores to query sequence IDs --> known as the seq_scores file  
Joins seqs_scores file to blout_queries file  --> known as the final_results file  
Checks assignment  
