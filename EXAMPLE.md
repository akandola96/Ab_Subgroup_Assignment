# Example code run

#### This file contains a simple code run based on Mus musculus and the data files present in the data folder. The placeholder approach, detailed in the ReadMe is employed.

## Core Analysis

### data_extraction.py 
1) Run the *xml_parser* function on an input XML file: `xml_parser('kabat.xml','mus musculus','abysis.csv')` 

This will produce a CSV file containing extracted sequence data. 
If multiple XML files are present, each should be parsed individually and the resultant CSV files concatenated on the command line. 

2) Run the *make_fasta* function on the CSV file produced from step 1 to extract FASTA formatted query sequences: `make_fasta('abysis.csv','T6','raw_queries.fasta)`  
-Extracts queries into FASTA format    
-'T6' employs placeholder approach

3) Run the *remove_spaces* function on the output file of the previous step: `remove_spaces('raw_queries.fasta')`    
-Removes spaces within the file 

4) Run *remove_short_sequences*: `remove_short_sequences('raw_queries.fasta','raw_queries_b.fasta')`  
-Removes short sequences of length < 21  

5) Run *remove_seqs_missing_residues*: `remove_seqs_missing_residues('raw_queries_b.fasta','T6','raw_queries_c.fasta')`  
-Removes sequences missing residues in their N-terminus   
-Nature of sequences removed depends on whether version is set to 'EA' or 'T6'. Needs to be consistent.  

6) Run *seqkit_clean*: `seqkit_clean('raw_queries_c.fasta','windows','raw_queries_d.fasta')`  
Cleans queries for duplicates using SeqKit.  

7) Run *convert_seqkit*: `convert_seqkit('raw_queries_d.fasta',queries.fasta')`      
-Converts SeqKit output into more readable format.  

  
-N.B. I tried to integrate the above steps into one function however sequence loss was observed. Will be amended in a later release.
  
8) Optional - Run the *tidy_up* function to remove redundant files: `tidy_up()`  

### blast.py
9) Run the *blast_steps* function using the query sequence file created in step 2 and the IMGT reference data found in the *data* folder to perform a tBLASTn search: `blast_steps('queries.fasta','imgtrefseqs.fasta','mus musculus','ref_db')`  
-Extracts reference sequences from IMGT data  
-Makes a BLAST database based on this  
-Runs tBLASTN  
This step will output a reference database in addition to an XML file containing BLAST results.    
This step can take some time depending on the number of queries being evaluated.    

10) Run the *blast_output_formatting* function to organise the results of step 4 and create the blout_queries file: `blast_output_formatting('queries.fasta')`  
-Parses the XML formatted BLAST output and converts it to CSV.  
-Extracts the first 21 residues of query sequences and stores in a CSV.  
-Combines these two CSV files.  
For each query sequence the top hit will be extracted.  
The result of these steps is the *blout_queries* file, the format of which can be seen in the ReadMe file.   

### derive_profiles.py

11) Run the master_derive_profiles function to generate subgroup profiles: `master_derive_profiles('Mus musculus','normal','2line')`  
'normal' refers to normal residue frequencies i.e. not 'log'  
'2line' refers to the profile style; the top 2 most commonly occuring residues at each position rather than a 'full' matrix  

### core_analysis.py

12) Run the *fasta2pir* function to convert queries into hsubgroup compatible form: `fasta2pir('queries.fasta','PIR_queries.pir')`  

13) Run *run_hsubgroup* to score input sequences: `run_hsubgroup('PIR_queries.pir','profiles.txt','2line','hsub_scores.csv')`    

14) Run the *master_core_analysis* function to produce an output file containing subgroup MCC values: `master_core_analysis('queries.fasta','hsub_scores.csv','Mus musculus','MCC.csv')`  
-Joins hsubgroup scores to query sequence IDs --> known as the seq_scores file  
-Joins seqs_scores file to blout_queries file  --> known as the final_results file  
-Checks assignment  



#### The results of running these commands can be seen in the 'example' folder. The BLAST output (.xml) was deleted due to the size of the file. 
#### Note that IGLV5 does not have an MCC value associated due to division by 0.
