#%%
def phrases(phrase,text):
    """
    Summary:
    Replaces using 'in', which is non-specific.
    
    Args:
    phrase = text to find (str)
    text = text to search in (str)
    
    Desc:    
    Returns true if EXACT match between strings is found, false oteherwise. 
    Case insensitive.
    Can differentiate between IGHV1 and IGHV10, which 'in' cannot.
    """    
    import re
    return re.search(r"\b{}\b".format(phrase), text, re.IGNORECASE) is not None
def extract_misclassified_subgroup(in_file, fasta_file, subgroup, sentence,out_file,misclassified_type,pull_assignment):
                    
    """
    Summary:
    Extracts misclassified sequences for a given subgroup.
    
    Args:
    in_file = full results file (.csv)
    fasta_file = fasta queries file (.fasta)
    subgroup e.g. IGKV1 (str)
    sentence e.g. Mus musculus Kappa Chain 1 (str)
    out_file = named output file (.txt)
    misclassified_type = 'FP' or 'FN' (str)
    pull_assignment = True or False (bool)
        If true returns actual assignment to console.
    
    Desc:
    Function works by determining if a particular sequence is misclassified for 
    the subgroup and if so adds its seq.id to a list. Uses FASTA queries 
    file and extracts the relevant sequence. For use in phylo tree creation.    
    """
    import csv
    import sys
    from Bio import SeqIO
    sys.stdout = open(out_file,'a+')
    with open(in_file) as csv_in:
        reader = csv.reader(csv_in)
        
        seq_ids_FP = []                     # False Positives IDs
        seq_ids_FN=[]                       # False negatives IDs
        hsubgroup_assignment_FP = []        
        hsubgroup_assignment_FN = []
        
        for row in reader:
            seq_id = str(row[0])
            assigned_subg = str(row[1])
            actual_subg = str(row[6]) 
            

            # Need to add functionality to deal with intrasubgroup subgs
            # False negative
            # Does not allow IGHV1S1 to be counted as IGHV1.
            if phrases(subgroup, actual_subg) == True and \
            phrases(sentence,assigned_subg) == False:
                seq_ids_FN.append(seq_id)
                hsubgroup_assignment_FN.append(row[1])

            # False positive
            elif phrases(subgroup, actual_subg) == False and \
            phrases(sentence,assigned_subg) == True:
                seq_ids_FP.append(seq_id)
                hsubgroup_assignment_FP.append(row[6])
                
        if misclassified_type == 'FP':
            for seq_id in seq_ids_FP:
                seq_id = seq_id[:-1]    # Remove second space after ID
                for record in SeqIO.parse(fasta_file,'fasta'):
                    if record.id == seq_id:
                        record.id = record.id[-10:] # Keep acc code only
                        record_seq = str(record.seq)
                        header = '>' + subgroup + '_' + misclassified_type \
                                + '|' +record.id
                        sys.stdout.write(header+'\n'+ record_seq+'\n')
            # Returns hsubgroup's prediction to console
            if pull_assignment == True:
                return hsubgroup_assignment_FP
                        
        elif misclassified_type == 'FN':
            for seq_id in seq_ids_FN:
                seq_id = seq_id[:-1] # Remove second space after ID (comes from list)
                for record in SeqIO.parse(fasta_file,'fasta'):
                    if record.id == seq_id:
                        record.id = record.id[-10:] # Keep acc code only
                        record_seq = str(record.seq)
                        header = '>' + subgroup + '_' + misclassified_type + '|' +record.id                       
                        sys.stdout.write(header+'\n'+ record_seq+'\n')
            # Returns BLAST's assignment to console 
            if pull_assignment==True:
                 return hsubgroup_assignment_FN         
    sys.stdout.close
    
def extract_first_21_residues(input_file, output_file):
    """
    Summary:
    For the sequences extracted by extract_misclassified_subgroup, extracts 
    the first 21 residues only.
    
    Args:
    input_file = file containing sequences from previous function (.fasta)
    out_file = named output file (.fasta)
    
    Desc:
    For use in phylo trees.
    """
    from Bio import SeqIO
    import sys 
    
    sys.stdout = open(output_file,'w+')
    with open(input_file,'r') as fasta_in, open(output_file,'w+') as fasta_out:
        for record in SeqIO.parse(fasta_in,'fasta'):
            print('>' + record.id)
            seq = str(record.seq)
            seq = seq[:21]
            print(seq)
                     
def extract_random_TPs(in_file,subgroup,sentence, number,fasta_file,
                        out_file,pull_sequences):

    """
    Summary:
    Extracts random TP sequences for a given subgroup. Used for phylo trees.
    
    Args:
    in_file = full results file (.csv)
    subgroup e.g. IGKV1 (str)
    sentence e.g. Mus musculus Kappa Chain 1 (str)
    number = number of TP sequences desired (int)
    fasta_file = fasta queries file (.fasta)
    out_file = named output file (.txt)
    pull_sequence = 'Y' or 'N' (str)
        If Y, extracts sequence to file.
        
    Desc:
        For use in phylo trees
"""   
    import csv
    import sys
    import random
    from Bio import SeqIO
    sys.stdout = open(out_file,'a+')
    with open(in_file) as csv_in:
        reader = csv.reader(csv_in)
        
        TP_ids_list = []
        TP_ids_list_sample = []
        
        for row in reader:
            seq_id = str(row[0])
            assigned_subg = str(row[1])
            actual_subg = str(row[6]) 
            
            
            # If TP then extract ID to list 
            if phrases(subgroup,actual_subg) == True and \
            phrases(sentence,assigned_subg) == True:
                TP_ids_list.append(seq_id)
        # Pick random SeqIDs from list         
        TP_ids_list_sample = random.sample(TP_ids_list,number)
         
        # Extract sequences from FASTA queries file
        if pull_sequences == 'Y':
            for seq_id in TP_ids_list_sample:
                seq_id = seq_id[:-1]
                for record in SeqIO.parse(fasta_file,'fasta'):
                    if seq_id == record.id:
                        record.id = record.id[-10:]
                        record_seq = str(record.seq)
                        header = '>' + subgroup + '_' + record.id
                        sys.stdout.write(header + '\n' + record_seq + '\n')

        # Dont extract sequences                
        elif pull_sequences == 'N':
            print(TP_ids_list_sample)
        sys.stdout.close 