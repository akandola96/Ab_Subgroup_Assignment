def extract_misclassified_subgroup(in_file, fasta_file, subgroup, sentence,
                                out_file,misclassified_type,pull_assignment):
                    
    """Extracts misclassified sequences for a given subgroup.
    in_file = full results file (.csv)
    fasta_file = fasta queries file (.fasta)
    subgroup e.g. IGKV1 (str)
    sentence e.g. Mus musculus Kappa Chain 1 (str)
    out_file = named output file (.txt)
    misclassified_type = 'FP' or 'FN' (str)
    pull_assignment = True or False (bool)
        If true returns actual assignment to console.
    
    Function works by determinig if a particular sequence is misclassified for 
    the subgroup and if so adds its seq.id to a list. Then goes FASTA queries 
    file and extracts the sequences. For use in phylo tree creation.    
    """
    import csv
    import sys
    from Bio import SeqIO
    sys.stdout = open(out_file,'a+')
    with open(in_file) as csv_in:
        reader = csv.reader(csv_in)
        
        seq_ids_FP = []                     #False Positives IDs
        seq_ids_FN=[]                       #False negatives IDs
        hsubgroup_assignment_FP = []        
        hsubgroup_assignment_FN = []
        
        for row in reader:
            seq_id = str(row[0])
            assigned_subg = str(row[1])
            actual_subg = str(row[6]) 
            

            #need to add functionality to deal with intrasubgroup subgs
            #False negative
            if phrases(subgroup, actual_subg) == True and \
            phrases(sentence,assigned_subg) == False:
                seq_ids_FN.append(seq_id)
                hsubgroup_assignment_FN.append(row[1])

            #False positive
            elif phrases(subgroup, actual_subg) == False and \
            phrases(sentence,assigned_subg) == True:
                seq_ids_FP.append(seq_id)
                hsubgroup_assignment_FP.append(row[6])
                
        if misclassified_type == 'FP':
            for seq_id in seq_ids_FP:
                seq_id = seq_id[:-1]    #remove second space after ID
                for record in SeqIO.parse(fasta_file,'fasta'):
                    if record.id == seq_id:
                        record.id = record.id[-10:] #keep acc code only
                        record_seq = str(record.seq)
                        header = '>' + subgroup + '_' + misclassified_type \
                                + '|' +record.id
                        sys.stdout.write(header+'\n'+ record_seq+'\n')
            #Returns hsubgroup's prediction to console
            if pull_assignment == True:
                return hsubgroup_assignment_FP
                        
        elif misclassified_type == 'FN':
            for seq_id in seq_ids_FN:
                seq_id = seq_id[:-1] #remove second space after ID (comes from list)
                for record in SeqIO.parse(fasta_file,'fasta'):
                    if record.id == seq_id:
                        record.id = record.id[-10:] #keep acc code only
                        record_seq = str(record.seq)
                        header = '>' + subgroup + '_' + misclassified_type + '|' +record.id                       
                        sys.stdout.write(header+'\n'+ record_seq+'\n')
            #Returns BLAST's assignment to console 
            if pull_assignment==True:
                 return hsubgroup_assignment_FN         
    sys.stdout.close
    
def extract_first_21_residues(input_file, output_file):
    """For the sequences extracted by extract_misclassified_subgroup, extracts 
    the first 21 residues only.

    input_file = file containing sequences from previous function (.fasta)
    out_file = named output file (.fasta)

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
               
#%%
def extract_random_sequences(subgroup,number):
    #extracts random sequences by BLAST assignment. No preference to TP, FP etc.
    import sys
    import csv        
    import random
    from Bio import SeqIO
    sys.stdout = open('IGKV1__single_seq_evolution_alignment.txt','a+')
    with open('Mus_musculus_full_results.csv') as csv_in:
            reader = csv.reader(csv_in)
            ids = []
            test_ids = []
            
            count = 0 
            for row in reader:
                blast_assignment = str(row[6])
                seq_id = row[0]
                if subgroup in blast_assignment:
                    ids.append(seq_id)
                    
                    
            for x in range(1,number):
                test_seq_id = random.choice(ids)
                test_ids.append(test_seq_id)
            
                
            for test_id in test_ids:
                test_id = test_id[:-1] #remove extra space
                for record in SeqIO.parse('mus_musculus_queries_e.fasta','fasta'):
                    if test_id == record.id:
                        record.id = record.id[-10:]
                        record_seq = str(record.seq)
                        sys.stdout.write('>'+record.id+'|'+subgroup+'\n'+record_seq+'\n')
#%%            
def extract_random_TPs(in_file,subgroup,sentence, number,fasta_file,out_file,pull_sequences):
#e.g. full_results, IGHV1-, Mus musculus Kappa Chain 1, 10
#picks 10 random true positive sequences for a given subgroup     
    
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
            actual_subg = str(row[6]) #will need to change to row 7 in later runs
            
            
            # If TP then extract 
            if phrases(subgroup,actual_subg) == True and \
            phrases(sentence,assigned_subg) == True:
                TP_ids_list.append(seq_id)
                
        TP_ids_list_sample = random.sample(TP_ids_list,number) #picks number random choices
         
        
        if pull_sequences == 'Y':
            for seq_id in TP_ids_list_sample:
                seq_id = seq_id[:-1]
                
                for record in SeqIO.parse(fasta_file,'fasta'):
                    if seq_id == record.id:
                        record.id = record.id[-10:]
                        record_seq = str(record.seq)
                        
                        header = '>' + subgroup + '_' + record.id
                        sys.stdout.write(header + '\n' + record_seq + '\n')
                        
        elif pull_sequences == 'N':
            print(TP_ids_list_sample)
            
        sys.stdout.close
        
    