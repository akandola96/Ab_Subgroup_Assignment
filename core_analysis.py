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

def fasta2pir(infile,outfile):
    """
    Summary:
    Converts a FASTA formatted file to PIR format. Required format of 
    hsubgroup.   
    
    Args:
    infile = input file (.fasta)
    outfile = output file (.pir)
    """
    from Bio import SeqIO
    import sys 
    sys.stdout = open(outfile,'a')
    file = SeqIO.parse(infile,'fasta')
    for record in file:
        print('>P1;'+ record.id) 
        print(record.id)            # Used to make PIR format work
        print(record.seq + '*')   

def run_hsubgroup(query_file,matrix_file,matrix_type,out_file):
    """
    Summary:
    Takes PIR formatted query sequences and profiles and runs hsubgroup
    
    Args:
    query_file = PIR formatted queries (.pir)
    matrix_file = file containing subgroup profile (.txt)
    matrix_type = '2line' or 'full' (str)
    score_type = 'normal' or 'log' (str)
    out_file = named out file containing hsubgroup scores (.csv)
    
    Desc:
    Creates a command string and runs using subproces.
    """

    # Having some trouble running. Need use cygwin to do so on windows
    import subprocess
    
    # Profile version
    if matrix_type == 'full':
        matrix_type == '-f'
    elif matrix_type =='2line':
        matrix_type == ''
    
    cmd = ('hsubgroup -d ' + matrix_file + ' -v ' + matrix_type +' ' 
           + query_file + ' > ' + out_file)
       
    hsubgroup_run = subprocess.Popen(cmd,shell=True) 
    hsubgroup_run.wait()
 
        
def attach_scores_to_queries(in_file,score_file,out_file):
    """
    Summary:
    Matches fasta IDs to their score from hsubgroup.
    
    Args:
    infile = queries.fasta
    score_file = csv file containing scores from hsubgroup
    out_file = named filed containing seq IDs matched to scores (.csv)
    
    Desc:
    hsubgroup scores sequences in same order they are entered in. Thus, can 
    use a zip function to match query sequences to their scores. 

    Output file format (known as seqs_scores file):
    Query_ID, Mus musculus Heavy Chain x, 85.2222, Mus musculus Heavy Chain y,
    67.3333.
    """
    from Bio import SeqIO
    import csv 
    import sys
    sys.stdout = open(out_file,'a+')
    with open(score_file,'r') as scores, open(out_file) as csv_out:
        fasta_file = SeqIO.parse(in_file,'fasta')
        for row, record in zip(scores,fasta_file):
            r = row
            x = record.id
            joined = str(x+ ' ' + ',' + ' ' + r)
            sys.stdout.write(joined)
        
def make_final_results(in_file,blout_file,organism,out_file):
    """
    Summary:
    Creates final results file for downstream analysis.
    
    Args:    
    in_file = file containing fasta IDs matched to their hsubgroup results 
    blout_file = blout queries file
    organism = name of query organism (str)
    out_file = named out file containing full results (.csv)
    
    Desc:    
    Uses a zip function to join the in_file and blout_file together.
    Organism is also attached to the end of the row for species comparison
    later on. 

    Output file format:

    QueryID, Assign 1, Score 1, Assign 2, Score 2, QueryID, BLAST record,
    Residues, Organism.
    
    Notes:
    Redundancy in this file of both QueryID being present twice and residues 
    being included. Not a big issue but could be amended in later release.
    """ 
    import csv 
    import sys 
    sys.stdout = open(out_file,'a+')
    
    with open(out_file,'a+') as full_results, open(blout_file,'r') as blout, \
    open(in_file,'r') as seqs_scores:
        
        # Read input rows 
        seqs_scores_reader = csv.reader(seqs_scores)
        blout_reader = csv.reader(blout)
        
        # Write output rows     
        writer = csv.writer(full_results,lineterminator='\n')
        for seqs_scores_row,blout_row in zip(seqs_scores_reader,blout_reader):
            seqs_scores_row = ','.join(seqs_scores_row)
            blout_row = ','.join(blout_row)            
            print(seqs_scores_row + ',' + blout_row + ',' + organism)
                                 

def check_assignment(in_file,blout_queries_file,organism,out_file):
    """
    Summary:
    Determines the MCC of each subgroup and outputs to file.
    
    Args:    
    in_file = full results file (.csv)
    blout_queries_file = blout_queries_file (.csv)
            Used to generate subgroup dictionary
    organism = current query organism (str)
    out_file = named file to contain results of MCC calculations (.txt)
    
    Desc:    
    Function works by comparing hsubgroup's allocation (e.g. Mouse Heavy 1) to
    BLAST's allocation (e.g. IGHV1). Uses a dictionary to link subgroups e.g.
    {Mouse Heavy Chain 1: IGHV1, Mouse Heavy Chain 2: IGHV2}. Dictionary
    derived from derive_profiles.py functions. 
    
    For each row in the full results file, determines whether it is TP, TN, FP
    or FN. Then uses MCC formula to determine the MCC. Repeats for all 
    subgroups.
    """   
    import csv
    import sys
    sys.stdout = open(out_file,'a+')
    with open(in_file) as csv_in:
        reader = csv.reader(csv_in)
        
        print('Subgroup,','MCC,','TPR,','TNR,','PPV,' , 'TP,','FN,','FP,',
        'TN')
            
        av_MCC_num = 0
        av_MCC_count = 0
        MCC_list = []
    
        # Determine subgroup names to evaluate
        my_dict = determine_subgroups(blout_queries_file,organism)

        # Iterate through subgroups using dictionary      
        for x in my_dict: 
            csv_in.seek(0) #remove
            TP = 0 
            FP = 0 
            FN = 0 
            TN = 0 
            misc = 0
            

            alt_x = my_dict[x] + 'S'   # IGHV1S. Deals with intrasubgroup subgs
            
            for row in reader:            
                prediction = str(row[1])  
                actual = str(row[6])
                
                
                # Assign entry to confusion matrix variable
                
                # If Heavy 1 is prediction and IGHV1 is actual TP + 1
                if phrases(x, prediction) == True and \
                phrases(my_dict[x],actual) == True:
                    TP +=1 
                    
                # If Heavy 1 is prediction and IGHV1S1 is actual TP + 1
                elif phrases(x, prediction) == True and \
                alt_x in actual:
                    TP +=1 
                 
                    
                 # If Heavy 1 not prediction and IGHV1 is actual FN + 1   
                elif phrases(x, prediction) == False and \
                phrases(my_dict[x],actual) == True:
                    FN +=1
                    
                # If Heavy 1 not prediction and IGHV1S1 is actual FN + 1
                elif phrases(x, prediction) == False and \
                alt_x in actual:
                    FN +=1
                
        
                # if Heavy 1 is prediction and IGHV1 is not actual FP + 1
                elif phrases(x,prediction) == True and \
                phrases(my_dict[x],actual) == False:
                    FP +=1
                    
                # If Heavy 1 is prediction and IGHV1S1 is not actual FP + 1    
                elif phrases(x,prediction) == True and \
                alt_x not in actual:
                    FP +=1
                    
                    
                # If Heavy 1 is not prediction and IGHV1 is not actual TN + 1
                elif phrases(x,prediction) == False and \
                phrases(my_dict[x],actual) == False:
                    TN +=1     
                # If Heavy 1 is not prediction and IGHV1S1 is not actual TN+1
                elif phrases(x,prediction) == False and \
                alt_x not in actual:
                    TN +=1     
                    
                    
                #Should always be 0
                else:
                    misc +=1
                 
            # Assess performance
            # Try/except used to catch division by 0 (can sometimes occur)
            try:
                MCC = float((TP * TN) - (FP * FN))/ float(((TP+FP)*(TP+FN)*(TN+FP)*
                           (TN+FN)) ** 0.5)
                MCC = round(MCC,5)
                
                #Other performance measures                         
                TPR = round((TP/(TP+FN)),5)
                TNR = round((TN/(TN+FP)),5)
                PPV = round((TP/(TP+FP)),5)
                                
                # Average calc variables
                av_MCC_num += MCC
                av_MCC_count +=1
                MCC_list.append(MCC)
                print(x, ',' ,MCC, ',', TPR, ',',TNR, ',',PPV,',',TP,',',FN,','
                ,FP,',',TN)
            except Exception as e:
                MCC = 'N/A'
                print(x, MCC, 'Error:', e)
                           
        average_MCC = av_MCC_num/av_MCC_count 
        average_MCC = str(average_MCC)
        print('Simple MCC Average:',',',average_MCC) 

        # Makes MCC values available for other functions
        return MCC_list
