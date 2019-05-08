#%%
def master_derive_profiles(query_organism,freq_type, matrix_type):
    """
    Summary:
    Master function to derive subgroup profiles
    
    Args:
    query_organism = query organism (str)
    freq_type = 'log' OR 'normal' (str)
    matrix_type = '2line' OR 'full' (str)
    """
    
    loci = ['Heavy','Kappa','Lambda']
    for locus in loci:
        get_profiles(locus, 'blout_queries.csv', 'profiles.txt',query_organism,freq_type,matrix_type)
    


def phrases(phrase,text):
    """
    Summary:
    Alternative to using 'in', which is non-specific.
    
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

def determine_subgroups(in_file,organism):
    """
    Summary:
    Master function to determine subgroups to make profiles for.
    
    Args:
    in_file = input file. Should be blout_queries. (.csv)
    organism = organism name 
    
    Desc:    
    Utilises get_sentences and get_numerics functions to create a dictionary 
    of the form below:
    {Mus Musculus Heavy Chain 1:IGHV1, Mus musculus Heavy Chain 2:IGHV2}.
    Heavy, Lambda, Kappa all present in this dictionary.
    """
    loci = ['Heavy','Kappa','Lambda']
    subgroup_dict = {}
    # Get locus dictionary and append to overall dictionary
    for locus in loci:
        locus_dict = get_sentences(locus,in_file,organism) # Function call
        subgroup_dict.update(locus_dict)     
    return subgroup_dict        # Makes dict available for other functions

def get_sentences(locus, in_file, organism):
    """
    Summary:
    Converts subgroup numeric representation to a sentence and creates 
    dictionary based on this.
    
    Args:
    locus = locus type to examine (str)
    in_file = blout_queries file (.csv)
    organism = query organism (str)
    
    Desc:    
    E.g. IGHV1 --> Mus musculus Heavy Chain 1.  
    
    Dictionary of form:
        {Mus Musculus Heavy Chain 1:IGHV1, Mus musculus Heavy Chain 2:IGHV2}.
    Single chain type only.
    """
    numeric_list = get_numerics(locus,organism,in_file) # Function call
    sentences = []    
    for subgroup_numeric in numeric_list:           # e.g.[IGHV1, IGHV2, IGHV3]
        subgroup_numeric = str(subgroup_numeric)
        subgroup_numeric = subgroup_numeric[4:]     # Gets subgroup number                           
        sentence = organism + ' ' + locus + ' Chain ' + subgroup_numeric 
        sentences.append(sentence)
        
    locus_dict = dict(zip(sentences,numeric_list))  # Combine lists
    
    # Makes dict available for other functions
    return locus_dict
  
def get_numerics(locus,organism, in_file):
    """
    Summary:
    Determines the 'numeric' of each subgroups in each locus.
    
    Args:    
    locus = Kappa, Lambda or Heavy (str)
    organism = organism name e.g. Mus musculus (str)
    in_file = file containing query seqeunces matched to their BLAST record
        (this is the blout_queries file) (.csv)
        
    Desc:    
    BLAST record contains subgroup assignment e.g. IGHV1, IGHV2...
    Function works by creating a string e.g. IGHV1, IGHV2 ... IGHV30 for each
    possible subgroup (up to 30) and determines whether this string is present 
    in the BLAST record. If string is present then this subgroup is added to a 
    list (numeric_list).
    
    O.mykiss and O.cuniculus records are all labelled as intrasubgroup 
    subgroups. There are no simple 'IGHV1' subgroups present. This function 
    handles this by adding an 'S' character to the query e.g. IGHV1 -> IGHV1S1.
    If IGHV1S1 is found in an O.cuniculus or O.mykiss record, it is counted as
    IGHV1. Code is duplicated but will be fixed in later release.
    """
    import csv
    with open(in_file) as csv_in:
        reader = csv.reader(csv_in)
    
        if locus == 'Heavy':
            root = 'IGHV'
        elif locus == 'Kappa':
            root = 'IGKV'
        elif locus == 'Lambda':
            root = 'IGLV'
        numeric_list = []
        
        for x in range(1,31):               # Each chain can have up to 30 subgs
            query = root + str(x)           # e.g. IGHV1, IGHV2...
            csv_in.seek(0)                  # Move to top of BLAST file
            
            if organism != 'Oryctolagus cuniculus' and \
            organism!= 'Oncorhynchus mykiss':            
                for row in reader:
                    blast = str(row[1])
                    if phrases(query,blast) == False:   # If cant find query 
                        continue
                    elif phrases(query,blast) == True:  # If finds query
                        numeric_list.append(query)      # Update list
                        break
            
            # Handles inconsistent labelling
            # If it finds IGHV1S1, this is counted as IGHV1, which is added to
            # the numeric list.
            elif organism == 'Oryctolagus cuniculus':
                sub_query = query+ 'S'
                for row in reader:
                    blast = str(row[1])
                    if sub_query not in blast:
                        continue
                    elif sub_query in blast:
                        numeric_list.append(query)
                        break     
            elif organism == 'Oncorhynchus mykiss':
                sub_query = query+ 'S'
                for row in reader:
                    blast = str(row[1])
                    if sub_query not in blast:
                        continue
                    elif sub_query in blast:
                        numeric_list.append(query)
                        break
                        
        return numeric_list
             
def get_profiles(locus, infile, out_file,query_organism,freq_type,matrix_type):
    # Backup comment 
    """
    Summary:
    Master function to generate the profiles of a locus' subgroups.
    
    Args:    
    locus = 'Heavy', 'Lambda' or 'Kappa'. Func has to be run for each one (str).
    infile = file containing query seqeunces matched to their BLAST record 
            (blout_queries file) (.csv)
    out_file = named output file cotnaining subgroup profiles (.txt)
    query_organism = organism e.g. Mus musculus (str)
    freq_type = 'log' OR 'normal' (str)
    matrix_type = '2line' OR 'full' (str)
    
    Desc:
    Uses the get_numerics function to determine how many subgroups should be
    evaluated. 
    Then uses either derive_profiles_2line or derive_profiles_full to determine
    a subgroup's profile.
    Adds formatting and subgroup titles.
    """    
    import sys 
    sys.stdout = open(out_file,'a')
    
    subgroups = get_numerics(locus,query_organism,infile) 
    organism = query_organism    
    for subgroup in subgroups:             
        chain_type = locus       
        upper_chain_type = chain_type.upper()
        subg_num = subgroup[4:]             # Get subgroup number 
        try:
            # Write titles
            sys.stdout.write('>' + upper_chain_type + ',' + '  ' 
                             + subg_num +',' + '\n')
            sys.stdout.write('"' + organism + ' ' + chain_type + ' ' 
                             + 'Chain' + ' ' + subg_num + '"' + ',')
            
            # Get profiles
            if matrix_type == '2line':
                derive_profiles_2line(infile,subgroup,freq_type,out_file)
            elif matrix_type == 'full':
                derive_profiles_full(infile,subgroup,freq_type,out_file)
            
            # Write terminator
            sys.stdout.write('//' + '\n')
            
        # Error catcher
        except Exception as e:
            print('ERROR OCCURRED')
            print(e)            

          
def derive_profiles_2line(in_file, query_subgroup,freq_type,out_file):      
    """
    Summary:
    Generates a two line profile for a subgroup.
    
    Args:    
    in_file = file containing query seqeunces matched to their BLAST record
                (blout_queries file)
    query_subgroup = subgroup being evaluated (str)
    freq_type = frequency type, either 'log' or 'normal' (str)
    out_fiile = output file (.txt)
    
    Desc:
    tops = list of most common residues 
    freqs = list of most common residue frequencies 
    tops2 = list of 2nd most common residues
    freqs2 = list of 2nd most common residue frequencies
    
    For each subgroup, a list for each position is made. Every sequence's 
    residue at this position is recorded in it's respective list.
    Unknown  residues (X) are removed from each position list.
    Iterates through each position list to get the most common residue. 
    Builds a profile based on this. 
    
    Lambda function removes the most common residue from each position list 
    and the second most common residue for each position is then calculated.

    The intra_subgroups variable handles IGHV1S1 etc. Treats IGHV1S1 the same 
    a IGHV1.   
    
    Frequency of each residue is recorded to 3 d.p.
    """             
    import csv 
    import sys
    import math
    sys.stdout = open(out_file, 'a')
    with open(in_file,'r') as cinput:
        reader = csv.reader(cinput)
        
        # Initialise position lists
        
        N = 21
    
        list_names = [[] for i in range(N)]
        
        tops = []       
        tops2 = []
        freqs = []      
        freqs2 = []
        intra_subgroups = query_subgroup + 'S'
        
        # Fill position lists        
        for row in reader:
            subgroup = str(row[1])         # Contains BLAST record e.g. IGHV1
            
            # If the query subgroup (or intra) is found in the blast record
            if phrases(query_subgroup,subgroup) == True \
            or intra_subgroups in subgroup:
                x = 0      
                
            # Add 1st residue to pos1 list, 2nd residue to pos2 list...
                for i in list_names:
                    residues = list(row[2])         
                    i.append(residues[x])       # Residue 1 into list 1 etc.        
                    x+=1                        # Update for next residue
                               
     # Determine profiles                   
        for i in list_names:
            i = [x for x in i if x != 'X']  # Remove X's from list
            list_length = len(i)            # Determine number of residues
            
                        
            # Calculate most common residue + frequency
            mode = max(set(i), key= i.count) # Get primary mode as str                
            freq_mode = i.count(mode)/list_length # Calc freq of primary mode
            

            # Calculate second most common residue + frequency
            # If the primary mode is not invariant 
            if freq_mode != 1:      
                i = list(filter(lambda a:a != mode,i))  # Remove primary mode                               
                mode2 = max(set(i), key= i.count)       # Get 2nd
                freq_mode2 = i.count(mode2)/list_length # Calc 2nd mode freq
             
            # If primary mode is invariant 
            elif freq_mode == 1: 
                mode2 = mode
                freq_mode2 = freq_mode
            
        
            # Perform logs after
            if freq_type == 'log':
                freq_mode +=1
                freq_mode = math.log(freq_mode,10)
                freq_mode2 +=1
                freq_mode2 = math.log(freq_mode2,10)
            
            
            # Append to consensus profile (primary)
            tops.append(mode)                   
            freqs.append(freq_mode) 
            formatted_freqs = ['%.3f' % x for x in freqs] # Rounding 
            
            # Append to consensus profile (secondary)
            tops2.append(mode2)
            freqs2.append(freq_mode2) 
            formatted_freqs2 = ['%.3f' % x for x in freqs2] # Rounding
        
        # After iterating through all position lists, write and format
        sys.stdout.write('\n' + "".join(tops) + ',' +'\n') 
        sys.stdout.write(",".join(str(x) for x in formatted_freqs) + ',' +'\n')                 
        sys.stdout.write("".join(tops2) + ',' +'\n') 
        sys.stdout.write(",".join(str(x) for x in formatted_freqs2) +'\n')         
        
def derive_profiles_full(in_file, query_subgroup,freq_type,out_file): 
    """Generates a full profile for a subgroup.
    
    in_file = file containing query seqeunces matched to their BLAST record
                (blout_queries file)
    query_subgroup = subgroup being evaluated (str)
    out_fiile = output file (.txt)
    
    For each subgroup a list for each position is made. Every sequence's 
    residue at this position is recorded in the list.
    
    Iterates through rows. Each position list is updated with the residue at
    its position in each sequence.
    
    Calculates the frequency of every amino acid at each position by iterating 
    through list of amino acids.
    """           
    import csv  
    import sys
    import math
    sys.stdout = open(out_file, 'a+')
    with open(in_file,'r') as cinput:
        reader = csv.reader(cinput)
        
    # Creates  lists
        N = 21
        list_names = [[] for i in range(N)]
        
        intra_subgroups = query_subgroup + 'S'
        
    # Populate lists  
        for row in reader:
            subgroup = str(row[1])
            
            # If subgroup is present in blast record
            if phrases(query_subgroup,subgroup) == True \
            or intra_subgroups in subgroup:
                x = 0
                
                # Res1 to list 1, res2 to list 2 ... res21 to list 21
                for position_list in list_names:
                    residues = list(row[2])               # Gets seq
                    position_list.append(residues[x])     
                    x+=1                     
     # Get profiles 
        counter = 0 
        total_store = []
        for position in list_names:
            counter += 1
            position = [x for x in position if x != 'X'] # Remove x's
            length = len(position)                       # How many res?
            aminos = ['A','C','D','E','F','G','H','I','K','L','M','N',
                      'P','Q','R','S','T','V','W','Y','X']
            value_store = []            # Stores freqs for given pos list
            
            # Calculate frequency of each residue within current position list
            for amino in aminos:
                freq = position.count(amino)/length
                if freq_type == 'log':
                    freq +=1
                    freq = math.log(freq,10)
                freq= '%.3f' % freq                                 
                value_store.append(freq)            # Append pos list          
            total_store.append(value_store)         # Overall profile
        
        sys.stdout.write('\n')                      # Line separates 1st entry        
        amino_counter = -1                          
        for amino in aminos:
            amino_counter +=1                       #look at list position 0 
            sys.stdout.write(amino + ' ')           #write residue and space
            
            for l in total_store:                   #for each list 
                value = l[amino_counter]            #get res freq 
                sys.stdout.write(value + ' ')            
            sys.stdout.write('\n')  
#%%
            
def phrases(phrase,text):
    """
    Summary:
    Alternative to using 'in', which is non-specific.
    
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


in_file = 't6_musculus_blout_queries.csv'
query_subgroup = 'IGHV1'
profile_type = 'full'
import numpy as np
import csv  
with open(in_file,'r') as cinput:
    reader = csv.reader(cinput)
    
    # Creates  lists
    N = 21
    list_names = [[] for i in range(N)]
    
    intra_subgroups = query_subgroup + 'S'
    
    # Populate lists  
    for row in reader:
        subgroup = str(row[1])
        
        # If subgroup is present in blast record
        if phrases(query_subgroup,subgroup) is True \
        or intra_subgroups in subgroup:
            x = 0
            
            # Res1 to list 1, res2 to list 2 ... res21 to list 21
            for position_list in list_names:
                residues = list(row[2])               # Gets seq
                position_list.append(residues[x])     
                x+=1                     
 # Get profile 
    counter = 0 
    profile_list_of_lists = []
    for position in list_names:
        counter += 1
        position = [x for x in position if x != 'X'] # Remove x's
        length = len(position)                       # How many res?
        
        aminos_dict = {0:'A',1:'C',2:'D',3:'E',4:'F',5:'G',6:'H',7:'I',8:'K',9:'L',10:'M',11:'N',12:'P',13:'Q',14:'R',15:'S',16:'T',17:'V',18:'W',19:'Y',20:'X'}
        amino_pos_freqs = []
        for amino in aminos_dict:
            amino_value = aminos_dict[amino]
            amino_freq = position.count(amino_value)/length
            amino_pos_freqs.append(amino_freq)
        profile_list_of_lists.append(amino_pos_freqs)
        profile_array = np.column_stack((profile_list_of_lists))
        
 
  # Write profile
    if profile_type == '2line':
        primary_residues = []
        frequencies = []
        for column in profile_array.T:
            column = list(column)
            frequency = max(column)
            index = column.index(frequency)
            residue = aminos_dict[index]
            
            primary_residues.append(residue)
            frequencies.append(frequency)
    
        
        
        print(primary_residues)
        print(frequencies)
        
    elif profile_type == 'full':
        for amino,row in zip(aminos_dict,profile_array):
            amino_value = aminos_dict[amino]
            print(amino_value,row)
            

#%%

aminos = {'A':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,'K':9,'L':10,
                  'M':11,'N':12,'P':13,'Q':14,'R':15,'S':16,'T':17,'V':18,
                  'W':19,'Y':20,'X':21}
    
        