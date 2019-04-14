def phrases(phrase,text):
    """Replaces using 'in', which is non-specific.
    
    Returns true if EXACT match between strings is found, false oteherwise. 
    Case insensitive.
    Can differentiate between IGHV1 and IGHV10, which 'in' cannot.
    """    
    import re
    return re.search(r"\b{}\b".format(phrase), text, re.IGNORECASE) is not None
#%%
def determine_subgroups(in_file,organism):
    """Master function to determine subgroups to make profiles for.
    
    Utilises get_sentences and get_numerics functions to create a dictionary 
    of the form below:
    {Mus Musculus Heavy Chain 1:IGHV1, Mus musculus Heavy Chain 2:IGHV2}
    Heavy, Lambda, Kappa all present in this dictionary.
    """
    loci = ['Heavy','Kappa','Lambda']
    subgroup_dict = {}
    for locus in loci:
        locus_dict = get_sentences(locus,in_file,organism)
        subgroup_dict.update(locus_dict)     
    return subgroup_dict 

#%%     
def get_sentences(locus, in_file, organism):
    """Converts subgroup numeric representation to a sentence and creates 
    dictionary based on this.
    
    E.g. IGHV1 --> Mus musculus Heavy Chain 1.  
    
    Dictionary of form:
        {Mus Musculus Heavy Chain 1:IGHV1, Mus musculus Heavy Chain 2:IGHV2}
    Single chain type only.
    """
    numeric_list = get_numerics(locus,organism,in_file)
    sentences = []    
    for subgroup_numeric in numeric_list:           #e.g. [IGHV1, IGHV2, IGHV3]
        subgroup_numeric = str(subgroup_numeric)
        subgroup_numeric = subgroup_numeric[4:]     #gets subgroup number                           
        sentence = organism + ' ' + locus + ' Chain ' + subgroup_numeric 
        sentences.append(sentence)
        
    locus_dict = dict(zip(sentences,numeric_list))  #combine lists
    
    return locus_dict
#%%
#21/02 added an if statement to deal with o_cuniculus    
def get_numerics(locus,organism, in_file):
    """Determines the number of subgroups in each chain
    
    locus = Kappa, Lambda or Heavy
    organism = organism name e.g. Mus musculus 
    in_file = file containing query seqeunces matched to their BLAST record
    
    BLAST record contains subgroup assignment e.g. IGHV1, IGHV2...
    Function works by creating a string e.g. IGHV1, IGHV2 ... IGHV30 for each
    possible subgroup (up to 30) and determines whether this string is present 
    in the BLAST record. If string is present then this subgroup is added to a 
    list (numeric_list).
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
        
        for x in range(1,31):               #each chain can have up to 30 subgs
            query = root + str(x)           #e.g. IGHV1, IGHV2...
            csv_in.seek(0)                  #move to top of BLAST file
            
            if organism != 'Oryctolagus cuniculus':            
                for row in reader:
                    blast = str(row[1])
                    if phrases(query,blast) == False:   #if cant find query 
                        continue
                    elif phrases(query,blast) == True:  #if finds query
                        numeric_list.append(query)      #update list
                        break
            
            #handles inconsistent labelling in o.cuniculus
            elif organism == 'Oryctolagus cuniculus':   
                for row in reader:
                    blast = str(row[1])
                    if query not in blast:
                        continue
                    elif query in blast:
                        numeric_list.append(query)
                        break
                
        return numeric_list
#%%               
def get_matrices(locus, infile, out_file,query_organism,freq_type,matrix_type):
    """Master function to derive the profiles of a loci's subgroups.
    
    locus = Heavy, Lambda or Kappa
    infile = file containing query seqeunces matched to their BLAST record (bloutqueries)
    out_file = named output file cotnaining subgroup profiles (.txt)
    query_organism = organism 
    freq_type = 'log' or 'normal'
    matrix_type = '2line' OR 'full'
    
    Uses the get_numerics function to determine how many subgroups should be
    evaluated. 
    Adds formatting and subgroup titles.
    
    05/04 Need to amend so accepts 2 line or full as arguments.
    
    """    
    import sys 
    sys.stdout = open(out_file,'a')
    
    subgroups = get_numerics(locus,query_organism,infile) 
    organism = query_organism    
    for subgroup in subgroups:             
        chain_type = locus       
        upper_chain_type = chain_type.upper()
        subg_num = subgroup[4:]             #get subgroup number 
        try:
            sys.stdout.write('>' + upper_chain_type + ',' + '  ' 
                             + subg_num +',' + '\n')
            sys.stdout.write('"' + organism + ' ' + chain_type + ' ' 
                             + 'Chain' + ' ' + subg_num + '"' + ',')
            if matrix_type == '2line':
                derive_matrices_2line_updated(infile,subgroup,freq_type,out_file)
            elif matrix_type == 'full':
                derive_matrices_full(infile,subgroup,freq_type,out_file)
            sys.stdout.write('//' + '\n')
        except Exception as e:
            print('ERROR OCCURRED')
            print(e)            

def derive_matrices_2line_updated(in_file, query_subgroup,freq_type,out_file):      
    """Generates a two line profile for a subgroup
    
    in_file = file containing query seqeunces matched to their BLAST record
    query_subgroup = subgroup being evaluated 
    out_fiile = output file (.txt)
    
    tops = list of most common residues 
    freqs = list of most common residue frequencies 
    tops2 = list of 2nd most common residues
    freqs2 = lisf of 2nd most common residue frequencies
    
    For each subgroup a list for each position is made. Every sequence's 
    residue at this position is recorded in the list.
    Missing  residues (X) are removed from each position list, however the 
    size of the list is maintained.
    Determines the most commmon residue in each position list. Iterate through
    each position list to get the most common residue. Build a profile based
    on this. 
    
    Lambda function removes the most common residue from each position list 
    and the second most common residue for each position is then calculated.    
    
    Frequency of each residue is recorded to 3 d.p.
    
    """             
    import csv 
    import sys
    import math
    sys.stdout = open(out_file, 'a')
    with open(in_file,'r') as cinput:
        reader = csv.reader(cinput)
        
        #initialise position lists
        pos1 = []
        pos2 = []
        pos3 = []
        pos4 = []
        pos5 = []
        pos6 = []
        pos7 = []
        pos8 = []
        pos9 = []
        pos10 = []
        pos11 = []
        pos12 = []
        pos13 = []
        pos14 = []
        pos15 = []
        pos16 = []
        pos17 = []
        pos18 = []
        pos19 = []
        pos20 = []
        pos21 = []
        
        #initialise list of lists
        list_names = [pos1,pos2,pos3,pos4, pos5, pos6, pos7, pos8, pos9, pos10,                     
                      pos11, pos12, pos13, pos14, pos15, pos16, pos17, pos18, 
                      pos19, pos20, pos21]
        
        #initialise lists 
        tops = []       
        tops2 = []
        freqs = []      
        freqs2 = []
        
        #fill position lists        
        for row in reader:
            subgroup = str(row[1])
            
            #if the query subgroup is found in the blast record
            if phrases(query_subgroup,subgroup) == True:
                x = 0                 
            #add 1st residue to pos1 list, 2nd residue to pos2 list...
                for i in list_names:
                    residues = list(row[2])         #residues
                    i.append(residues[x])       #res 1 into list 1 etc.        
                    x+=1   
                               
        #determine profiles        
        #Iterate through position lists            
        for i in list_names:
            list_length = len(i)            #determine number of seqs used
            i = [x for x in i if x != 'X']  #remove X's from list
                        
            #calculate most common residues + frequencies
            mode = max(set(i), key= i.count) #get primary mode as str                
            freq_mode = i.count(mode)/list_length #calc freq of primary mode
            

            #calculate second most common residues + frequencies
            #if the primary mode is not invariant 
            if freq_mode != 1:      #WILL NEED TO FIX!!! WONT LIKE LOGS
                if len(list(set(i))) > 1:               #if other residues present
                    i = list(filter(lambda a:a != mode,i)) #remove primary mode                               
                    mode2 = max(set(i), key= i.count)    #get 2nd
                elif len(list(set(i))) < 2:           #if no other residues present
                    mode2 = mode 
                        
                freq_mode2 = i.count(mode2)/list_length       #calc 2nd mode freq
             
            #if primary mode is invariant 
            elif freq_mode == 1: 
                mode2 = mode
                freq_mode2 = freq_mode
            
            
            if freq_type == 'log':
                freq_mode +=1.003
                freq_mode = math.log(freq_mode,10)
                freq_mode2 +=1.003
                freq_mode2 = math.log(freq_mode2,10)
            
            
            #write primary lists
            tops.append(mode)                   
            freqs.append(freq_mode) 
            formatted_freqs = ['%.3f' % x for x in freqs] #rounding 
            
            #write secondary lists
            tops2.append(mode2)
            freqs2.append(freq_mode2) 
            formatted_freqs2 = ['%.3f' % x for x in freqs2] #rounding
        
        #write lists and format
        sys.stdout.write('\n' + "".join(tops) + ',' +'\n') 
        sys.stdout.write(",".join(str(x) for x in formatted_freqs) + ',' +'\n')                 
        sys.stdout.write("".join(tops2) + ',' +'\n') 
        sys.stdout.write(",".join(str(x) for x in formatted_freqs2) +'\n')       
#%%      
def derive_matrices_full(in_file, query_subgroup,freq_type,out_file): 
    """Generates a full profile for a subgroup
    
    in_file = file containing query seqeunces matched to their BLAST record
    query_subgroup = subgroup being evaluated 
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
        pos1 = []
        pos2 = []
        pos3 = []
        pos4 = []
        pos5 = []
        pos6 = []
        pos7 = []
        pos8 = []
        pos9 = []
        pos10 = []
        pos11 = []
        pos12 = []
        pos13 = []
        pos14 = []
        pos15 = []
        pos16 = []
        pos17 = []
        pos18 = []
        pos19 = []
        pos20 = []
        pos21 = []
        list_names = [pos1,pos2,pos3,pos4, pos5, pos6, pos7, pos8, pos9, pos10, 
                      pos11, pos12, pos13, pos14, pos15, pos16, pos17, pos18, 
                      pos19, pos20, pos21]
        
        
        
        for row in reader:
            subgroup = str(row[1])
            
            #if subgroup is present in blast record
            if phrases(query_subgroup,subgroup) == True: 
            #if query_subgroup in subgroup:  #nonspecific variant. oryctolagus
                x = 0
                #for each sequence, what residue is at position x. Build a 
                #list of all residues at this position.
                for position_list in list_names:
                    residues = list(row[2])               #gets the sequence
                    position_list.append(residues[x])     #build the list
                    x+=1                     
     #each position list now contains every sequence's residue at that position
        counter = 0 
        total_store = []
        for position in list_names:
            counter += 1
            length = len(position)      #length of the list 
            aminos = ['A','C','D','E','F','G','H','I','K','L','M','N',
                      'P','Q','R','S','T','V','W','Y','X']
            value_store = []            #stores frequencies
            
            #calculate frequency of each amino acid at each position
            for amino in aminos:
                freq = position.count(amino)/length #calculates frequncy
                if freq_type == 'log':
                    freq +=1.003
                    freq = math.log(freq,10)
                freq= '%.3f' % freq                                 
                value_store.append(freq)            #appends frequency list          
            total_store.append(value_store)         #appends a list of lists
                
        amino_counter = -1                          
        for amino in aminos:
            amino_counter +=1                       #look at list position 0 
            sys.stdout.write(amino + ' ')           #write the residue and a spcae
            
            for l in total_store:                   #for each list in the list of lists
                value = l[amino_counter]            #get the freq of the residue at each position
                sys.stdout.write(value + ' ')            
            sys.stdout.write('\n')  