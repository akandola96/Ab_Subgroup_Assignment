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
#%%
def derive_matrices_2line(in_file, query_subgroup,freq_type,out_file):      
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
    Determine the most commmon residue in each position list. Iterate through
    each position list to get the most common residue. Build a profile based
    on this. 
    
    Lambda function removes the most common residue from each position list 
    and the second most common residue for each position is then calculated.    
    
    Frequency of each residue is recorded.
    
    Fixed bug that appeared 04/04 --> incorrect movement along residues
    """ 

    #would be good to tidy up after exams etc.#tag repo with release so differentiate between marking and not marking code            
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
                    i.append(residues[x])           
                    x+=1   
                               
        #determine profiles
        
        #Iterate through position lists            
        for i in list_names:            
            #calculate most common residues + frequencies
            mode = max(set(i), key= i.count) #get primary mode as str
            freq_mode = i.count(mode)/len(i) #calc freq of primary mode
            if freq_type == 'log':
                freq_mode +=1
                freq_mode = math.log(freq,10)
                
            #write lists
            tops.append(mode)                   
            freqs.append(freq_mode) 
            formatted_freqs = ['%.3f' % x for x in freqs] #round
            
            #calculate second most common residues + frequencies
            #if the primary mode is not invariant 
            if freq_mode != 1: 
                i_2 = list(filter(lambda a:a != mode,i)) #remove primary mode
                mode2 = max(set(i_2), key= i_2.count)    #get 2nd mode
                freq_mode2 = i.count(mode2)/len(i)       #calc 2nd mode freq
            
            #if primary mode is invariant 
            elif freq_mode == 1: 
                mode2 = mode
                freq_mode2 = freq_mode
            #write lists
            tops2.append(mode2)
            freqs2.append(freq_mode2) 
            formatted_freqs2 = ['%.3f' % x for x in freqs2] 
        
        #write lists and format
        sys.stdout.write('\n' + "".join(tops) + ',' +'\n') 
        sys.stdout.write(",".join(str(x) for x in formatted_freqs) + ',' +'\n')                 
        sys.stdout.write("".join(tops2) + ',' +'\n') 
        sys.stdout.write(",".join(str(x) for x in formatted_freqs2) +'\n')

#%%
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
    Determine the most commmon residue in each position list. Iterate through
    each position list to get the most common residue. Build a profile based
    on this. 
    
    Lambda function removes the most common residue from each position list 
    and the second most common residue for each position is then calculated.    
    
    Frequency of each residue is recorded.
    
    Fixed bug that appeared 04/04 --> incorrect movement along residues
    """ 

    #would be good to tidy up after exams etc.#tag repo with release so differentiate between marking and not marking code            
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
                freq_mode +=1
                freq_mode = math.log(freq_mode,10)
                freq_mode2 +=1
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
                    freq +=1
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

#%%
def fasta2pir(infile,outfile):
    """Converts a FASTA formatted file to PIR format
    
    PIR formatted sequences are inputted to the hsubgroup program   
    """
    from Bio import SeqIO
    import sys 
    sys.stdout = open(outfile,'a')
    file = SeqIO.parse(infile,'fasta')
    for record in file:
        print('>' + 'P1' + ';' + record.id) #change to one string
        print(record.id)
        print(record.seq + '*')        #PIR strictly ID is not supposed to be more than 6 or 8 chars
        #next line PIr is meant to be a descriptor 
#%%
def run_hsubgroup(query_file,matrix_file,matrix_type,score_type,out_file):
    """Takes PIR formatted query sequences and profiles and runs hsubgroup"""
    import subprocess
    
    if matrix_type == 'full':
        matrix_type == '-f'
    elif matrix_type =='2line':
        matrix_type == ''
        
        
    if score_type == 'product':
        score_type == ' -p'
    elif score_type == 'sum':
        score_type = ''
    
    
    cmd = ('hsubgroup -d ' + matrix_file + ' -v ' + matrix_type + score_type +' ' + query_file + ' >> ' + out_file)
    
    
    hsubgroup_run = subprocess.Popen(cmd,shell=True)
    
    hsubgroup_run.wait()
    
    
#%%
#attach queries to their corresponding score         
def fasta2csv(in_file,score_file,out_file):
    """Matches fasta IDs to their score from hsubgroup
    
    infile = queries.fasta
    score_file = csv file containing scores from hsubgroup
    out_file = named filed containing seq IDs matched to scores (.csv)
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
   
#%%        
def attach_blast_record(in_file,blout_file,organism,out_file):
    """Creates final results file for downstream analysis.
    
    in_file = file containing fasta IDs matched to their hsubgroup results 
    blout_file = blout queries file
    organism = name of query organism 
    out_file = named out file containing full results (.csv)
    
    Uses a zip function to join the in_file and blout_file together.
    Organism is also attached to the end of the row for species comparison
    later on. 
    """ 
    import csv 
    import sys 
    sys.stdout = open(out_file,'a+')
    
    with open(out_file,'a+') as full_results, open(blout_file,'r') as blout, \
    open(in_file,'r') as seqs_scores:
        
        #read input rows 
        seqs_scores_reader = csv.reader(seqs_scores)
        blout_reader = csv.reader(blout)
        
        #write output rows     
        writer = csv.writer(full_results,lineterminator='\n')
        for seqs_scores_row,blout_row in zip(seqs_scores_reader,blout_reader):
            seqs_scores_row = ','.join(seqs_scores_row)
            blout_row = ','.join(blout_row)            
            print(seqs_scores_row + ',' + blout_row + ',' + organism)
                                 
#%%
def check_assignment(in_file,blout_queries_file,organism,out_file):
    """Determines the MCC of each subgroup and outputs to file
    
    in_file = full results file 
    blout_queries_file = blout_queries_file
    organism = current query organism
    out_file = named file to contain results of MCC calculations (.txt)
    
    Function works by comparing hsubgroup's allocation (e.g. Mouse Heavy 1) to
    BLAST's allocation (e.g. IGHV1). Uses a dictionary to link subgroups e.g.
    {Mouse Heavy Chain 1: IGHV1, Mouse Heavy Chain 2: IGHV2}
    
    For each row in the full results file, determines whether it is TP, TN, FP
    or FN. Then uses MCC formula to determine the MCC. Repeats for all 
    subgroups.
    """   
    import csv
    import sys
    sys.stdout = open(out_file,'a+')
    with open(in_file) as csv_in:
        reader = csv.reader(csv_in)
        
        print('Subgroup,','MCC,','TPR,','TNR,','PPV')
            
        
        #initialise variables
        av_MCC_num = 0
        av_MCC_count = 0
        MCC_list = []
    
        #determine subgroup names to evaluate
        my_dict = determine_subgroups(blout_queries_file,organism)
        #will have to repeat this for multiple organisms
        #you probably dont need a blout_queries_file? full_results has same
        #info that blout_queries_does
        
        #iterate through subgroups using dictionary      
        for x in my_dict: 
            csv_in.seek(0)
            #initialise confusion matrix variables
            TP = 0 #
            FP = 0 
            FN = 0 
            TN = 0 
            misc = 0
            training = 0 
            
            #iterate through input file
            for row in reader:            
                prediction = str(row[1])  
                actual = str(row[6])
                #for good code have vbl name which says what
                
                if organism != 'Oryctolagus cuniculus':
                                
                    #assign entry to confusion matrix variable
                    if phrases(x, prediction) == True and \
                    phrases(my_dict[x],actual) == True:
                        TP +=1 
                        
                    elif phrases(x, prediction) == False and \
                    phrases(my_dict[x],actual) == True:
                        FN +=1
                        
                    elif phrases(x,prediction) == True and \
                    phrases(my_dict[x],actual) == False:
                        FP +=1
                        
                    elif phrases(x,prediction) == False and \
                    phrases(my_dict[x],actual) == False:
                        TN +=1
                    
                    #should always be 0
                    else:
                        misc +=1
                        
                elif organism == 'Oryctolagus cuniculus':
                
                    if phrases(x, prediction) == True and \
                    my_dict[x] in actual:
                        TP +=1 
                        
                    elif phrases(x, prediction) == False and \
                    my_dict[x] in actual:
                        FN +=1
                        
                    elif phrases(x,prediction) == True and \
                    my_dict[x] not in actual:
                        FP +=1
                        
                    elif phrases(x,prediction) == False and \
                    my_dict[x] not in actual:
                        TN +=1
                    
                    else:
                        misc +=1
            
            #assess performance
            #get rid of int label
            try:
                MCC = float((TP * TN) - (FP * FN))/ float(((TP+FP)*(TP+FN)*(TN+FP)*
                           (TN+FN)) ** 0.5)
                #dont int sqrt for definite. losing accuracy
                MCC = round(MCC,3)
                
                #other performance measures                         
                TPR = round((TP/(TP+FN)),3)
                TNR = round((TN/(TN+FP)),3)
                PPV = round((TP/(TP+FP)),3)
                                
                #variables for average calculation 
                av_MCC_num += MCC
                av_MCC_count +=1
                MCC_list.append(MCC)
                print(x, ',' ,MCC, ',', TPR, ',',TNR, ',',PPV)
            except Exception as e:
                MCC = 'N/A'
                print(x, MCC, 'Error:', e)
                
                
        average_MCC = av_MCC_num/av_MCC_count 
        average_MCC = str(average_MCC)
        print('Simple MCC Average:',',',average_MCC) 
        return MCC_list