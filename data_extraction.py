#%%
def xml_parser(input_file,query_organism,output_file_name):
    #this is a test
    """
    Summary:
    Extracts raw sequence data from abysis.
    
    Args:
    input_file = abysis data file (.xml)
    query_organism = organism name as found in abysis. Can also do genus (str)
    output_file_name = named output file (.csv)
    
    Desc: 
    Uses element tree to derive the XML tree and iterates through data by
    antibody chain. Extracts all non-pseudogene sequences beloning to query
    organism.
    
    Extracts residue, position number, and chothia numbering for each sequence. 
    If chothia numbering is unavailable, uses Kabat numbering instead.
    
    Each sequence extracted in the format of:
    > ID Acc
    H 1 H1
    L 2 H2
    V 3 H3
    * * *    
    """    
    import xml.etree.ElementTree as ET
    import sys
    sys.stdout = open(output_file_name, 'a')
    tree = ET.parse(input_file)                         # Builds the XML tree
    root = tree.getroot()                               # Gets the tree root
    query_organism = query_organism.lower()
    for antibody in root.iter('antibody'):
        for chain in antibody.iter('chain'):            # Iterate chains
            warning_present = chain.findall('warning')  # Check for warning 
            if warning_present:                          
                
                for warning in chain.iter('warning'):
                    wtype = warning.attrib['type']
                    if wtype != 'pseudogene':           # If not pseudogene        
                        accession = chain.find('accession').text 
                        name = chain.find('name').text
                        organism = chain.find('organism').text
                        if organism in query_organism:  # Write ID and acc                     
                            name = name.replace(',','') 
            
            
            # Will rewrite in later release to avoid code duplication
            elif not warning_present:                   # If no warning present
                accession = chain.find('accession').text 
                name = chain.find('name').text
                organism = chain.find('organism').text
                if organism in query_organism:  # Write ID and acc                     
                    name = name.replace(',','') 
                    
def xml_write(query_organism):
    sys.stdout.write('>' + ',' + name + ',' 
                     + accession + '\n')
    
    # Residues for-loop
    for residue in chain.iter('residue'):
        aa = residue.attrib['aa']
        pos = residue.attrib['pos']
        if 'chothia' in residue.attrib:    
            chothia = residue.attrib['chothia']
            sys.stdout.write(aa+','+pos+','
                             + chothia + '\n')
        elif 'kabat' in residue.attrib:      
            kabat = residue.attrib['kabat']                                
            sys.stdout.write(aa+','+ pos +',' 
                             + kabat + '\n') 
        else:
            continue                                    
    sys.stdout.write('*' + ',' + '*' +',' + '*' + '\n')
    
                
    
def count_abysis(in_file):
    """Counts number of records extracted from abYsis into csv file"""
    import csv   
    with open(in_file,'r') as csv_in:
        reader = csv.reader(csv_in)
        count = 0 
        for row in reader:
            if row[0] == '>':
                count +=1
        return(count)    
    
def original_make_fasta(input_file,output_file_name):
    """Extracts all sequences from abysis data"""
    import csv
    import sys
    csv_f = csv.reader(open(input_file,newline=''))
    sys.stdout=open(output_file_name,'a')
    for row in csv_f:
        if row[0] == '>':
            sys.stdout.write('\n' + row[0] + row[1] + '|' + row[2] + '\n')
            continue
        elif row[0] != '>':
            sys.stdout.write(row[0])
         
def make_fasta(input_file,version,output_file_name):
    """
    Summary:
    Converts extracted abYsis data into FASTA format. Two versions:
    placeholders or remove all.
    
    Args:
    input_file = extracted abysis data file (.csv)
    version = 'EA' or 'T6' (str)
        EA = Excludes all sequences for which Chothia/Kabat numbering begins
        past 1.
        T6 = Allows sequences with up to 6 missing N-terminal residues. X's 
        inserted instead.
    output_file_name = named output file (.fasta)
    
    Desc:
    Outer for loop scans through abYsis document and looks for '>', the start 
    of a sequence. Once this is a found a nested for loop iterates through the 
    first residue only. If this is found to begin on Chothia 1 then sequence is
    kept.

    Keeps sequences w/ insertions provided sequence does not begin on an 
    insertion.
    """
    import csv
    import sys
    csv_f = csv.reader(open(input_file))
    sys.stdout=open(output_file_name,'a')
    
    for row in csv_f:     
        if row[0] == '>':               # Finds start of a sequence
            count = 0 
            # Writes > ID | accession
            sys.stdout.write('\n' + row[0] + row[1] + '|' + row[2] + '\n') 
    
            # Nested for loop determines if sequence should be kept
            for row in csv_f:           # Checks within sequence          
               count +=1
               if row[0] == '*':        # If end of seq, exit loop
                   break
               if count < 2:            # Evaluate 1st residue only
                   numbering = str(row[2])
                   if '6A' in numbering:   # If first residue is an insertion
                       numbering = 7
                   elif '6A' not in numbering:
                       numbering = int(numbering[1:]) # Converts str(H1) to int(1)
                   if numbering == 1:       # If seq starts on 1
                       keep_marker = True 
                       sys.stdout.write(row[0])  # Write current residue
                       break
                   elif numbering != 1:     # If seq doesnt start on 1, evaluate
                       # Exclude all approach
                       if version == 'EA':                       
                           keep_marker = False
                           break
                       # Insert placeholder of up to 6X
                       elif version == 'T6':
                           if numbering < 8: # Allows up to a 6X insertion
                               missing_residues = numbering - 1
                               sys.stdout.write('X' * missing_residues) 
                               sys.stdout.write(row[0])
                               keep_marker = True
                               break
                           else:
                               keep_marker = False
                               break
        # If sequence is marked as correct then write residues
        elif row[0] == '*':                 
            continue
        else:
            if not keep_marker:                    # If incorrect
                continue
            elif keep_marker:                      # If correct
                sys.stdout.write(row[0]) # Writes remaining residues                  
                               
def remove_spaces(in_file):
    """
    Summary:
    Removes spaces in a fasta file.
    
    Args:
    in_file = input file (.fasta)
    
    Desc:
    Used to handle inconsistent labelling in IDs.
    
    """   
    with open(in_file,'r') as file:
        lines = file.readlines()
        
    # Replace spaces with empty string 
    lines = [line.replace(' ','') for line in lines]
    
    with open(in_file,'w') as file:
        file.writelines(lines)
   
def remove_short_sequences(infile,outfile):
    """
    Summary:
    Removes FASTA sequences of less than 21 residues
    
    Args:
    infile = input file (.fasta)
    outfile = output file (.fasta)
    """    
    from Bio import SeqIO
    import sys 
    sys.stdout = open(outfile,'a')
    for record in SeqIO.parse(infile,'fasta'):
        if len(record.seq) > 20: 
            print('>' + record.id)
            print(record.seq)    

def remove_seqs_missing_residues(in_file,version,out_file):
    """
    Summary:
    Removes sequences with unknown residues in their N-terminus
    
    Args:
    in_file = input file (.fasta)
    version = 'EA' or 'T6' (str)
        EA = removes all sequences with an 'X' present in first 21 residues 
        T6 = removes all sequences with more than 6 X residues in 
        first 21 residues
    out_file = output file (.fasta)
    
    Notes:
    Will be merged with remove_short_sequences function in later release
    Whilst earlier functions ensure non-zeroed sequences receive placeholder X
    residues, this function removes those sequences with placeholder residues 
    already present (i.e. fundamentally in the data, rather than added by this
    program)
    
    """
    from Bio import SeqIO
    import sys
    sys.stdout = open(out_file,'a')
    for record in SeqIO.parse(in_file,'fasta'):
        seq = str(record.seq)
        seq_21 = seq[:21]      # First 21 residues only
        if version == 'EA':
            if 'X' not in seq_21:
                print('>' + record.id)
                print(record.seq)  
        elif version == 'T6':
            if seq_21.count('X') <7:        # Maximum of 6
                print('>' + record.id)
                print(record.seq)
    
def seqkit_clean(in_file,os,out_file):  
    """
    Summary:
    Removes duplicates in a FASTA file by seq identified and seq record.
    
    Args:
    in_file = input file (.fasta)
    os = 'windows' or 'linux' (str)
        Need to differentiate between using 'type' or 'cat' on command line
    out_file = output file (.fasta)
    
    Desc:
    FASTA sequences of duplicate ID or Sequence are deleted.
    Uses subprocess.call rather than subprocess.Popen as .call is blocking.
    """

    if os == 'windows':
        import subprocess 
        cmd = ('type ' + in_file +' | seqkit rmdup -n | seqkit rmdup -s -o ' 
               + out_file)
        subprocess.Popen(cmd,shell=True)
    elif os == 'linux':
        cmd = ('cat ' + in_file +' | seqkit rmdup -n | seqkit rmdup -s -o ' 
               + out_file)
        subprocess.Popen(cmd)
       
def convert_seqkit(infile,outfile):
    """
    Summary:
    Converts ouput of seqkit into more readable format
    
    Args:
    infile = input file (.fasta)
    outfile = output file (.fasta)
    
    Desc:
    Output of SeqKit is a FASTA file with no spaces or new line chars. This 
    func simply reformats it to make it more readable.
    
    """
    outfile = 'queries.fasta'
    from Bio import SeqIO
    import sys 
    sys.stdout = open(outfile,'a')
    file = SeqIO.parse(infile,'fasta')
    for record in file:
        print('>' + record.id)  # Print auto inserts '\n'
        print(str(record.seq))
               
def count_num_queries(in_file):
    """
    Summary:
    Counts number of sequences in a fasta file
    
    Args:
    in_file = input file (.fasta)
    
    Desc:
    Used for figure making. Not essential. 
    """
    from Bio import SeqIO
    count = 0 
    for record in SeqIO.parse(in_file,'fasta'):
        if len(str(record.seq)) > 1:
            count+=1
    print(count)

def tidy_up():
    """Deletes redundant files"""
    import os
    os.remove('raw_queries_b.fasta')
    os.remove('raw_queries_c.fasta') 
    os.remove('raw_queries_d.fasta')
                        