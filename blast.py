#%%
def extract_ref_data(in_file,organism, region, out_file):
    """
    Summary:
    Extracts data from IMGT file.
    
    Args:
    in_file = IMGT data file (.fasta)
    organism = organism name as found in IMGT data (str)
    region = antibody chain region e.g. V, D, J (str)
    out_file = named output file (.fasta)
    
    Desc:
    Excludes sequences that are not functional genes.  
    """  
    from Bio import SeqIO
    import sys 
    sys.stdout = open (out_file, 'a')
    # Possible ways IMGT data denotes functional genee
    functional_gene = ['|F|', '|[F]|', '|(F)|']  
    region = region.upper() + '-REGION'     
    organism = organism.capitalize()        
    for seq_record in SeqIO.parse(in_file,'fasta'):
        if any(x in seq_record.description for x in functional_gene):
            if organism in seq_record.description:
                if region in seq_record.description:
                     print('>' + seq_record.description)    
                     print(seq_record.seq)                    
                 
def make_ref_blastdb(in_file,db_name):
    """
    Summary:
    Creates a reference nucleotide BLAST database from extracted IMGT 
    sequences.
    
    Args:
    in_file = input IMGT data file (.fasta)
    db_name = name of resultant database produced by this func (str)
    
    Desc:
    Makes a nucleotide database, also parses SeqIDs. NO masking.

    Will not work if duplicate sequences present in in_file. Use SeqKit to 
    clean if this happens.
    Runs a command line process via Python's subprocess module.   
    """
    import subprocess
    cmd = ('makeblastdb -in ' +in_file + ' -out '+db_name 
           +' -dbtype nucl -parse_seqids')
    subprocess.Popen(cmd)


def tBLASTn_full(queries,db_name,out_file):
    """
    Summary:
    Runs tBLASTn using IMGT reference data and abysis query sequences.
    
    Args:
    queries = query file (.fasta)
    db_name = name of BLAST database
    out_file = named output file (.xml)
    
    Desc:
    soft_masking false --> masking off  
    seg no --> masking off

    Runs a full tBLASTn search, all hits are included. Parsed later. 
    """    
    import subprocess 
    with open(out_file,'a') as out:
        cmd = ('tblastn -query ' + queries + ' -db ' + db_name 
               + ' -soft_masking false -seg no -outfmt 5')
        subprocess.Popen(cmd,stdout=out) 
        # Need to add a line that makes it wait till completion 

        
def count_xml_blast_records(in_file):
    """
    Summary:
    Counts the number of BLAST records in a BLAST output file
    
    Args:
    in_file = input file (.xml)
    
    Desc:
    Not essential.
    
    """
    from Bio.Blast import NCBIXML
    blast_record_count = 0 
    handle = open(in_file)    
    blast_records = NCBIXML.parse(handle)
    
    # Count records
    for blast_record in blast_records:
        blast_record_count +=1
    print(blast_record_count)     
    

#%%
def blast_output_xml2csv(in_file,hits,out_file):
    """
    Summary:
    Extracts top hit from XML formatted BLAST output. Converts to
    CSV format. 
    
    Args:
    in_file = raw blast xml output (.xml)
    alignments = number of alignments desired from each record (int)
    out_file = csv file containing a user defined number of alignments for each 
    blast record (csv)
    
    Desc:
    Records in XML file are ordered by E-value. Parsing the first record thus 
    gives the best hit by e-value. 
    BLAST alignment description contains subgroup assignment. 
    """   
    
    from Bio import SearchIO
    import sys
    sys.stdout=open(out_file,'a+')
    
    # q_result is list obj. Top hit stored q_result[0]
    hits -= 1
    blast_qresults = SearchIO.parse(in_file,'blast-xml')
    for blast_qresult in blast_qresults:  
        # Get query ID
        query_id = blast_qresult.id
        
        # Get defined hit + ID and description
        top_hit = blast_qresult[hits]
        top_hit_id = top_hit.id
        top_hit_desc = top_hit.description      
        
        # Write
        print(query_id,',', top_hit_id, top_hit_desc)
    
                    
#%%                    
                    
def convert_queries2csv(in_file,out_file):
    """
    Summary:
    Extracts sequences from FASTA file into csv file (first 21 residues only).
    
    Args:
    in_file = input file (.fasta)
    out_file = output file (.fasta)
    
    Desc:    
    Slices sequnce string to extract the first 21 residues only.  
    """
    from Bio import SeqIO
    import sys
    import csv
    sys.stdout = open(out_file,'a+')
    with open(out_file,'a+') as csv_output:
        writer = csv.writer(csv_output,lineterminator = '\n')
        for record in SeqIO.parse(in_file,'fasta'):
            seq = str(record.seq)
            seq = seq[:21]  # First 21 residues only
            print(seq)                    
                    
def join_queries2blout(blout_file,queries_file,out_file):
    """
    Summary:
    Attaches query sequence to its corresponding BLAST record.
    
    Args:    
    blout_file = Converted BLAST output file (.csv)
    queries_file = CSV file containing first 21 res of each query (.csv)
    out_file = Named output file (.csv)
    
    Desc:    
    Goes through rows in both blout_file and queries_file and writes both to 
    a new row in a new file (output file). Uses a zip function to do so.

    Output file formatting:
    ID, BLAST_record, Residues 
    MOPC|00001,RefSeq5115|IGHV1|Mus musculus,QVQLQVS...
    Anti-A|00002,RefSeq4123|IGKV1|Mus musculus,EVQLST...
    """   
    import csv
    import sys
    sys.stdout = open(out_file,'a+')
    
    with open(out_file,'a+') as csv_output, open(blout_file,'r') as blout, \
    open(queries_file,'r') as queries:
        
        # Read input rows
        writer = csv.writer(csv_output,lineterminator = '\n')
        blout_reader = csv.reader(blout)
        queries_reader = csv.reader(queries)
        
        # Write output row
        for blout_row, queries_row in zip(blout_reader,queries_reader):
            blout_row = ','.join(blout_row)
            queries_row = ''.join(queries_row)
            print(blout_row + ',' + queries_row)                    