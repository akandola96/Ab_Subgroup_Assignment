#%%
def extract_ref_data(in_file,organism, region, out_file):
    """Extracts data from IMGT file.
    
    in_file = IMGT data file (.fasta)
    organism = organism name as found in IMGT data
    region = antibody chain region e.g. V, D, J
    out_file = named output file (.fasta)
    
    Excludes sequences that are not functional genes.  
    """  
    from Bio import SeqIO
    import sys 
    sys.stdout = open (out_file, 'a')
    #possible ways IMGT data denotes functional gene
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
    """Creates a reference nucleotide BLAST database from extracted IMGT 
    sequences.
    
    Nucleotide database. Parses SeqIDs.

    Should be no duplicate sequences due to SeqKit steps. If duplicate 
    sequences are present, database construction will fail. 

    Runs a command line process.    
    """
    import subprocess
    cmd = ('makeblastdb -in ' +in_file + ' -out '+db_name 
           +' -dbtype nucl -parse_seqids')
    subprocess.Popen(cmd)


def tBLASTn_full(queries,db_name,out_file):
    """Runs tBLASTn using IMGT reference data and abysis query sequences.
    
    queries = query file (.fasta)
    db_name = name of BLAST database
    out_file = named output file (.xml)
    
    soft_masking false --> masking off  

    Runs a full tBLASTn search, all hits are included. 
    """    
    import subprocess 
    with open(out_file,'a') as out:
        cmd = ('tblastn -query ' + queries + ' -db ' + db_name 
               + ' -soft_masking false -outfmt 5')
        subprocess.Popen(cmd,stdout=out) 
        #need to add a line that makes it wait till completion 
        
def count_xml_blast_records(in_file):
    """Counts the number of BLAST records in a BLAST output file"""
    from Bio.Blast import NCBIXML
    blast_record_count = 0 
    handle = open(in_file)    
    blast_records = NCBIXML.parse(handle)
    
    #count records
    for blast_record in blast_records:
        blast_record_count +=1
    print(blast_record_count)     
    
 
def blast_output_xml2csv(in_file,alignments, out_file):
    """Extracts top hit (E-value) from XML formatted BLAST output. Converts to
    CSV format.
   
    in_file = raw blast xml output.
    alignments = number of alignments desired from each record.
    out_file = csv file containing a user defined number of alignments for each 
    blast record.

    Records in XML file are ordered by E-value. Parsing the first record thus 
    gives the best hit.
    """    
    from Bio.Blast import NCBIXML
    import sys
    handle = open(in_file)
    sys.stdout = open(out_file,'a+') 
    blast_records = NCBIXML.parse(handle)
    alignments+=1
    for record in blast_records:        
        if record.alignments:
            count = 0 
            for description in record.alignments:
                count +=1                
                if count < alignments:
                    print(record.query, ',', description.title)    
                    
def convert_queries2csv(in_file,out_file):
    """Extracts sequences from FASTA file into csv file.
    
    Splices sequnce string to extract the first 21 residues only.  
    """
    from Bio import SeqIO
    import sys
    import csv
    sys.stdout = open(out_file,'a+')
    with open(out_file,'a+') as csv_output:
        writer = csv.writer(csv_output,lineterminator = '\n')
        for record in SeqIO.parse(in_file,'fasta'):
            seq = str(record.seq)
            seq = seq[:21]  #first 21 residues only
            print(seq)                    
                    
def join_queries2blout(blout_file,queries_file,out_file):
    """Attaches query sequence to its corresponding BLAST record
    
    blout_file = BLAST output file (.csv)
    queries_file = Queries containing file (.csv)
    out_file = Named output file (.csv)
    
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
        
        #read input rows
        writer = csv.writer(csv_output,lineterminator = '\n')
        blout_reader = csv.reader(blout)
        queries_reader = csv.reader(queries)
        
        #write output row
        for blout_row, queries_row in zip(blout_reader,queries_reader):
            blout_row = ','.join(blout_row)
            queries_row = ''.join(queries_row)
            print(blout_row + ',' + queries_row)                    