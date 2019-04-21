def make_random_fasta_file(number_seqs,out_file):
    """
    Summary:
    Generates random number of FASTA sequences of length 200.
    
    Args:
    number_seqs = number of sequences desired (int)
    out_file = named output file (.fasta)
    
    Desc:
    Was originally used to benchmark BLAST speed vs hsubgroup. 
    """
    
    import sys  
    import random 
    sys.stdout = open(out_file,'a+')
    number_seqs +=1
    
    
    aminos = ['A','C','D','E','F','G','H','I','K','L','M','N',
                          'P','Q','R','S','T','V','W','Y','X']   
    aminos = aminos*200
    x=1
    while x < number_seqs:
       
        seq = random.sample(aminos,75)
        seq = ''.join(seq)
        y = str(x)
        sys.stdout.write('>'+y+'\n')
        print(seq)
        
        x+=1