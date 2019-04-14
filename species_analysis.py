#%%
#works as of 22/02. Does what it needs to do #outputs to console      
def check_species_assignment(in_file,out_file):
    """Works out MCC on a species basis rather than subgroup"""
    #e.g. combined_full_results.csv
    import csv
    #import sys
    #sys.stdout = open(out_file,'a+')
    with open(in_file) as csv_in:
        reader = csv.reader(csv_in)
        
        organisms = ['Homo sapiens','Mus musculus','Oryctolagus cuniculus','Macaca mulatta','Gallus']
        
        for organism in organisms:
            csv_in.seek(0)
            TP = 0 #confusion matrix variables
            FP = 0 
            FN = 0 
            TN = 0 
            my_count = 0 
        
            for row in reader:
               actual_org = str(row[8]) #changed from 7 to 8 on 04/04
               predicted_org = str(row[1])
               
                   
               if organism in predicted_org and \
               organism in actual_org:
                   TP+=1
               elif organism not in predicted_org and \
               organism in actual_org:
                   FN+=1
               elif organism in predicted_org and \
               organism not in actual_org:
                   FP+=1
               elif organism not in predicted_org and \
               organism not in actual_org:
                   TN+=1
               else:
                   continue
               
               try:
                   org_MCC = int((TP * TN) - (FP * FN))/ \
                   int(((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)) ** 0.5)
               except:
                   continue
                   
            print(organism,'MCC:', org_MCC)
            print('TP:',TP,'FP:',FP,'FN:', FN ,'TN:', TN )
            print(TP+FP+FN+TN,'\n')
#%%
def determine_species_misassignment_master(in_file):
    """Works out which species were assigned as which"""
    query_organisms =['Homo sapiens', 'Mus musculus','Macaca mulatta','Oryctolagus cuniculus','Gallus']
    
    for organism in query_organisms:
        determine_species_misassignment(in_file,organism)

def determine_species_misassignment_FNs(in_file, organism):
    #works out which species gets misclassified as which most often (FNs only)
    import csv
    #import sys
    #sys.stdout = open(out_file,'a+')
    with open(in_file) as csv_in:
        reader = csv.reader(csv_in)
        macaca = 0 
        mouse = 0 
        oryctolagus = 0 
        homo = 0
        gallus = 0 
        
        
        for row in reader:
            assigned_org = str(row[1])
            actual_org = str(row[8]) #changed from 7 to 8 on 04/04
                        
            if organism not in assigned_org and organism in actual_org:
                if 'Macaca mulatta' in assigned_org:
                    macaca +=1
                elif 'Mus musculus' in assigned_org:
                    mouse+=1
                elif 'Oryctolagus cuniculus' in assigned_org:
                    oryctolagus+=1
                elif 'Homo sapiens' in assigned_org:
                    homo+=1
                elif 'Gallus' in assigned_org:
                    gallus+=1
                    
    
    print('Query organism:', organism)
    print('Macaca:',macaca,'Mouse:', mouse,'Oryctolagus:', oryctolagus, 'Homo:',homo,'Gallus: ',gallus)
    print('Sum:',macaca+mouse+oryctolagus+homo+gallus)
    print('\n')