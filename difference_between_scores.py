#%%
def difference_in_scores(in_file,blout_queries_file,organism):
    """
    Summary:
    Calculates the difference between the primary and secondary scores of 
    correctly and incorrectly assigned sequences.
    
    Args:
    in_file = full results file for an organism (.csv)
    blout_queries_file = blout queries file for the organism (.csv)
    organism = organism (str)
    
    """
    
    import csv
    import statistics
    import matplotlib.pyplot as plt
    import numpy as np
    import scipy.stats
    with open(in_file) as csv_in:
        reader = csv.reader(csv_in)

        my_dict = determine_subgroups(blout_queries_file,organism)
        correct_diff = []
        incorrect_diff = []
                        
        for x in my_dict: # this comes from the determine_subgroups function
            
            TP = 0 #confusion matrix variables
            FP = 0 
            FN = 0 
            TN = 0 
            misc = 0
            alt_x = my_dict[x] + 'S'
            
            for row in reader:
                score1 = float(row[2])
                score2 = float(row[4])
            
                prediction = str(row[1])   
                actual = str(row[6])     
    
                # Assign entry to confusion matrix variable
                
                # If Heavy 1 is prediction and IGHV1 is actual TP + 1
                if phrases(x, prediction) == True and \
                phrases(my_dict[x],actual) == True:
                    TP+=1
                    conf_matrix_type = 'TP' 
                    
                # If Heavy 1 is prediction and IGHV1S1 is actual TP + 1
                elif phrases(x, prediction) == True and \
                alt_x in actual:
                    TP+=1
                    conf_matrix_type = 'TP'
                 
                    
                 # If Heavy 1 not prediction and IGHV1 is actual FN + 1   
                elif phrases(x, prediction) == False and \
                phrases(my_dict[x],actual) == True:
                    FN+=1
                    conf_matrix_type = 'FN'
                    
                # If Heavy 1 not prediction and IGHV1S1 is actual FN + 1
                elif phrases(x, prediction) == False and \
                alt_x in actual:
                    FN+=1
                    conf_matrix_type = 'FN'
                
        
                # if Heavy 1 is prediction and IGHV1 is not actual FP + 1
                elif phrases(x,prediction) == True and \
                phrases(my_dict[x],actual) == False:
                    FP+=1
                    conf_matrix_type = 'FP'
                    
                # If Heavy 1 is prediction and IGHV1S1 is not actual FP + 1    
                elif phrases(x,prediction) == True and \
                alt_x not in actual:
                    FP+=1
                    conf_matrix_type = 'FP'
                    
                    
                # If Heavy 1 is not prediction and IGHV1 is not actual TN + 1
                elif phrases(x,prediction) == False and \
                phrases(my_dict[x],actual) == False:
                    TN+=1
                    conf_matrix_type = 'TN'    
                # If Heavy 1 is not prediction and IGHV1S1 is not actual TN+1
                elif phrases(x,prediction) == False and \
                alt_x not in actual:
                    TN+=1
                    conf_matrix_type = 'TN'    
                    
                # Should always be 0
                else:
                    misc +=1   
                    
                if conf_matrix_type == 'TP':
                    diff = score1 - score2
                    correct_diff.append(diff)
                elif conf_matrix_type == 'FP':
                    diff = score1-score2
                    incorrect_diff.append(diff)
                elif conf_matrix_type == 'FN':
                    diff = score1 - score2
                    incorrect_diff.append(diff) 
                    
        # Calculations            
        incorrect_mu = sum(incorrect_diff)/len(incorrect_diff)
        incorrect_sigma = statistics.stdev(incorrect_diff)
        
        correct_mu = sum(correct_diff)/len(correct_diff)
        correct_sigma = statistics.stdev(correct_diff)
        
        print(incorrect_mu,incorrect_sigma)
        print(correct_mu,correct_sigma)
        
        # Normal Distributions
        
        # Norm dist for correct
        x = np.linspace(correct_mu - 3*correct_sigma, correct_mu +3*correct_sigma,100)
        plt.plot(x,scipy.stats.norm.pdf(x,correct_mu,correct_sigma),linewidth=3.5, label = 'Correct',color = 'Red')
        # Norm dist for incorrect
        z = np.linspace(incorrect_mu - 3*incorrect_sigma, incorrect_mu + 3*incorrect_sigma,100)
        plt.plot(z,scipy.stats.norm.pdf(z,incorrect_mu,incorrect_sigma),linewidth = 3.5,label = 'Incorrect',color = 'Blue')
        
        # Formatting
        plt.xticks(fontsize = 24)
        plt.yticks(fontsize = 24)
        plt.show
    
