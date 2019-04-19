#%%
def difference_in_scores2(in_file,blout_queries_file,organism):
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
            
            for row in reader:
                score1 = float(row[2])
                score2 = float(row[4])
            
                prediction = str(row[1])   
                actual = str(row[6])     
    
                if phrases(x, prediction) == True and \
                phrases(my_dict[x],actual) == True:
                    TP +=1 
                    conf_matrix_type = 'TP'
                    
                elif phrases(x, prediction) == False and \
                phrases(my_dict[x],actual) == True:
                    FN +=1
                    conf_matrix_type = 'FN'
                    
                elif phrases(x,prediction) == True and \
                phrases(my_dict[x],actual) == False:
                    FP +=1
                    conf_matrix_type = 'FP'
                    
                elif phrases(x,prediction) == False and \
                phrases(my_dict[x],actual) == False:
                    TN +=1
                    conf_matrix_type = 'TN'
                
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
                    
                    
        incorrect_mu = sum(incorrect_diff)/len(incorrect_diff)
        incorrect_sigma = statistics.stdev(incorrect_diff)
        
        correct_mu = sum(correct_diff)/len(correct_diff)
        correct_sigma = statistics.stdev(correct_diff)
        
        print(incorrect_mu,incorrect_sigma)
        print(correct_mu,correct_sigma)
        
        #norm dist for correct
        x = np.linspace(correct_mu - 3*correct_sigma, correct_mu +3*correct_sigma,100)
        plt.plot(x,scipy.stats.norm.pdf(x,correct_mu,correct_sigma),linewidth=3.5, label = 'Correct',color = 'Red')
        #norm dist for incorrect
        z = np.linspace(incorrect_mu - 3*incorrect_sigma, incorrect_mu + 3*incorrect_sigma,100)
        plt.plot(z,scipy.stats.norm.pdf(z,incorrect_mu,incorrect_sigma),linewidth = 3.5,label = 'Incorrect',color = 'Blue')
        
        #formatting
        plt.xticks(fontsize = 24)
        plt.yticks(fontsize = 24)
        plt.show
    
