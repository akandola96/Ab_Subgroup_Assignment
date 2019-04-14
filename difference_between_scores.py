def difference_in_scores2(in_file,blout_queries_file,organism):
    import csv
    import statistics
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.mlab as mlab
    import scipy.stats
    import math
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
            training = 0 
            
            for row in reader:
                score1 = float(row[2])
                score2 = float(row[4])
            
                prediction = str(row[1])   #my assignment is in the 2nd col
                actual = str(row[6])     #blast results are in the 7th col 
    
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
    
#%%
#works as of 22/02    
def difference_in_scores(in_file,blout_queries_file,organism): 
    #e.g. h_sapiens_full_results.csv,h_sapiens_blout_queries, Homo sapiens 
    #only works on an individual organism basis. Doesnt check by subgroup
    #simply checks if assignment was correct or not and averages it out
    import csv
    import sys
    #sys.stdout = open(out_file,'a+')
    with open(in_file) as csv_in:
        reader = csv.reader(csv_in)
        count = 0 
        sum_diff = 0 
        count_diff = 0 
        wrong_sum_diff = 0 
        wrong_diff_count = 0 
        correct_assignment = []     #can use these to calculate the s.d.
        incorrect_assignment = []   #which can then be used in a t-test for signific
        
        my_dict = determine_subgroups(blout_queries_file,organism)
    
        for row in reader:
                                        
            assigned_subg = str(row[1])
            assigned_subg = assigned_subg[1:] #removes space before e.g. ' Mouse Heavy Chain 1'
            actual_subg = str(row[6])
            score1 = float(row[2])
            score2 = float(row[4])
            x = my_dict[assigned_subg]
            
            
            
            if phrases(x,actual_subg)==True: #does predictes
                diff = score1 - score2 #if so calc diff between scores
                sum_diff +=diff
                count_diff+=1
            elif phrases(x,actual_subg) == False: #was assignment incorrect?
                diff = score1-score2
                wrong_sum_diff += diff
                wrong_diff_count+=1
                
                
        av_diff_correct = sum_diff/count_diff
        av_diff_wrong = wrong_sum_diff/wrong_diff_count
        
        print('Organism:',organism)
        print('Average difference of correct assignment',av_diff_correct,'\n')
        print('Average difference of wrong assignment',av_diff_wrong,'\n')