#!/usr/bin/env python3


import re, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# MAIN:
if __name__ == "__main__":
    #path_logFile = "../3-1_deter/SRR4333334_margin.log"
    path_logFile = sys.argv[1]
    my_regex = re.compile("([0-9]+)(?:\siteration).+(?:likelihood:)\s+(-{0," +
                          "1}\d*\.{0,1}\d+).+(trial_[0-9].hmm)")
    data_frame  = pd.DataFrame(data=None)
    #print(data_frame) ; sys.exit()

    with open(path_logFile) as log_file:
        for line in log_file:
            matches = my_regex.search(line)
            
            if matches:
                iteration = int(matches.group(1))
                likeli = float(matches.group(2))
                trial = matches.group(3).split('.')[0]                
                data_frame.loc[iteration, trial] = likeli
                
    
    # Now we can display the plots:
    fig = plt.figure()
    axis = plt.subplot(111)
    plt.title("Likelihodd values along the different iterations")
    
    # A list of plots to add the different trials on the same plot
    plts = []
    col_df = data_frame.columns        
    
    for trial in col_df:
        plts.append(plt.plot(data_frame[trial]))
    #plt.plot(data_frame['trial_0'])
    plt.legend([plot[0] for plot in plts], col_df)
    plt.xlabel('Iterations')
    plt.ylabel("Likelihood (without 'e+')")
    
    plt.show()
    
    
    


