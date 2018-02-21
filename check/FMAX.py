# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 13:38:18 2018

@author: mcgerten
"""
import numpy



'''
F maximum
'''

def output(info, mode):
    ''' 
    Calculate the Fmax 
    
    Input:
    info --
    mode -- 
    '''
    # Intialize Variables
    fmax = 0.0
    fmax_threshold = 0.0
    PR = []
    RC = []
    
    # Run over all threshold values from 0 to 1, two signifigant digits
    for threshold in numpy.arange(0.00, 1.01, 0.01, float):
        
        threshold = numpy.around(threshold, decimals = 2)
        # Run PRRC on given threshold
        pr, rc = PRRC_average(info, threshold, mode)
        if pr is None:
            # No prediction above this threshold 
            break
        else:
            PR.append(pr)
            RC.append(rc)
            # Find the F-value for this particular threshold
            try:
                fval = f(pr, rc)
            except ZeroDivisionError:
                fval = None
                
        if fval is not None and fval >= fmax: ###########QUESTION##############
            fmax = fval
            fmax_threshold = threshold
    #Have found the Fmax at this point       
    return ([PR, RC, fmax, fmax_threshold])
    
    
def f(precision, recall):
    ''' Calculate F function '''
    
    try:
        f = (2*precision*recall)/(precision+recall)
    except ZeroDivisionError:
        f = None
    return f


def PRRC(info, threshold, protein):
    '''
    Calculate the PRRC of a single protein
    
    Input:
    info --
    threshold --
    protein --
    '''
    
    # Initalize Variables
    TP = 0.0      # True positive
    count = 0   # Count how many terms are above the threshold
    TT_length = len(info.true_terms[protein]) # Number of True terms
    
    if(threshold == 0):
        TP = TT_length
        count = info.obocount
    else:
        # For every term related to the protein
        for term in info.predicted_bench[protein]:
         # If it is above the threshold, increment the count
            if info.predicted_bench[protein][term][0] >= threshold:
                count += 1
                # If it is actually True, increment TP
                if info.predicted_bench[protein][term][1] :
                    TP += 1
    # Find PR: TP / (TP + FP)
    try:
        precision = TP / count 
    except ZeroDivisionError:
        precision = None
    # Find RC: TP (TP + FN)
    recall = TP/TT_length
    return (precision,recall)
        
        
def PRRC_average(info, threshold, mode):
    '''
    Calculate the overall PRRC of file
    
    Input:
    info --
    threshold --
    mode -- 
    '''
    
    # Initialize Variables
    PR = 0.0
    RC = 0.0
    info.count_above_threshold[threshold] = 0

    for protein in info.predicted_bench:
        pr,rc = PRRC(info, threshold, protein)
        if pr is not None:
            PR += pr
            info.count_above_threshold[threshold] += 1
        if rc is not None:
            RC += rc   
            
    if mode == 'partial':
        try:
            recall = RC/info.count_predictions_in_benchmark
        except ZeroDivisionError:
            recall = 0
            print("No protein in this predicted set became benchmarks\n")
            
    elif mode == 'full':
        try:
            recall = RC/info.count_true_terms
        except ZeroDivisionError:
            recall = 0
            print("No protein in this benchmark set\n")
            
    try:
        precision = PR/info.count_above_threshold[threshold]   
    except ZeroDivisionError:
        precision = None
        print("No prediction is made above the %.2f threshold\n" % threshold)
       
    return (precision, recall)     