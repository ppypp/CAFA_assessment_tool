# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 13:49:36 2018

@author: mcgerten
"""
import numpy
import FMAX

 # Will use parts of Fmax just using IC values
'''
Weighted F maximun
'''
def output(info, mode):
    '''
    Calculate WFmax
    
    Input:
    self :
    info :
    mode :
    
    Output:
    [0]  :
    '''
    
    # Intialize Variables
    wfmax = 0.0
    wfmax_threshold = 0.0
    WPR = []
    WRC = []
    
    # Run over all threshold values from 0 to 1, two signifigant digits
    for threshold in numpy.arange(0.00, 1.01, 0.01, float):
        #print(threshold)
        # print("\n")
        threshold = numpy.around(threshold, decimals = 2)
        # Run PRRC on given threshold
        wpr, wrc = WPRRC_average(info, threshold, mode)
        if wpr is None:
            # No prediction above this threshold 
            break
        else:
            WPR.append(wpr)
            WRC.append(wrc)
            # Find the F-value for this particular threshold
            try:
                wf = FMAX.f(wpr, wrc)
            except ZeroDivisionError:
                wf = None
                
        if wf is not None and wf >= wfmax: ###########QUESTION##############
            wfmax = wf
            wfmax_threshold = threshold
    #Have found the Fmax at this point       
    return ([WPR, WRC, wfmax, wfmax_threshold])
    
    
def WPRRC_average( info, threshold, mode):
    '''
    Calculate the overall PRRC of file
    
    Input:
    info      :
    threshold :
    mode      :
    
    Output:
    [0]       :
    '''
    
    # Initialize Variables
    WPR = 0.0
    WRC = 0.0
    info.count_above_threshold[threshold] = 0

    for protein in info.predicted_bench:
        wpr, wrc = WPRRC(info, threshold, protein)
        if wpr is not None:
            WPR += wpr
            info.count_above_threshold[threshold] += 1
        if wrc is not None:
            WRC += wrc   
            
    if mode == 'partial':
        try:
            recall = WRC/info.count_predictions_in_benchmark
        except ZeroDivisionError:
            recall = 0
            print("No protein in this predicted set became benchmarks\n")
            
    elif mode == 'full':
        try:
            recall = WRC/info.count_true_terms
        except ZeroDivisionError:
            recall = 0
            print("No protein in this benchmark set\n")
            
    try:
        precision = WPR/info.count_above_threshold[threshold]   
    except ZeroDivisionError:
        precision = None
        print("No prediction is made above the %.2f threshold\n" % threshold)
       
    return (precision, recall) 
    
    
def WPRRC(info, threshold, protein):
    '''
    Calculate the PRRC of a single protein
    
    Input:
    info      :
    threshold :
    protein   :
    
    Ouput:
    [0]       :
    '''
    
    # Initalize Variables
    total = 0.0        
    TP_total = 0.0      # True positive IC sum
    
    # Get the value of TT terms    
    TT_sum = 0.0
    for term in info.true_terms[protein]:
        try: 
            #print("TT" + term)
            #print("TT IC val" + str(info.ic[term][1]))
            #print("TT_sum" + str(TT_sum))
            TT_sum += info.ic[term][1] 
        except KeyError:
            TT_sum = TT_sum           
        
        

    # For every term related to the protein
    for term in info.predicted_bench[protein]:
        if info.predicted_bench[protein][term][0] >= threshold:
            #Add IC value to total
            try:
                total += info.ic[term][1]
            except KeyError:
                pass
                #print("KEY ERROR Total")
                #print(total)
            # If it is actually True, add its IC to TP_total
            if info.predicted_bench[protein][term][1] :
                try:
                    TP_total += info.ic[term][1]
                except KeyError:
                    pass
                    #print("KEY ERROR TP_total")
                    #print(TP_total)
                    
    # Find PR: TP / (TP + FP)
    try:
        precision = TP_total / total 
    except ZeroDivisionError:
        precision = None
    # Find RC: TP / (TP + FN) OR TP / TT (True Terms)
    try:
        recall = TP_total/TT_sum 
    except ZeroDivisionError:
        recall = None

    return (precision,recall)