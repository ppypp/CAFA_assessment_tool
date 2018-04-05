# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 13:38:42 2018

@author: mcgerten
"""

'''
S minimum
'''

import numpy
import helper


def output(info, ontology, Type, mode):
    '''
    Calculate Smin
    
    k = 2 by convention
    
    Input:
    info : Object defined in check.py
    mode : String  (Full, Partial)
    
    Output:
    [0]  : List[RU, MI, smin, smin_threshold]
    
    '''
    
    k = 2.0
    # Intialize Variables
    smin = float("inf")
    smin_threshold = 0.0
    RU   = []
    MI   = []
    S    = []
    # Run over all threshold values from 0 to 1, two signifigant digits
    for threshold in numpy.arange(0.00, 1.01, 0.01, float):
        
        threshold = numpy.around(threshold, decimals = 2)
        # Run S on given threshold
        ru, mi = rumi_average(info, k, threshold, ontology, Type, mode)
        #print("Protein Totals")
        #print(ru)
        #print(mi)
        RU.append(ru)
        MI.append(mi)
        # Find the S-value for this particular threshold
        sval = s(k, ru, mi)
        #print(sval)
        S.append(sval)
        if (sval is not None and sval <= smin): 
            smin = sval
            smin_threshold = threshold
                        
            #print("New Low")
    # Have found the Smin at this point       
    return ([RU, MI, S, smin, smin_threshold])


def ru(info, T, P):
    '''
    Calculate Remaining Uncertainity for a particular protein
    
    Input:averaged
    T   : Set [ Truth      ]
    P   : Set [ Prediction ]
    
    Output:
    [0] : Float 
    '''
    
    # Sum of IC values of terms meeting criteria   
    total = 0.0
    # Sum the IC values of every element in T not in P
    for term in T:
        if term not in P:
            try:
                total += info.ic[term][1]
            except KeyError:
                pass
            
    return total
    
    
def mi(info, T, P):
    '''
    Calculate Misinformation for a particular protein
    
    Input:    
    T   : Set [ Truth      ]
    P   : Set [ Prediction ]
    
    Output:
    [0] : Float
    '''

    # Sum of IC values of terms meeting criteria   
    total = 0.0
    # Sum the IC values of every element in P not in T
    for term in P:
        if term not in T:
            try:
                total += info.ic[term][1]
            except KeyError:
                pass
    return total            
    

def s(k, ru, mi):
    ''' 
    Semantic Distance 
    
    Input:
    k   : Integer      k-value
    ru  : Float        remaining uncertainity
    mi  : Float        misinformation
    
    Output:
    [0] : Float        s-value
    '''

    s = 0.0
    
    if(ru is None or mi is None):
        s = None
    else:
        s = (ru**k + mi**k)**(1 / k)
    return s
    

def rumi_average(info, k, threshold, ontology, Type, mode):
    '''
    Remaing Uncertainity, Misinformation   
    
    At a particular threshold, averaged over all proteins
    
    Input:
    info       : 
    k          : Integer  (Typically 2)
    threshold  : Float    A cutoff value
    mode       : String   Mode to run
    
    Output:
    [0]        : Float    Remaining Uncertainity
    [1]        : Float    Misinformation
    '''
    
    RU = 0.0  
    MI = 0.0
#    # Ne 
#    countPredicted = 0 
#    # N    
#    countTotal     = 0
    count          = 0
    for protein in info.predicted_bench:
        T = set()
        P = set()
        # Find T
        T = info.true_terms[protein]
#        for term in info.true_terms[protein]:
#            countTotal += 1
        # Find P (within threshold)
        for term in info.predicted_bench[protein]:
            # If it is above the threshold, add to the P set
            if info.predicted_bench[protein][term][0] >= threshold:
                P.add(term) 
#                countPredicted += 1
        # Calculate ru & mi     
        r = ru(info, T, P)
        #print("RU sum per protein")
        #print(r)
        m = mi(info, T, P)
        #print("MI sum per protein")
        #print(m)
        # Check both to ensure calculation worked
        if r is not None and m is not None:
            RU += r
            MI += m 
            #info.count_above_threshold[threshold] += 1
#        if (mode == "full"):
#            # Count is the number of benchmark proteins
#            count = countTotal
#        else: # "partial"
#            # Count is the number of predicted proteins
#            count = countPredicted
        if mode == 'partial':
            count = info.count_predictions_in_benchmark
        elif mode == 'full':
            count = info.count_true_terms
    try:
        remain  = RU / count
        misinfo = MI / count
        
    except ZeroDivisionError:
        remain  = None
        misinfo = None
        print("No prediction is made above the %.2f threshold\n" % threshold)
            
    return remain, misinfo        