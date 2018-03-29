# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 13:38:18 2018

@author: mcgerten
"""
import numpy
import helper


'''
AUC By GOterm apparently
'''

def output(info, mode):
    ''' 

    
    Input:
    info   : Object
    mode   : Sting       {partial, full}
    
    Output:
    [0]    : List[List[Float], List[Float], List[Float], Float]
    '''
    # Intialize Variables
    CP = 0
    CN = 0
    
    Truth     = set()
    Predicted = set()
    
    terms = []   
    areas = []
    
    gained = {} # The set  of terms that gained 10 Positive Annoted Sequences
    total = {} # The set of all terms \
    
    
    for term in total:
        if term is not gained:
            CN += 1
        if term is gained:
            CP += 1
            
            

    for term in gained:
        # Calculate Area
        area = numpy.trapz(TP)
        # Add to list
        terms.append(term)
        areas.append(area)
    return [terms, areas]
    
    
def TP(info, T, P):
    '''
    Calculate True Positives
    
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
                total += 1
                #total += info.ic[term][1]
            except KeyError:
                pass
            
    return total
    
    
def FP(info, T, P):
    '''
    Calculate False Positives
    
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
                total += 1
                #total += info.ic[term][1]
            except KeyError:
                pass
    return total            