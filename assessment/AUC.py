import numpy
import helper
from collections import defaultdict


'''
AUC By GOterm NEEDS A LOT OF HELP
'''

def output(info, ontology, Type, mode):
    ''' 

    
    Input:
    info   : Object
    mode   : Sting       {partial, full}
    
    Output:
    [0]    : 
    '''
    TermList, TrueList = configureLists(info, mode)
    # Intialize Variables
    
    Truth = set()
    
    
    
    # Build truth set
    for term in TrueList:    
        Truth.add(term)
    
    gained = {} # The set  of terms that gained 10 Positive Annoted Sequences
    outThreshold = {}
    for threshold in numpy.arange(0.00, 1.01, 0.01, float):
        threshold = numpy.around(threshold, decimals = 2)        
        # Reset for each threshold
        CP = 0 # TP + FN
        # Effectively true positives
        
        CN = 0 # FP + TN
        # Effectively true negatives
        
        for term in TermList:
            # For all terms
            if term in Truth:
                CP += 1
            else: 
                CN += 1
        
        terms = []
        areas = []
        
        Predicted = set()
        for term in TermList:
        # If it is above the threshold, add to the P set
            if TermList[term][1] >= threshold:
                Predicted.add(term) 
        # Now have A predicted set 
        tp = TP(info, Truth, Predicted)
        fp = FP(info, Truth, Predicted)                    
                            
                            
        
                
        tpr = tp / CP # Sensitivity
        fpr = fp / CN # Fall out
            

        for term in gained:
            # Calculate Area
            area = numpy.trapz(TP)
            # Add to list
            terms.append(term)
            areas.append(area)
            
        out = {} # {Go:term, Area}
        outThreshold[threshold].append(out) 
    # Return all thresholds     
    return outThreshold
    
    
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
    
def configureLists(info, mode):
    '''
    Want to first go through all predictions and make a list of Terms to Proteins
    '''
   
    termList = defaultdict(list)
    trueList = defaultdict(list)
    
    
    for protein in info.data:
        for term in protein:
            confidence = 0 # FIX THIS#################################################
            termList[term].append({protein, confidence})
        if (protein in info.true_terms[protein]):
            for term in info.termterms[protein]:
                trueList[term].append({protein, True})
    # When done here, I will have a dictionary 
    # KEY: GO term,  Value: LIST [ PROTEIN, CONFIDENCE]
    # KEY: GO term,  Value: LIST [ PROTEIN, Truth]   
    return (termList, trueList)