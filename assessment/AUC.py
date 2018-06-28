import numpy
import helper
from collections import defaultdict


'''
AUC By GOterm 
NEEDS LOTS OF HELP -> USE NEW ONE 
'''

def output(info, ontology, Type, mode):
    ''' 

    
    Input:
    info   : Object
    mode   : Sting       {partial, full}
    
    Output:
    [0]    : 
    '''
    # Build the prediction and truth term lists from the protein based predictions
    TermList, TrueList = configureLists(info, mode)
    # Reduce list to only those who gained 10+ Positive Annoted Sequences (Benchmark has 10+)
    TermList, TrueList = reduceLists(info, mode, TermList, TrueList)
    # Intialize Variables
    Truth = set()
    # Build truth set
    for term in TrueList:    
        Truth.add(term)
    # For every Threshold
    for threshold in numpy.arange(0.00, 1.01, 0.01, float):
        threshold = numpy.around(threshold, decimals = 2)  
        Predicted = set()
        tempList = defaultdict(list)        
        # Build the predicted set
        for term in TermList:
        # If it is above the threshold
            for t in TermList[term]:
                protein    = t['protein']
                confidence = t['confidence']
                if confidence >= threshold:
                    tempList[term].append({'protein':protein, 'confidence':confidence})
        # Now have list of only predictions above threshold
        for term in tempList:
            Predicted.add(term)
            
        # Reset for each threshold
        CP = 0 # TP + FN
        # Effectively true positives
        
        CN = 0 # FP + TN
        # Effectively true negatives
            
        # Now have a predicted set
        if mode == "partial":
            # For all terms predicted
            for term in Predicted: # Use only those predicted at this threshold
                tp = TP(info, Truth, Predicted)
                fp = FP(info, Truth, Predicted)   
                if term in Truth:
                    CP += 1
                else: 
                    CN += 1
            
        else: #full
            for term in TermList: # Use all terms
                tp = TP(info, Truth, Predicted)
                fp = FP(info, Truth, Predicted)   
                if term in Truth:
                    CP += 1
                else: 
                    CN += 1
        
        tpr = tp / CP # Sensitivity
        fpr = fp / CN # Fall out
        
        terms = []
        areas = []
        

        for term in TermList:
            # Calculate Area
            area = numpy.trapz(TP)
            # Add to list
        
        # Want an object that per TERM, returns TP, FP, CP, CN, TPR, FPR, Area, Threshold
            
        #out = {} # {Go:term, Area}
        #outThreshold[threshold].append(out) 
    # Return all thresholds     
    # TPR, FPR, AREA are lists         
    return TPR, FPR, AREA, AUC, AUCThreshold
    
    
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
    
    # For every protein
    for protein in info.data:
        # Grab Term / confidence lists
        for t in info.data[protein]:
            # Seperate Term and confidence
            term       = t['term']
            confidence = t['confidence']
            #
            #print (t)
            #print (t['term'])
            #print (t['confidence'])
            #print (confidence)
            #print(term)
            # Append to a new dictionary Term : {protein, confidence}
            # termList is the prediction
            termList[term].append({'protein':protein, 'confidence':confidence})
            #print (termList[term])
        # If protein is in the benchmark   
        if info.true_terms[protein]:
            # For every term that is true, add it to the trueList
            for term in info.true_terms[protein]:
                trueList[term].append({'protein':protein, 'turth':True})
    # When done here, I will have a dictionary 
    # KEY: GO term,  Value: LIST [ PROTEIN, CONFIDENCE]
    # KEY: GO term,  Value: LIST [ PROTEIN, Truth]   
    return (termList, trueList)
    
def reduceLists(info, mode, TermList, TrueList):
    '''
    Reduce Lists to only those with 10+ proteins
    * what if term based benchmarks would be different??
    
    '''
        
    
    
    return TermList, TrueList