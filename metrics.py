











def f (precision, recall):
    ''' 
    Calculate Harmonic Mean 
    
     Input:
     precision : Float
     recall    : Float
     
     Output:
     [0]       : Float
    '''
    
    try:
        fval = (2 * precision * recall) / (precision + recall)
    except ZeroDivisionError:
        fval = None
    return fval
    
    
def pr(TP, CountPredictedTrue):
    '''
    Per protein
    Find PR: TP / (TP + FP)
    '''
    # 
    try:
        precision = TP / CountPredictedTrue 
    except ZeroDivisionError:
        precision = None
    return precision    
        
def rc(TP, CountTrueTerms):
    '''
    Per protein
    Find RC: TP / (TP + FN)
    '''
    try:
        recall = TP / CountTrueTerms
    except ZeroDivisionError:
        recall = None
    return recall
    
def s(k, ru, mi):
    ''' 
    Semantic Distance 
    
    Input:
    k   : Integer      k-value (Normally 2)
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
    
    
def ru(IC, T, P):
    '''
    Calculate Remaining Uncertainity for a particular protein
    
    Input:
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
                total += IC[term]
            except KeyError:
                pass
            
    return total
    
    
def mi(IC, T, P):
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
                total += IC[term]
            except KeyError:
                pass
    return total  