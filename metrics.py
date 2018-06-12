











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
    
def s(ru, mi):
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
    k = 2
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
    
def rumi_allProteins(data, threshold):
        '''
    Remaing Uncertainity, Misinformation   
    
    At a particular threshold, averaged over all proteins
    
    Input:
    info       : 
    k          : Integer    
    threshold  : Float      
    ontology   : String     {bpo, cco, mfo}
    Type       : String     {type1, type2}
    mode       : String     {partial, full}
    
    Output:
    [0]        : Float    Remaining Uncertainity
    [1]        : Float    Misinformation
    '''
    
    RU    = 0.0  
    MI    = 0.0
    count = 0
    # Set up mode
    if   data.mode == 'partial':
        count = info.count_predictions_in_benchmark
    elif data.mode == 'full':
        count = info.count_true_terms
            
    for protein in data.Prediction:
        r, m = rumi(info, threshold, protein, ontology, Type, mode)
        # Check both to ensure calculation worked
        if r is not None and m is not None:
            RU += r
            MI += m 
            
        
    try:
        remain  = RU / count
        misinfo = MI / count
        
    except ZeroDivisionError:
        remain  = None
        misinfo = None
        print("No prediction is made above the %.2f threshold\n" % threshold)
            
    return remain, misinfo 


def rumi_protein(data, threshold):
    '''
    Calculate the RUMI of a single protein
    
    Input:
    data       : Object
    threshold  : Float      {0.0 -> 1.0}
    
    Output:
    [0]        : Float      Remaining Uncertainity
    [1]        : Float      Misinformation
    '''
    
    data = open(info.path + "/SMIN/{}/{}/{}/{}/{}.txt".format(ontology, Type, mode, threshold, protein), 'w')
    
    T = set()
    P = set()
    # Find T
    T = info.true_terms[protein]
    # Find P (within threshold)
    for term in info.predicted_bench[protein]:
        # If it is above the threshold, add to the P set
        if info.predicted_bench[protein][term][0] >= threshold:
            try:                
                data.write('{}\t {}\t {}\n'.format(term, info.predicted_bench[protein][term][0], info.ic[term][1]))
            except KeyError:
                data.write('{}\t {}\t {}\n'.format(term, "ERROR", "ERROR"))            
            P.add(term) 
    data.close()
    # Calculate ru & mi     
    r = ru(info, T, P)
    m = mi(info, T, P)
    
    return r, m    