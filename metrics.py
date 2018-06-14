


def FMAX(data):
    return
def WFMAX(data):
    return
def SMIN(data):
    return
def NSMIN(data):   
    return





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
   
   
def prrc_allProteins (data, threshold):
    '''
    Calculate the overall PRRC of file
    
    Input:
    data      : Object
    threshold : Float      {0.00 -> 1.00}
    
    Output:
    [0]       : Float     Precision Value
    [1]       : Float     Recall value
    '''
    
    # Initialize Variables
    PR = 0.0
    RC = 0.0
    ProteinAboveThreshold = 0
    # Sum Values for all proteins
    for protein in data.Info:
        pr, rc = prrc_protein(data, threshold, protein)
        if pr is not None:
            PR += pr
            info.count_above_threshold[threshold] += 1
        if rc is not None:
            RC += rc   
            
    if mode == 'partial':
        try:
            recall = RC / data.ProteinInBenchmark
        except ZeroDivisionError:
            recall = 0
            print("No protein in this predicted set became benchmarks\n")
            
    elif mode == 'full':
        try:
            recall = RC / data.TrueTermsCount
        except ZeroDivisionError:
            recall = 0
            print("No protein in this benchmark set\n")
    else:
        print("Invalid Mode")    
        
    try:
        precision = PR / ProteinAboveThreshold  
    except ZeroDivisionError:
        precision = None
        print("No prediction is made above the %.2f threshold\n" % threshold)
       
    return (precision, recall)     
    
    
def prrc_protein (info, threshold, protein, ontology, Type, mode):
    '''
    Calculate the PRRC of a single protein
    
    Input:
    info       : Object
    threshold  : Float      {0.0 -> 1.0}
    protein    : 
    
    Output:
    [0]        : Float
    [1]        : Float
    '''
    
    # Initalize Variables
    TP = 0.0     # True positive
    count = 0    # Count how many terms are above the threshold
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
        data.write('Term, Count, TP')                 
        data.write('{}\t {}\t {}\n'.format(term, count, TP))             
        data.close() 
    # Find PR: TP / (TP + FP)
    try:
        precision = TP / count 
    except ZeroDivisionError:
        precision = None
    # Find RC: TP / (TP + FN)
    try:
        recall = TP / TT_length
    except ZeroDivisionError:
        recall = None
    
    return (precision,recall)
    
    
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
                vprint("{} has a KeyError".format(term), 4)
            
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
                vprint("{} has a KeyError".format(term), 4)
    return total  
    
def rumi_allProteins(data, threshold):
        '''
    Remaing Uncertainity, Misinformation   
    
    At a particular threshold, averaged over all proteins
    
    Input:
    data       : 
    threshold  : Float      
    
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
            
    for protein in data.Info:
        r, m = rumi(data, threshold, protein)
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


def rumi_protein(data, threshold, protein):
    '''
    Calculate the RUMI of a single protein
    
    Input:
    data       : Object
    threshold  : Float      {0.0 -> 1.0}
    
    Output:
    [0]        : Float      Remaining Uncertainity
    [1]        : Float      Misinformation
    '''
    

    T = set()
    P = set()
    # Find T
    T = data.ProteinTrueTerms[protein]
    # Find P (within threshold)
    for term in data.Info[protein]:
        # If it is above the threshold, add to the P set
        if data.Info[protein][term]['confidence']  >= threshold:
            #onfidence = data.Info[protein][term]['confidence'] 
            #ic = info.ic[term][1]           
            P.add(term) 
    # Calculate ru & mi     
    r = ru(info, T, P)
    m = mi(info, T, P)
    
    return r, m    