'''
Class to calculate the Harmonic Mean Maximum (F-MAX)
'''

import numpy


def output (info, ontology, Type, mode):
    ''' 
    Main Method 
    
    Input:
    info     : Object
    ontology : String     {bpo, cco, mfo}
    Type     : String     {type1, type2}
    mode     : String     {partial, full}
    
    Output:
    [0]      : List[List[Float], List[Float], Float, Float]
    '''
    
    # Intialize Variables
    fmax           = 0.0
    fmax_threshold = -1.0
    PR             = []
    RC             = []
    F              = []
    # Run over all threshold values from 0 to 1, two signifigant digits
    for threshold in numpy.arange(0.00, 1.01, 0.01, float):
        threshold = numpy.around(threshold, decimals = 2)
        # Run PRRC on given threshold
        pr, rc = PRRC_average(info, threshold, ontology, Type, mode)
        if pr is None:
            # No prediction above this threshold 
            fval = None
            #break
        else:
            PR.append(pr)
            RC.append(rc)
            # Find the F-value for this particular threshold
            try:
                fval = f(pr, rc)
            except ZeroDivisionError:
                fval = None
        F.append(fval)        
        if fval is not None and fval >= fmax:
            fmax = fval
            fmax_threshold = threshold
    # Have found the Fmax at this point       
    return ([PR, RC, F, fmax, fmax_threshold])
    
    
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

 
def PRRC_average (info, threshold, ontology, Type, mode):
    '''
    Calculate the overall PRRC of file
    
    Input:
    info      : Object
    threshold : Float      {0.00 -> 1.00}
    ontology  : String     {bpo, cco, mfo}
    Type      : String     {type1, type2}
    mode      : String     {partial, full}
    
    Output:
    [0]       : Float     Precision Value
    [1]       : Float     Recall value
    '''
    
    # Initialize Variables
    PR = 0.0
    RC = 0.0
    info.count_above_threshold[threshold] = 0

    for protein in info.predicted_bench:
        pr, rc = PRRC(info, threshold, protein, ontology, Type, mode)
        if pr is not None:
            PR += pr
            info.count_above_threshold[threshold] += 1
        if rc is not None:
            RC += rc   
            
    if mode == 'partial':
        try:
            recall = RC / info.count_predictions_in_benchmark
        except ZeroDivisionError:
            recall = 0
            print("No protein in this predicted set became benchmarks\n")
            
    elif mode == 'full':
        try:
            recall = RC / info.count_true_terms
        except ZeroDivisionError:
            recall = 0
            print("No protein in this benchmark set\n")
    else:
        print("Invalid Mode")    
        
    try:
        precision = PR / info.count_above_threshold[threshold]   
    except ZeroDivisionError:
        precision = None
        print("No prediction is made above the %.2f threshold\n" % threshold)
       
    return (precision, recall)     
    
    
def PRRC (info, threshold, protein, ontology, Type, mode):
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
    POS = 0    # Count how many terms are above the threshold
    TRUE = len(info.true_terms[protein]) # Number of True terms
    
    
    if(threshold == 0):
        TP = TRUE
        POS = info.obocount
        
    else:
        data  = open(info.path + "/FMAX/{}/{}/{}/{}/{}.txt".format(ontology, Type, mode, threshold, protein), 'w')
        # For every term related to the protein
        for term in info.predicted_bench[protein]:
         # If it is above the threshold, increment the count
            if info.predicted_bench[protein][term][0] >= threshold:
                POS += 1
                # Write to file
                try:
                    data.write('{}\t {}\t {}\n'.format(term, POS, info.predicted_bench[protein][term][1])) 
                except KeyError:
                    data.write('{}\t {}\t {}\n'.format(term, POS, "KEY ERROR")) 
                
                # If it is actually True, increment TP
                if info.predicted_bench[protein][term][1] :
                    TP += 1
            data.write('Term, POS, TP, TRUE\n')                 
            data.write('{}\t {}\t {}\t {}\n'.format(term, POS, TP, TRUE))             
        data.close() 
    # Find PR: TP / (TP + FP)
    try:
        precision = TP / POS 
    except ZeroDivisionError:
        precision = None
    # Find RC: TP / (TP + FN)
    try:
        recall = TP / TRUE
    except ZeroDivisionError:
        recall = None
    
    return (precision,recall)