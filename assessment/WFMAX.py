'''
Weighted Harmonic Mean (WFMAX)
'''
import numpy
# Will use parts of Fmax just using IC values
import assessment.FMAX as F

def output(info, ontology, Type, mode):
    '''
    Calculate the weighted harmonic mean.
    
    Input:
    info     : Object
    ontology : String     {bpo, cco, mfo}
    Type     : String     {type1, type2}
    mode     : String     {partial, full}
    
    Output:
    [0]    : List[ List[Float], List[Float], Float, Float]
    '''
    
    # Intialize Variables
    wfmax           = 0.0
    wfmax_threshold = -1.0
    WPR             = []
    WRC             = []
    WF              = []
    # Run over all threshold values from 0 to 1, two signifigant digits
    for threshold in numpy.arange(0.00, 1.01, 0.01, float):
        threshold = numpy.around(threshold, decimals = 2)
        # Run PRRC on given threshold
        wpr, wrc = WPRRC_average(info, threshold, ontology, Type, mode)
        if wpr is None:
            # No prediction above this threshold 
            wf = None
            #break 
        else:
            WPR.append(wpr)
            WRC.append(wrc)
            # Find the F-value for this particular threshold
            try:
                wf = F.f(wpr, wrc)
            except ZeroDivisionError:
                wf = None
        WF.append(wf)   
        if wf is not None and wf >= wfmax: 
            wfmax = wf
            wfmax_threshold = threshold
    # Have found the WFmax at this point       
    return ([WPR, WRC, WF, wfmax, wfmax_threshold])
    
    
def WPRRC_average( info, threshold, ontology, Type, mode):
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
    WPR = 0.0
    WRC = 0.0
    info.count_above_threshold[threshold] = 0
    
    data  = open(info.path + "/WFMAX/{}/{}/{}/{}/summary.txt".format(ontology, Type, mode, threshold), 'w')
    for protein in info.predicted_bench:
        wpr, wrc = WPRRC(info, threshold, protein, ontology, Type, mode)
        data.write('{}\t {}\t {}\n'.format(protein, wpr, wrc)) 
        if wpr is not None:
            WPR += wpr
            info.count_above_threshold[threshold] += 1
        if wrc is not None:
            WRC += wrc   
            
    if mode == 'partial':
        try:
            recall = WRC / info.count_predictions_in_benchmark
        except ZeroDivisionError:
            recall = 0
            print("No protein in this predicted set became benchmarks\n")
            
    elif mode == 'full':
        try:
            recall = WRC / info.count_true_terms
        except ZeroDivisionError:
            recall = 0
            print("No protein in this benchmark set\n")
    else:
        print("Invalid Mode")         
    
    try:
        precision = WPR / info.count_above_threshold[threshold]   
    except ZeroDivisionError:
        precision = None 
        print("No prediction is made above the %.2f threshold\n" % threshold)
    # Write to file info to help debug
    data.write('{}\t {}\t {}\n'.format("WPR", "", WPR)) 
    data.write('{}\t {}\t {}\n'.format("WRC", "", WRC)) 
    data.write('{}\t {}\t {}\n'.format("Prediction in bench", "", info.count_predictions_in_benchmark)) 
    data.write('{}\t {}\t {}\n'.format("True Term count", "", info.count_true_terms)) 
    data.write('{}\t {}\t {}\n'.format("Count above threshold", "", info.count_above_threshold[threshold] )) 
    data.write('{}\t {}\t {}\n'.format("PR", "", precision)) 
    data.write('{}\t {}\t {}\n'.format("RC", "", recall)) 
    
    return (precision, recall) 
    
    
def WPRRC(info, threshold, protein, ontology, Type, mode):
    '''
    Calculate the PRRC of a single protein
    
    Input:
    info      : Object
    threshold : Float      {0.0 -> 1.0}
    protein   : 
    
    Ouput:
    [0]       : Float     Precision Value
    [1]       : Float     Recall value
    '''
    data  = open(info.path + "/WFMAX/{}/{}/{}/{}/{}.txt".format(ontology, Type, mode, threshold, protein), 'w')
    truth = open(info.path + "/WFMAX/{}/{}/{}/{}/{}_truth.txt".format(ontology, Type, mode, threshold, protein), 'w')
    # Initalize Variables
    total = 0.0        
    TP_total = 0.0      # True positive IC sum
    
    # Get the value of TT terms    
    TT_sum = 0.0
    for term in info.true_terms[protein]:
        try: 
            #print (info.ic[term][1] )
            TT_sum += info.ic[term][1] 
        # When prediction has newer terms than IC 
        except KeyError:
            data.write('{}\t {}\t {}\n'.format(term, "KEY ERROR", "TT_SUM"))  
            #print ("KEY ERROR")
        
        
    

    # For every term related to the protein
    for term in info.predicted_bench[protein]:
        # If the term is above the threshold
        if info.predicted_bench[protein][term][0] >= threshold:
            # Add IC value to total
            try:
                total += info.ic[term][1]
                data.write('{}\t {}\t {}\n'.format(term, info.predicted_bench[protein][term][0], info.ic[term][1])) 
            except KeyError:
                print("There was a key error on {}".format(term))
                data.write('{}\t {}\t {}\n'.format(term, "KEY ERROR", "KEY ERROR"))  
            # If it is actually True, add its IC to TP_total
            if info.predicted_bench[protein][term][1] :
                try:
                    TP_total += info.ic[term][1]
                    truth.write('{}\t {}\t {}\n'.format(term, info.predicted_bench[protein][term][0], info.ic[term][1])) 
                except KeyError:
                    truth.write('{}\t {}\t {}\n'.format(term, "KEY ERROR", "KEY ERROR"))  
    # Print out sums                
    data.write('{}\t {}\t {}\n'.format(protein, "Overall Prediction", total))  
    data.write('{}\t {}\t {}\n'.format(protein, "True Positive Total", TP_total))
    data.write('{}\t {}\t {}\n'.format(protein, "True Total", TT_sum))              
    data.close()

    truth.close()         
    # Find PR: TP / (TP + FP)
    try:
        precision = TP_total / total 
    except ZeroDivisionError:
        precision = None
    # Find RC: TP / (TP + FN) OR TP / TT (True Terms)
    try:
        recall = TP_total / TT_sum 
    except ZeroDivisionError:
        recall = None

    return (precision,recall)